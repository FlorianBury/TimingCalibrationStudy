import os
import sys
import copy
import yaml
import math
import glob
import datetime
import itertools
import argparse
import subprocess
import numpy as np
import ROOT
from pprint import pprint 

from .environment import getEnv

default_params = {'N'                  : 100,
                  'threshold'          : 5800, 
                  'thresholdsmearing'  : 0., 
                  'tofsmearing'        : 0., 
                  'mode'               : 'scan',
                  'offset'             : -1.,
               }

CMSSW_DIR = os.path.join(getEnv()['paths']['cmssw'])
SETUP_CMSSW = (
               "module --force purge;"
               "module load cp3;"
               "module load grid/grid_environment_sl7;"
               "module load crab/crab3;"
               "module load cms/cmssw;"
               "module load slurm/slurm_utils;"
               "eval `scramv1 runtime -sh`"
               )
HARVESTER_SCRIPT = 'Harvester_cfg.py'
SLURM_OUTPUT_DIR = getEnv()['paths']['production']

class Scan:
    def setParameterDict(self,params):
        self.params = {}
        for pName,pVal in params.items():
            if not isinstance(pVal,(list,tuple,np.ndarray)):
                pVal = [pVal]
            self.params[pName] = pVal
        print ('Parameters for scan :')
        pprint (self.params)
        self.parameterCombinations()

    def parameterCombinations(self):
        self.inputParamsNames = list(self.params.keys())
        self.inputParams = []
        for prod in itertools.product(*self.params.values()):
            inputParam = []
            for p in prod:
                if isinstance(p,float):
                    p = round(p,10) # Truncation of 0.000..00x problems
                inputParam.append(str(p))
            self.inputParams.append(inputParam)
        print (f'Generated {len(self.inputParams)} combinations')

    def getParameters(self):
        return [{pName:pVal for pName,pVal in zip(self.inputParamsNames,params)} for params in self.inputParams]

    @staticmethod
    def run_command(command,return_output=False,**kwargs):
        process = subprocess.Popen(command,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,**kwargs)
        # Poll process for new output until finished #
        output = []
        while True:
            try:
                nextline = process.stdout.readline()
            except UnicodeDecodeError:
                continue
            if nextline == '' and process.poll() is not None:
                break
            print(nextline.strip())
            if return_output:
                output.append(nextline)
        process.communicate()
        exitCode = process.returncode
        if return_output:
            return exitCode,output
        else:
            return exitCode

    @staticmethod
    def in_virtualenv():
        # Get base/real prefix, or sys.prefix #
        # Check if matches with sys.prefix #
        if hasattr(sys, "base_prefix"):
            return sys.base_prefix != sys.prefix
        elif hasattr(sys, "real_prefix"): # old versions 
            return sys.real_prefix != sys.prefix
        else:
            return False

    @classmethod
    def format_command(cls,cmd):
        full_cmd = ""
        # Check if virtual env and disable it in case #
        #if cls.in_virtualenv():
        #    path_virtual_env = os.environ.get('VIRTUAL_ENV')
        #    path_activate = os.path.join(path_virtual_env,'bin','activate')
        #    if not os.path.exists(path_activate):
        #        raise RuntimeError(f'Weird, could not find {path_activate}')
        #    full_cmd += f"source {path_activate}; deactivate; export PYTHONPATH=''; "
        #    # Need to first source the activation to be able to deactivate
        #    # Also need to clear python path to remove the LCG 
        #    # So far not working though ...
        #    raise RuntimeError('Recommended for now to not start from virtual env ...')
        cwd = os.getcwd()
        full_cmd += f"cd {CMSSW_DIR} ; {SETUP_CMSSW} ; cd {cwd}"
        if isinstance(cmd,str):
            full_cmd += f" && {cmd}"
        elif isinstance(cmd,list):
            full_cmd += f" && {' '.join(cmd)}"
        else:
            return ValueError
        return full_cmd

    @classmethod
    def run(cls,args):
        args_dict = copy.deepcopy(default_params)
        args_dict.update(args)
        if 'script' not in args_dict.keys():
            raise RuntimeError('Missing `script` entry')
        else:
            script = args_dict['script']
            del args_dict['script']
        args = [f"{k}={v}" for k,v in args_dict.items()]

        print ('Starting the DQM file production')
        print ('Arguments : '+' '.join(args))
        dqm_cmd = cls.format_command(['cmsRun',os.path.join(CMSSW_DIR,script)] + args)

        rc,output = cls.run_command(dqm_cmd,return_output=True,shell=True)
        dqm_file = None
        for line in output:
            print (line)
            if 'BXHist' in line and '.root' in line:
                for l in line.split():
                    if '.root' in l:
                        dqm_file = l
                if dqm_file is not None:
                    break
        print (f'... exit code : {rc}')
        if rc != 0:
            raise RuntimeError("Failed to produce the DQM root file")
        if dqm_file is None or not '.root' in dqm_file:
            raise RuntimeError(f"Wrong output root file : {dqm_file}")

        print (f"DQM root file created as {dqm_file}")
        if not os.path.exists(dqm_file):
            raise RuntimeError("DQM root file not present")
            
        print ('Starting harvesting')
        harvest_cmd = cls.format_command(['cmsRun',os.path.join(CMSSW_DIR,HARVESTER_SCRIPT),f'input={dqm_file}'])
        rc = cls.run_command(harvest_cmd,shell=True)
        print (f'... exit code : {rc}')
        if rc != 0:
            raise RuntimeError("Harvesting failed")

        print ('Starting renaming')
        hist_file = dqm_file.replace('raw','harvested')
        dirname = os.path.dirname(dqm_file)
        rename_cmd = ['mv',os.path.join(dirname,'DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root'),hist_file]
        rc = cls.run_command(rename_cmd)
        if rc != 0:
            raise RuntimeError("Could not rename the harvested file")
        print (f'... exit code : {rc}')
        print (f'Renamed to {hist_file}')

        print ('Starting cleaning')
        clean_cmd = ['rm', dqm_file]
        rc = cls.run_command(clean_cmd)
        if rc != 0:
            raise RuntimeError("Could note clean intermediate root file")
        print (f'... exit code : {rc}')

        F = ROOT.TFile(hist_file,"UPDATE")
        for name,arg in args_dict.items():
            p = ROOT.TNamed(name,str(arg))
            p.Write()
        F.Close()

    def findMissingJobs(self,path):
        with open(os.path.join(path,'infiles.yml'),'rb') as handle:
            infiles = yaml.safe_load(handle,)
        N = len(infiles)
        print (f"Looking over directory {path}")

        import enlighten
        rootfiles = list(glob.glob(os.path.join(path,'output','*.root')))
        pbar = enlighten.Counter(total=len(rootfiles), desc='Progress', unit='Files')

        for rootfile in rootfiles:
            pbar.update()
            F = ROOT.TFile(rootfile,"READ")
            ROOT.SetOwnership(F,False)
            params = {}
            for key in F.GetListOfKeys():
                if key.GetClassName() == "TNamed":
                    params[key.GetName()] = F.Get(key.GetName()).GetTitle()
            if params not in infiles:
                print ('Warning : parameters not in infiles',params)
            else:
                infiles.remove(params) 
            F.Close()

        print (f"All jobs     = {N}")
        print (f"Missing jobs = {len(infiles)}")
        if len(infiles) > 0:
            print ("Do not forget to send the jobs with --submit")
            for infile in infiles:
                print (infile)

            self.inputParamsNames = list(infiles[0].keys())
            self.inputParams = [[infile[name] for name in self.inputParamsNames] for infile in infiles]
                    
    def submit_on_slurm(self,name,debug=False):
        # Slurm configuration
        from CP3SlurmUtils.Configuration import Configuration
        from CP3SlurmUtils.SubmitWorker import SubmitWorker
        from CP3SlurmUtils.Exceptions import CP3SlurmUtilsException

        config = Configuration()
        config.sbatch_partition = 'cp3'
        config.sbatch_qos = 'cp3'
        config.sbatch_time = '0-02:00:00'
        config.sbatch_memPerCPU= '3000'
        config.inputSandboxContent = ['*py']
        config.sbatch_additionalOptions = ["--export=All"]
        config.useJobArray = True
        config.inputParamsNames = []
        config.inputParams = []

        base_payload = "BXRun --run " + " ".join([n+"={"+n+"}" for n in self.inputParamsNames])
    
        config.payload = "${taskcmd}"
    
        timestamp = datetime.datetime.now().strftime(name+'_%Y-%m-%d_%H-%M-%S')
    
        slurm_config = copy.deepcopy(config)
        slurm_working_dir = os.path.join(SLURM_OUTPUT_DIR,timestamp)
    
        slurm_config.batchScriptsDir = os.path.join(slurm_working_dir, 'scripts')
        slurm_config.inputSandboxDir = slurm_config.batchScriptsDir
        slurm_config.stageoutDir = os.path.join(slurm_working_dir, 'output')
        slurm_config.stageoutLogsDir = os.path.join(slurm_working_dir, 'logs')
        slurm_config.stageoutFiles = ["*.root"]

        all_params = []
        for params in self.inputParams:
            paramSet = {name:str(params[i]) for i,name in enumerate(self.inputParamsNames)}
            paramSet.update({k:v for k,v in default_params.items() if k not in paramSet.keys()})
            all_params.append(paramSet)

        slurm_config.inputParamsNames = ['taskcmd']
        maxArr = 3000
        Njobs = int(math.ceil(float(len(self.inputParams))/maxArr))
        for job in range(Njobs):
            print ('Submitting batch of jobs %d/%d'%(job,Njobs))
            slurm_config.inputParams = []
            for inputParam in self.inputParams[job*maxArr:(job+1)*maxArr]:
                dParam = {n:p for n,p in zip(self.inputParamsNames,inputParam)}
                taskcmd = base_payload.format(**dParam)
                slurm_config.inputParams.append([taskcmd])
            # Submit job!
            print("Submitting job...")
            if not debug:
                submitWorker = SubmitWorker(slurm_config, submit=True, yes=True, debug=False, quiet=False)
                submitWorker()
                print("Done")
            else:
                print(f'Submitting {len(slurm_config.inputParams)} jobs')
                print(f'Payload : {slurm_config.payload}')
                print(f'Param names : {slurm_config.inputParamsNamea}')
                print ('Parameters :')
                for inputParam in slurm_config.inputParams:
                    print ("\t",inputParam)
                print('... don\'t worry, jobs not sent')

        # Save params #
        if not debug:
            all_params = []
            for params in self.inputParams:
                paramSet = {name:str(params[i]) for i,name in enumerate(self.inputParamsNames)}
                paramSet.update({k:str(v) for k,v in default_params.items() if k not in paramSet.keys()})
                if 'script' in paramSet:
                    del paramSet['script']
                all_params.append(paramSet)

            with open(os.path.join(slurm_working_dir,'infiles.yml'),'w') as handle:
                yaml.dump(all_params,handle)

        

def main():
    parser = argparse.ArgumentParser(description='Timing calibration setup job submission')
    parser.add_argument('--yaml',action='store',required=False,type=str,default=None,
                        help='Config to run several modes')
    parser.add_argument('--submit',action='store',required=False,type=str,default=None,
                        help='Name for job submission')
    parser.add_argument('--resubmit',action='store',required=False,type=str,default=None,
                        help='Slurm output dir to check for previous results and only resubmit the missing ones')
    parser.add_argument('--run',nargs='*',required=False,type=str,default=None,
                        help='Run parameters (with `=` between name of the arg and the value)')
    parser.add_argument('--debug',action='store_true',required=False,default=False,
                        help='Debug for job submission')

    args = parser.parse_args()


    instance = Scan()
    
    # Config #
    if args.yaml is not None:
        with open(args.yaml,'r') as handle:
            config = yaml.load(handle,Loader=yaml.FullLoader)
        instance.setParameterDict(config)
        params = instance.getParameters()
    else:
        params = [{}]

    # Serial run #
    if args.run is not None:
        runParams = {arg.split("=")[0]:arg.split("=")[1] for arg in args.run}
        if len(runParams) > 0 and len(params) > 0:
            print (f'[WARNING] Will overwrite following parameters from the config {runParams}')
        for param in params:
            param.update(runParams)
            Scan.run(param)
        sys.exit()

    # Resubmit #
    if args.resubmit is not None:
        instance.findMissingJobs(args.resubmit)
        if args.submit is not None:
            instance.submit_on_slurm(args.submit,args.debug)
        sys.exit()
    # Submit #
    if args.submit is not None:
        instance.submit_on_slurm(args.submit,args.debug)
        sys.exit()

if __name__ == "__main__":
    main()
