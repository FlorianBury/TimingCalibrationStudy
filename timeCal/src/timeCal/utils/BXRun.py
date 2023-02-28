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
                  'subdet'             : 'ALL',
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
VENV = getEnv()['paths']['venv']

MAX_ARR = 3000

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
            nextline = nextline.strip()
            if len(nextline) == 0:
                continue
            print (nextline)
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
    def get_env(cls):
        env = os.environ.copy()
        if cls.in_virtualenv():
            # If in virtual environment, gotta forge a copy of the environment, where we:
            # Delete the VIRTUAL_ENV variable.
            del(env['VIRTUAL_ENV'])

            # Delete the venv path from the front of my PATH variable.
            orig_path = env['PATH']
            virtual_env_prefix = sys.prefix + '/bin:'
            env['PATH'] = orig_path.replace(virtual_env_prefix, '')

        return env

    @classmethod
    def format_command(cls,cmd):
        cwd = os.getcwd()
        full_cmd = f"cd {CMSSW_DIR} ; {SETUP_CMSSW} ; cd {cwd}"
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

        rc,output = cls.run_command(dqm_cmd,return_output=True,shell=True,env=cls.get_env())
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
        rc = cls.run_command(harvest_cmd,shell=True,env=cls.get_env())
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

    @staticmethod
    def getFileParams(f):
        F = ROOT.TFile(f,"READ")
        params = {}
        for key in F.GetListOfKeys():
            if key.GetClassName() == "TNamed":
                params[key.GetName()] = F.Get(key.GetName()).GetTitle()
        F.Close()
        return params


    def findMissingJobs(self,path):
        with open(os.path.join(path,'infiles.yml'),'rb') as handle:
            infiles = yaml.safe_load(handle,)
        N = len(infiles)
        print (f"Looking over directory {path}")

        import enlighten
        rootfiles = list(glob.glob(os.path.join(path,'output','*.root')))
        pbar = enlighten.Counter(total=len(rootfiles), desc='Progress', unit='Files')

        idx_missing = [True] * N
        for rootfile in rootfiles:
            params = self.getFileParams(rootfile)
            if params not in infiles:
                raise RuntimeError(f'Cannot find {params} in infiles')
            idx_missing[infiles.index(params)] = False
            pbar.update()

        idx_resubmit = np.array([i for i in range(N) if idx_missing[i]])
        print (f"All jobs     = {N}")
        print (f"Missing jobs = {len(idx_resubmit)}")

        if len(idx_resubmit) > 0:
            slurm_scripts = glob.glob(os.path.join(path,'scripts','slurmSubmission*.sh'))

            print ('Missing parameters :')
            for i in idx_resubmit:
                print (f'{i}/{len(idx_resubmit)}',infiles[i])

            print ('Use following command(s) to resubmit missing jobs')

            # Run over all script bash script #
            for j in range(len(slurm_scripts)):
                # make and check name #
                script = os.path.join(path,'scripts',f'slurmSubmission_{j}.sh')
                if not os.path.exists(script):
                    raise RuntimeError(f'Cannot find {script}')

                # get all indices #
                idx_min = j*MAX_ARR
                idx_max = (j+1)*MAX_ARR
                idx = idx_resubmit[np.logical_and(idx_resubmit>=idx_min,idx_resubmit<idx_max)] - idx_min

                # Command #
                if len(idx) > 0:
                    cmd = f"sbatch --array={','.join([str(i+1) for i in idx])} {script}"
                    print (cmd)
                    print ()
        else:
            print ('All jobs have succeeded')

    def submit_on_slurm(self,name,debug=False,slurm_params={}):
        # Slurm configuration
        from CP3SlurmUtils.Configuration import Configuration
        from CP3SlurmUtils.SubmitWorker import SubmitWorker
        from CP3SlurmUtils.Exceptions import CP3SlurmUtilsException

        config = Configuration()
        config.sbatch_partition = 'cp3'
        config.sbatch_qos = 'cp3'
        config.sbatch_time = '0-04:00:00'
        config.sbatch_memPerCPU= '3000'
        config.inputSandboxContent = []
        config.sbatch_additionalOptions = ["--export=None"]
        config.useJobArray = True
        config.inputParamsNames = []
        config.inputParams = []

        # override #
        for key,val in slurm_params.items():
            setattr(config,f'sbatch_{key}',str(val))

        # base payload #
        base_payload = "BXRun --run " + " ".join([n+"={"+n+"}" for n in self.inputParamsNames])
        # TODO : figure out why when submitting on slurm, if the virtual env is passed through it then the subprocess with cmSRun fails

        config.payload = "${BX_cmd}"

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

        slurm_config.inputParamsNames = ['BX_cmd']
        Njobs = int(math.ceil(float(len(self.inputParams))/MAX_ARR))
        slurm_ids = []
        for job in range(Njobs):
            print ('Submitting batch of jobs %d/%d'%(job,Njobs))

            # Make name of bash job #
            slurm_config.batchScriptsFilename = f"slurmSubmission_{job}.sh"
            slurm_script = os.path.join(slurm_config.batchScriptsDir,slurm_config.batchScriptsFilename)

            # Change params #
            slurm_config.inputParams = []
            for inputParam in self.inputParams[job*MAX_ARR:(job+1)*MAX_ARR]:
                dParam = {n:p for n,p in zip(self.inputParamsNames,inputParam)}
                slurm_config.inputParams.append([base_payload.format(**dParam)])

            # Submit job!
            print("Submitting job...")
            if not debug:
                submitWorker = SubmitWorker(slurm_config, submit=False, yes=True, debug=False, quiet=True)
                submitWorker()
                if not os.path.exists(slurm_script):
                    raise RuntimeError(f'File {slurm_script} should have been created')

                content = []
                with open(slurm_script,'r') as handle:
                    for line in handle:
                        content.append(line)
                        if "Starting read of input parameters" in line:
                            content.append(f"\tsource {VENV}\n")

                with open(slurm_script,'w') as handle:
                    for line in content:
                        handle.write(line)
                rc,log = self.run_command(["sbatch",slurm_script],return_output=True)
                slurm_ids.append(log[-1].split(' ')[-1])
                print (f'Job submission exit code : {rc}, slurm ID = {slurm_ids[-1]}')
            else:
                print(f'Submitting {len(slurm_config.inputParams)} jobs')
                print(f'Payload : {slurm_config.payload}')
                print(f'Param names : {slurm_config.inputParamsNames}')
                print ('Parameters :')
                for inputParam in slurm_config.inputParams:
                    print ("\t",inputParam)
                print('... don\'t worry, jobs not sent')

        print ('All slurm IDs :')
        print (' '.join(slurm_ids))

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
    parser.add_argument('--slurm',action='store',required=False,nargs='+',default=None,
                        help='Slurm parameters, separated by spaces and specified with = signs (eg stime=04:00:00 mem=6000)')
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

    # Slurm parameters #
    if args.slurm is not None:
        slurmParams = {arg.split("=")[0]:arg.split("=")[1] for arg in args.slurm}
        print ('Overwriting default slurm values below')
        for slurmParam,slurmVal in slurmParams.items():
            print (f'... {slurmParam:20s} = {slurmVal}')
    else:
        slurmParams = {}


    # Resubmit #
    if args.resubmit is not None:
        instance.findMissingJobs(args.resubmit)
        if args.submit is not None:
            instance.submit_on_slurm(args.submit,args.debug)
        sys.exit()
    # Submit #
    if args.submit is not None:
        instance.submit_on_slurm(args.submit,args.debug,slurmParams)
        sys.exit()

if __name__ == "__main__":
    main()
