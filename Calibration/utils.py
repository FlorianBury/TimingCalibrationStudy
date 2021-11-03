import os
import sys
import yaml
import math
import glob
import shutil
import datetime
import enlighten
from copy import deepcopy
from itertools import product

# Slurm configuration
from CP3SlurmUtils.Configuration import Configuration
from CP3SlurmUtils.SubmitWorker import SubmitWorker
from CP3SlurmUtils.Exceptions import CP3SlurmUtilsException

SlurmJobStatus = ["PENDING", "RUNNING", "COMPLETED", "FAILED", "COMPLETING", "CONFIGURING", "CANCELLED", "BOOT_FAIL", "NODE_FAIL", "PREEMPTED", "RESIZING", "SUSPENDED", "TIMEOUT", "OUT_OF_MEMORY", "REQUEUED", "unknown"]
SlurmJobStatus_active = ["CONFIGURING", "COMPLETING", "PENDING", "RUNNING", "RESIZING", "SUSPENDED", "REQUEUED"]                                                                                            
SlurmJobStatus_failed = ["FAILED", "TIMEOUT", "CANCELLED", "OUT_OF_MEMORY"]
SlurmJobStatus_completed = "COMPLETED"

def submit(yamlFile,name,N_per_job,output_files=[],keys_to_split=[],debug=False):
    # Load yaml #
    with open(yamlFile,'r') as f:
        runConfig = yaml.load(f,Loader=yaml.FullLoader)

    # Make slurm directory #
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    out_dir = os.path.abspath(os.path.dirname(__file__))
    slurm_working_dir = os.path.join(out_dir,'slurm',f'{name}_{timestamp}')
    infiles_dir = os.path.join(slurm_working_dir,'infiles')

    if not os.path.exists(slurm_working_dir):
        os.makedirs(slurm_working_dir)
    if os.path.exists(infiles_dir):
        shutil.rmtree(infiles_dir)
    os.makedirs(infiles_dir)

    # Split dictionnary #
    keys_to_repeat = [key for key in runConfig.keys() if key not in keys_to_split] 

    subDict = {k:v for k,v in runConfig.items() if k not in keys_to_split}
    varyDict = {}
    for key in keys_to_split:
        if isinstance(runConfig[key],list):
            subDict[key] = []
            varyDict[key] = runConfig[key]
        elif isinstance(runConfig[key],dict):
            subDict[key] = {}
            varyDict[key] = list(runConfig[key].items())
        else:
            raise RuntimeError(f'Item {key} is not iterable')

    N = 0
    counter = 0
    list_files = []
    combinations = list(product(*list(varyDict.values())))
    N_jobs = math.floor(len(combinations)/N_per_job)
    print (f'Saving configs in {infiles_dir}')
    pbar = enlighten.Counter(total=N_jobs, desc='Progress', unit='configs')
    for comb in combinations:
        counter += 1
        for key,val in zip(keys_to_split,comb):
            if isinstance(val,tuple):
                subDict[key][val[0]] = val[1]
            else:
                subDict[key].append(val)
        if counter == N_per_job:
            filename = os.path.join(infiles_dir,'config_{}.yml'.format(N))
            list_files.append(filename)
            with open(filename,'w') as handle:
                yaml.dump(subDict,handle)
            pbar.update()
            for key in keys_to_split:
                if isinstance(runConfig[key],list):
                    subDict[key] = []
                if isinstance(runConfig[key],dict):
                    subDict[key] = {}
            counter = 0
            N += 1

    # Slurm config #
    config = Configuration()
    config.sbatch_partition = 'cp3'
    config.sbatch_qos = 'cp3'
    config.sbatch_chdir = os.path.dirname(os.path.abspath(__file__))
    config.sbatch_time = '0-02:00:00'
    config.sbatch_memPerCPU= '2000'
    config.sbatch_additionalOptions = ["--export=ALL"]
    config.useJobArray = True
    config.inputParamsNames = []
    config.inputSandboxContent = [sys.argv[0]]
    config.inputParams = []

    config.payload = f"python3 {sys.argv[0]} --yaml ${{yaml}}"

    slurm_config = deepcopy(config)

    slurm_config.batchScriptsDir = os.path.join(slurm_working_dir, 'scripts')
    slurm_config.inputSandboxDir = slurm_config.batchScriptsDir
    slurm_config.stageoutDir = os.path.join(slurm_working_dir, 'output')
    slurm_config.stageoutLogsDir = os.path.join(slurm_working_dir, 'logs')
    slurm_config.stageoutFiles = output_files 

    slurm_config.inputParamsNames = ['yaml']
    slurm_config.inputParams = [[f] for f in list_files]

    # Submit job!
    print("Submitting job...")
    if not debug:
        submitWorker = SubmitWorker(slurm_config, submit=True, yes=True, debug=False, quiet=False)
        submitWorker()
        print("Done")
    else:
        print('Submitting {} jobs'.format(len(slurm_config.inputParams)))
        print('Param names : '+ ' --'.join([str(p) for p in slurm_config.inputParamsNames]))
        for param in slurm_config.inputParams:
            print ('...',param)
        print('... don\'t worry, jobs not sent')
    


