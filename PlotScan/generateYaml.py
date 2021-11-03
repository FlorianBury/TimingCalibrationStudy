import os
import sys
import yaml
import glob
import ROOT
import enlighten
import numpy as np

import time


def loopOverFiles(path,paramConv):
    rootfiles = list(glob.glob(os.path.join(path,'*.root')))
    pbar = enlighten.Counter(total=len(rootfiles), desc='Progress', unit='Files')
    files = []
    param_arr = np.zeros((len(rootfiles),len(paramConv)))
    start_time = time.time()
    times = []
    for idx,rootfile in enumerate(rootfiles):
        pbar.update()
        F = ROOT.TFile(rootfile,"READ")
        ROOT.SetOwnership(F,False) 
        files.append(os.path.basename(rootfile))
        param_arr[idx,:] = np.array([float(F.Get(key).GetTitle()) for key in paramConv.keys()])
        F.Close()

    assert ((param_arr.sum(1)==0)*1).sum() == 0. # check that there are no pure zero line
    # Order by ascending order #
    idx = np.lexsort([param_arr[:,i] for i in range(param_arr.shape[1])])
    param_arr = param_arr[idx]
    files_arr = np.array(files)[idx]
    content = {'path': os.path.abspath(path),
               'parameters': list(paramConv.values()),
               'files':{str(files_arr[i]):[float(p) for p in param_arr[i,:]] for i in range(param_arr.shape[0])}}
        
    with open('blank.yml','w') as handle:
        yaml.dump(content,handle,sort_keys=False)

paramConv = {'threshold':           'Threshold',
             'thresholdsmearing':   'Threshold Smearing',
             'tofsmearing':         'ToF Smearing'}

loopOverFiles(sys.argv[1],paramConv)
