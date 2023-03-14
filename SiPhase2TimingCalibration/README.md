# CMSSW installation 

```
    # initialise cvmfs according to your system 
    cmsrel CMSSW_12_5_0
    cd CMSSW_12_5_0/src/
    git cms-addpkg SimTracker
    cd SimTracker
```

Add the current directory content to `SimTracker`

```
    cp -r * path_to_CMSSW/src/SimTracker/SiPhase2TimingCalibration/
```

then compile (either in `src` or in `SimTracker` or `SiPhase2TimingCalibration`, will just change the amount of things to compile)

```
    scram b -j 4
```

# Structure of the code 

```
 |-test
 | |-PythiaGunCalibration_cfg.py
 | |-Harvester_cfg.py
 | |-PUCalibration_cfg.py
 |-python
 | |-Phase2TrackerBXHistogram_cfi.py
 |-plugins
 | |-BuildFile.xml
 | |-Phase2TrackerBXHistogram.h
 | |-Phase2TrackerBXHistogram.cc
```



## plugins 

This is where the C++ classes lie, as well as the xml build file that links them. The instantiations are in the `.h` header, the methods and logics are in `.cc`

## python 

Contains the python config file that configures the plugin `Phase2TrackerBXHistogram`. Paremeters can be modified there, although for convenience the most important ones are overwritten in running configs (see `python`).

## test

This is where the `cmsRun` command can be run. The directory consists in several configs to run, more details follow.

### Files

- `PythiaGunCalibration_cfg.py` : config to run the production of DQM BX histograms, using a Pythia Gun of high-pT muons as input source for tracks and hits 
- `PUCalibration_cfg.py` : config to run the production of DQM BX histograms, using ttbar events with a pixing of minbias events
- `Harvester_cfg.py` : config to extract the histograms from the produced root files of the two other scripts (see below)

All these python configs can be run as (for options, see below)

```
    cmsRun <script>.py (options)
```

NB : for `PUCalibration_cfg.py` one needs to use root files for the hard scattering and the PU. One could use xrootd, but can be slow, best option is to use local copies of files. To avoid constant overwriting when collaborating, the absolute paths must be put in two txt files that are read by the config (and kept out of git synchronisation using gitignore) :
- `hard_scattering.txt` 
- `PU.txt`

Format within the file should be 
```
    file:<absolute_path_to_file>
    file:<absolute_path_to_file>
    [...]
```
No checks are done so make sure there are no spaces or empty lines

### Command and parameters

#### Running production 

The following parameters can be used in both production script (if not specified, a default value as written in the script is used) :
- `N` [int] : number of events to run over
- `pt` [float] : cut on the hit track pT in GeV
- `threshold` [float] : value fo the threshold in ADC for the hit detect logic (currently applied on both barrel and endcaps)
- `thresholdsmearing` [float] : value of smearing to be applied on the threshold (currently using uniform distribution) per detId (value is cached)
- `tofsmearing` [float] : value of smearing applied to the time-of-flight (currently using uniform distribution) per detId (value is cached)
- `mode` [string] : 'scan' [default] or 'emulate' [deprecated], specifies the mode (might be removed lated)
- `subdet` [string] : name of the part of the detector to use for the hits, allows to only target sub regions [default = 'ALL' = all tracker]
- `offset` [float] : specific offset value in 'emulate' mode, avoids the scan [default = -1 = do the scan]
- `verbose` [int] : different verbosity level for printing, from 0 (nothing) to 4 (super duper long printout)

How to use the parameters :

```
    cmsRun <script>.py parameter_name=parameter_value ...
```

example :

```
    cmsRun PUCalibration_cfg.py N=1 pt=2 threshold=6000
```

NB : the file will be created with the parameters used in the production in its name, with a suffix `raw` to indicate it has to go through the Harvester, otherwise the histograms are lost somewhere.

#### Runnign Harvester 

To run the Harvester, pass it the 'raw' file produced prior : 

```
    cmsRun Harvester_cfg.py inputs=<file>
```

NB : the file created will have the name `DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root`, better rename it to avoid overwriting.

# Content of the file

After the production and harvesting is done, here is the layout of the output file (`DQM...` or whatever it has been rename)

```
TFile**     DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root    
 TFile*     DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root    
  KEY: TDirectoryFile   DQMData;1   DQMData
```

The histograms are in `DQMData/Run 1/Ph2TkBXHist/Run summary`, can be accessed by ROOT TBrowser (faster on loca, after having downloaded the file from the server), or in the command prompt by

```
    root -l DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root
    _file0->cd("DQMData/Run 1/Ph2TkBXHist/Run summary")
``` 

Inside you will find two subdirectories : 
- `Hist1D` : the 1D BX histograms, for each delay (=offset) value and both Latched and Sampled mode. These are the simulated BX histograms that are the basis of this project, and also derived (in the future) from data
- `Hist2D` : several useful 2D histograms for easier visualisation (mode = Latched/Sampled)
    - `AttributionScan{mode}` : number of hist attribution to BX as a function of the delay (eg 0 = hit not identified in any BX, 2 = hit attributed to two different BX)
    - `HitsTrueNumberScan{mode}` : actual number of hits, useful for efficiency calculation
    - `OffsetScan{mode}` : aggregation of all 1D BX histograms as a 2D scan, for better visualisation
    - `HitsPositions3D` : (actually a 3D histogram) positions of all hits within the tracker
    - `HitsPositions2D` : same but projected on the y-z plane to show distribution in r and theta
    - `HitsPositions2DAbs` : same but using absolute values to restrict to a single quarter
