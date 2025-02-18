# photon+X(->qq) code

## Setup
Download the code (to do once): 
`git clone https://github.com/lathomas/GammaXqq`

The next commands must be issued every time:
```
cd GammaXqq 
source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_103 x86_64-centos7-gcc12-opt
```

## Running the code
```
cd macros 
python3 analysis.py  -o output.root --year 2023 --era 'C' --isData 0 -i /pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/ZGto2QG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2810000/9650ee14-6c75-4e22-bff3-595197189178.root
```

## MC samples/data sets

The relevant MC samples for this analysis are: 

"Signal" (gamma+Z(qq)):
```
/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/ZGto2QG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/*/*.root
/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/ZGto2QG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/*/*.root
```
"Background" (gamma+jets):
```
/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/GJ_PTG-100to200_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v14-v3/*/*.root
/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/GJ_PTG-200to400_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v14-v3/*/*.root
/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/GJ_PTG-400to600_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v14-v3/*/*.root
/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/GJ_PTG-400to600_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15_ext1-v2/*/*.root
/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/GJ_PTG-600_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v14-v3/*/*.root
/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/GJ_PTG-600_TuneCP5_13p6TeV_amcatnlo-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15_ext1-v2/*/*.root
``` 

Data (eras 2024C, D, E, F, G, H, I): 
```
/pnfs/iihe/cms/ph/sc4/store/data/Run2024C/EGamma*/NANOAOD/PromptReco-v*/000/*/*/00000/*.root
/pnfs/iihe/cms/ph/sc4/store/data/Run2024D/EGamma*/NANOAOD/PromptReco-v*/000/*/*/00000/*.root
/pnfs/iihe/cms/ph/sc4/store/data/Run2024E/EGamma*/NANOAOD/PromptReco-v*/000/*/*/00000/*.root
/pnfs/iihe/cms/ph/sc4/store/data/Run2024F/EGamma*/NANOAOD/PromptReco-v*/000/*/*/00000/*.root
/pnfs/iihe/cms/ph/sc4/store/data/Run2024G/EGamma*/NANOAOD/PromptReco-v*/000/*/*/00000/*.root
/pnfs/iihe/cms/ph/sc4/store/data/Run2024H/EGamma*/NANOAOD/PromptReco-v*/000/*/*/00000/*.root
/pnfs/iihe/cms/ph/sc4/store/data/Run2024I/EGamma*/NANOAOD/PromptReco-v*/000/*/*/00000/*.root
``` 

## Submitting jobs to condor
Condor submission scripts are in the`condorsubmission` folder
```
cd condorsubmission
```
The following commands will create a dedicate folder (`outputdir`) from which the jobs will be submitted and where the output will be stored. 

**For MC:** 
```
sh SubmitToCondor.sh outputdir '/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/ZGto2QG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/*/*.root' zg 2023 C 0 
```
The arguments are, in that order:
- the output directory name (**The name must contain the sample exact process ID string** e.g. `GJ_PTG-400to600_TuneCP5_13p6TeV_amcatnlo-pythia8` as it is used by other scripts downstream for cross section normalisation)
- the list of files to process
- the process name
- the year to consider (2023 for MC simulation, 2024 for data)
- the era period ("C" for data, "BCD"/"E"/"F"/"G"/"H"/"I" for data) 
- an integer specifying whether this is data (1) or not (0) 


**For data**:
(mind that each era needs to be submitted separately since it needs different calibrations (example here for era G): 
```
sh SubmitToCondor.sh outputdir '/pnfs/iihe/cms/ph/sc4/store/data/Run2024G/EGamma*/NANOAOD/PromptReco-v*/000/*/*/00000/*.root' data 2024 G 1
```

### To monitor the jobs on condor:
```
condor_q
```


## Merging and scaling output to a given integrated luminosity
All the outputfiles from a set of folders (`folder1`, `folder2`, ...) can be merged and rescaled according to their cross section and the desired integrated luminosity. From the `macros` folder, do; 
```
python3  hadd_scale_merge_prov.py -d folder1 folder2  -l 100000 -o  myoutputfile.root
```
where the argument `-l` is the integrated luminosity in /pb. 
**The script looks for a match between the folder name and the dictionary defined in `xsmc.py` to assign the sample cross section. As discussed above, the folder name must contain the process ID string** 
