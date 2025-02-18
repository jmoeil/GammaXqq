##XS in pb
xs_gjets = {
'GJ-4Jets_dRGJ-0p25_PTG-100to200_HT-40to200_TuneCP5_13p6TeV_madgraphMLM-pythia8':5.554e+02,
'GJ-4Jets_dRGJ-0p25_PTG-100to200_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8':2.043e+02,
'GJ-4Jets_dRGJ-0p25_PTG-100to200_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8':2.982e+01,
'GJ-4Jets_dRGJ-0p25_PTG-100to200_HT-600to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8':9.683e+00,
'GJ-4Jets_dRGJ-0p25_PTG-100to200_HT-1000_TuneCP5_13p6TeV_madgraphMLM-pythia8':1.629e+00,
'GJ-4Jets_dRGJ-0p25_PTG-10to100_HT-40to100_TuneCP5_13p6TeV_madgraphMLM-pythia8':1.226e+05,
'GJ-4Jets_dRGJ-0p25_PTG-10to100_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8':3.224e+04,
'GJ-4Jets_dRGJ-0p25_PTG-10to100_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8':5.535e+03,
'GJ-4Jets_dRGJ-0p25_PTG-200_HT-40to400_TuneCP5_13p6TeV_madgraphMLM-pythia8':4.367e+01,
'GJ-4Jets_dRGJ-0p25_PTG-200_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8':1.173e+01,
'GJ-4Jets_dRGJ-0p25_PTG-200_HT-600to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8':4.771e+00,
'GJ-4Jets_dRGJ-0p25_PTG-200_HT-1000_TuneCP5_13p6TeV_madgraphMLM-pythia8':1.025e+00,
'GJ_PTG-20to100_ETAG-2p0_TuneCP5_13p6TeV_amcatnlo-pythia8':1.951e+05,
'GJ_PTG-100to200_TuneCP5_13p6TeV_amcatnlo-pythia8':1.380e+03,
'GJ_PTG-200to400_TuneCP5_13p6TeV_amcatnlo-pythia8':8.792e+01,
'GJ_PTG-400to600_TuneCP5_13p6TeV_amcatnlo-pythia8':3.809e+00,
'GJ_PTG-600_TuneCP5_13p6TeV_amcatnlo-pythia8':5.741e-01,
'GJ-4Jets_dRGJ-0p25_PTG-30to100_ETAG-2p0_TuneCP5_13p6TeV_madgraphMLM-pythia8':2.576e+04, #it seems there's no upper cut on ptg despite the name says otherwise...
'GJ-4Jets_dRGJ-0p25_PTG-100_ETAG-2p0_TuneCP5_13p6TeV_madgraphMLM-pythia8':5.820e+02,
#'GJ-4Jets_dRGJ-0p25_PTG-30to100_ETAG-2p0_TuneCP5_13p6TeV_madgraphMLM-pythia8':1.030e+05, #before filter effcy 
#'GJ-4Jets_dRGJ-0p25_PTG-100_ETAG-2p0_TuneCP5_13p6TeV_madgraphMLM-pythia8':4.672e+03, #before filter effcy 
'WGto2QG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8':6.349e-01,
'WGto2QG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8':4.021e+00,
'ZGto2QG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8':1.952e+00,
'ZGto2QG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8':2.842e-01,
'WGto2QG-1Jets_PTG-10_TuneCP5_13p6TeV_amcatnloFXFX-pythia8':2.943e+02,
'ZGto2QG-1Jets_PTG-10_TuneCP5_13p6TeV_amcatnloFXFX-pythia8':1.424e+02, 
'TTtoLNu2Q':923.6*0.6832*(1-0.6832)*2,
'TTto2L2Nu':923.6*(1-0.6832)*(1-0.6832),
'TbarWplustoLNu2Q':87.9/2*0.6832*(1-0.6832)*2,
'TWminustoLNu2Q':87.9/2*0.6832*(1-0.6832)*2,
'WtoLNu-2Jets':9009.5*(1-0.6832),
'ZG2JtoG2L2J_EWK_MLL-50_MJJ-120_TuneCP5_13p6TeV_madgraph-pythia8':1.136e-01, #from XSGenAnalyzer
'VBFtoG_PTG-10to100_TuneCP5_13p6TeV_madgraph-pythia8':6.235e+02, #from XSGenAnalyzer 
'VBFtoG_PTG-100to200_TuneCP5_13p6TeV_madgraph-pythia8':7.573,#from XSGenAnalyzer
'VBFtoG_PTG-200_TuneCP5_13p6TeV_madgraph-pythia8':1.077,#from XSGenAnalyzer
'VBFZG':0.195,#from Si Hyun, 
'VBFHG':0.19/0.832*0.57 #from Si Hyun
}

#https://twiki.cern.ch/twiki/bin/viewauth/CMS/MATRIXCrossSectionsat13p6TeV
#https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopNNLORef
#https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
