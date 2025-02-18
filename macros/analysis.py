from datetime import datetime
import ROOT
import os
import sys
import argparse

#Importing stuff from other python files
sys.path.insert(0, '../helpers')

import helper_gammaztobb as h_gammaztobb

ROOT.gInterpreter.Declare('#include "../helpers/Helper.h"')

def main():
    ###Arguments 
    parser = argparse.ArgumentParser(
        description='''Photon+X analysis''',
        usage='use "%(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--max_events", dest="max_events", help="Maximum number of events to analyze. Default=-1 i.e. run on all events.", type=int, default=-1)
    parser.add_argument("-i", "--input", dest="inputFile", help="Input file", type=str, default='')
    parser.add_argument("-o", "--output", dest="outputFile", help="Output file", type=str, default='')
    parser.add_argument("--year", dest="year", help="Year considered (2022, 2023, 2024)", type=int, default=2023)
    parser.add_argument("--era", dest="era", help="Era", type=str, default='Cv4')
    parser.add_argument("--isData", dest="isData", help="is Data or MC", type=int, default=0)
    args = parser.parse_args()
    
    ###Define the RDataFrame from the input tree
    inputFile = args.inputFile
    if inputFile == '':
        inputFile = '/pnfs/iihe/cms/ph/sc4/store/mc/Run3Summer23NanoAODv12/ZGto2QG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/NANOAODSIM/130X_mcRun3_2023_realistic_v15-v2/2810000/9650ee14-6c75-4e22-bff3-595197189178.root'
        
    #Load the TTree and make a RDataFrame, see https://root.cern/doc/v628/classROOT_1_1RDataFrame.html
    df = ROOT.RDataFrame('Events', inputFile)
    
    #Event weight. If not defined (e.g. for data), set it to 1.
    if not 'LHEWeight_originalXWGTUP' in df.GetColumnNames():
        df = df.Define('LHEWeight_originalXWGTUP','return 1.0;')
    #Declare a few branches when they are missing
    if not 'Jet_btagPNetB' in df.GetColumnNames():
        df = df.Define('Jet_btagPNetB','Jet_btagDeepFlavB')
    if not 'Electron_mvaIso_WPHZZ' in df.GetColumnNames():
        df = df.Define('Electron_mvaIso_WPHZZ','Electron_mvaIso_WP90')
    if not 'HLT_Photon50EB_TightID_TightIso' in df.GetColumnNames():
        df = df.Define('HLT_Photon50EB_TightID_TightIso', 'return 1.;')
    if not 'HLT_Photon30EB_TightID_TightIso' in df.GetColumnNames():
        df = df.Define('HLT_Photon30EB_TightID_TightIso', 'return 1.;')
    if not 'HLT_Photon45EB_TightID_TightIso' in df.GetColumnNames():
        df = df.Define('HLT_Photon45EB_TightID_TightIso','HLT_Photon30EB_TightID_TightIso')
       
    
    #Example to make a histogram with the distribution of the number of vertices
    nvtx_histo = df.Histo1D(ROOT.RDF.TH1DModel("h_nvtx" , "Number of reco vertices;N_{vtx};Events"  ,    100, 0., 100.), "PV_npvs","LHEWeight_originalXWGTUP")
    nEvents = df.Count().GetValue()
    print('There are {} events'.format(nEvents))
    print('File is: ',inputFile )
    
    #Max events to run on 
    max_events = min(nEvents, args.max_events) if args.max_events >=0 else nEvents
    df = df.Range(0, max_events)
    #Next line to monitor event loop progress
    df = df.Filter('if(tdfentry_ %100000 == 0) {cout << "Event is  " << tdfentry_ << endl;} return true;')

    #Next few lines apply some cleaning to reject problematic events/data. Do not remove
    df = df.Filter('Flag_HBHENoiseFilter&&Flag_HBHENoiseIsoFilter&&Flag_goodVertices&&Flag_EcalDeadCellTriggerPrimitiveFilter&&Flag_BadPFMuonFilter&&Flag_BadPFMuonDzFilter')
    df = df.Filter('run<379344||run>379411') #Pixels off, at least for some of these runs

    #Output file 
    if args.outputFile == '':
        args.outputFile = 'output_'+args.channel+'.root'
    out = ROOT.TFile(args.outputFile, "recreate")
    ####The sequence of filters/column definition starts here

    #Everything is done in h_gammaztobb
    df, histos = h_gammaztobb.GammaZSelection(df, args.year, args.era, args.isData)
    df_report = df.Report()
    for i in histos:
        histos[i].GetValue().Write()

    df_report.Print()

    nvtx_histo.GetValue().Write()

if __name__ == '__main__':
    main()
