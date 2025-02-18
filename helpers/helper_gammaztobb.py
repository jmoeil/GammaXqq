import ROOT
import yaml
import json
import sys

from array import array
from math import floor, ceil
import numpy as np

#For JECs
from corrections_modified import *

def GammaZSelection(df, year=2023, era='C', isData=False):
    '''
    Select events with = 1 photon with pT>100 GeV.
    The event must pass a single photon trigger. 
    Requires exactly 2 jets
    '''
    histos = {}
    #Trigger
    df = df.Filter('HLT_Photon50EB_TightID_TightIso||HLT_Photon45EB_TightID_TightIso','HLT_Photon50EB_TightID_TightIso||HLT_Photon45EB_TightID_TightIso')


    #Offline photon, define a loose ID for pt>20. This will be used for vetoing additional photons 
    df = df.Define('Photon_LooseID_Pt20','Photon_mvaID_WP90&&abs(Photon_eta)<1.4442&&Photon_pt>20') 
    df = df.Define('Photon_LooseID_Pt20_pt','Photon_pt[Photon_LooseID_Pt20]') 
    df = df.Define('Photon_LooseID_Pt20_eta','Photon_eta[Photon_LooseID_Pt20]') 
    df = df.Define('Photon_LooseID_Pt20_phi','Photon_phi[Photon_LooseID_Pt20]') 
    
    #Offline photon, define a tight ID for the photon of interest, require pt > 100 GeV for now
    df = df.Define('Photon_TightID_Pt100','Photon_mvaID_WP80&&abs(Photon_eta)<1.4442&&Photon_pt>100&&Photon_electronVeto&&!Photon_pixelSeed') 
    df = df.Define('Photon_TightID_Pt100_pt','Photon_pt[Photon_TightID_Pt100]') 
    df = df.Define('Photon_TightID_Pt100_eta','Photon_eta[Photon_TightID_Pt100]') 
    df = df.Define('Photon_TightID_Pt100_phi','Photon_phi[Photon_TightID_Pt100]') 

    #Plot photon pt. Only filters called before this line applies for this plot
    histos['photon_pt_aftertrigger'] = df.Histo1D(ROOT.RDF.TH1DModel('photon_pt_aftertrigger', '', 1000, 0, 1000), 'Photon_LooseID_Pt20_pt', 'LHEWeight_originalXWGTUP')

    #Some simple event filters:
    df = df.Filter('Sum(Photon_LooseID_Pt20)==1','=1 loose photon with p_{T}>20 GeV')
    df = df.Filter('Sum(Photon_TightID_Pt100)==1','=1 tight photon with p_{T}>100 GeV')
    #Plot photon pt. Because of the two filters above, this distribution starts at 100 GeV.
    histos['photon_pt_afterptcut'] = df.Histo1D(ROOT.RDF.TH1DModel('photon_pt_afterptcut', '', 1000, 0, 1000), 'Photon_LooseID_Pt20_pt', 'LHEWeight_originalXWGTUP')
    
    #More filtering: electron/muon veto (to reject processes like ttbar) 

    #Electron veto
    df = df.Define('Electron_LooseID_Pt15','Electron_pt>15&&Electron_mvaIso_WPHZZ')    
    df = df.Filter('Sum(Electron_LooseID_Pt15)==0','0 good electron with p_{T}>15 GeV')
    
    #Muon veto
    df = df.Define('Muon_LooseID_Pt10','Muon_pfIsoId>=2&&Muon_mediumPromptId&&Muon_pt>10')
    df = df.Filter('Sum(Muon_LooseID_Pt10)==0','Sum(Muon_LooseID_Pt10)==0')

    #Jet selection
    # Apply the newest jet energy corrections: 
    JECfile, corrfile = JECsInit(year, era, isData)
    setupjecs(JECfile, corrfile)
    df = df.Redefine('Jet_pt', 'JetCorPt(Jet_area, Jet_eta, Jet_phi, Jet_pt, Jet_rawFactor, Rho_fixedGridRhoFastjetAll,'+str(isData)+')')
    #The following consideres only jets that are not pathological and do not have a large muon or "charged electromagnetic" energy fraction (i.e. the jet is not made mostly of a muon or an electron) 
    #Also remove jets mostly of neutral EM energy (= the photon)
    df = df.Define('Jet_TightID_Pt30', 'Jet_jetId>=4&&Jet_muEF<0.5&&Jet_chEmEF<0.5&&Jet_neEmEF<0.9&&Jet_pt>30')
    #The subset of these jets that are central 
    df = df.Define('Jet_TightID_Pt30_Central', 'Jet_jetId>=4&&Jet_muEF<0.5&&Jet_chEmEF<0.5&&Jet_neEmEF<0.9&&Jet_pt>30&&abs(Jet_eta)<2.4')
    df = df.Define('Jet_TightID_Pt30_Central_Pt', 'Jet_pt[Jet_TightID_Pt30_Central]')
    df = df.Define('Jet_TightID_Pt30_Central_Eta', 'Jet_eta[Jet_TightID_Pt30_Central]')
    df = df.Define('Jet_TightID_Pt30_Central_Phi', 'Jet_phi[Jet_TightID_Pt30_Central]')
    df = df.Define('Jet_TightID_Pt30_Central_Mass', 'Jet_mass[Jet_TightID_Pt30_Central]')
    
    #For now, consider only events with exactly 2 central jets and no other jet
    df = df.Filter('Sum(Jet_TightID_Pt30)==2&&Sum(Jet_TightID_Pt30_Central)==2','=2 central jets with pt>30 GeV, no additional jet')
    histos['photon_pt_2jselection'] = df.Histo1D(ROOT.RDF.TH1DModel('photon_pt_2jselection', '', 1000, 0, 1000), 'Photon_LooseID_Pt20_pt', 'LHEWeight_originalXWGTUP')

    #Compute the dijet invariant amss 
    df = df.Define('Mjj', 'InvariantMass(Jet_TightID_Pt30_Central_Pt[0], Jet_TightID_Pt30_Central_Eta[0], Jet_TightID_Pt30_Central_Phi[0], Jet_TightID_Pt30_Central_Mass[0], Jet_TightID_Pt30_Central_Pt[1], Jet_TightID_Pt30_Central_Eta[1], Jet_TightID_Pt30_Central_Phi[1], Jet_TightID_Pt30_Central_Mass[1])')
    histos['mjj'] = df.Histo1D(ROOT.RDF.TH1DModel('mjj', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')    

    return df, histos


