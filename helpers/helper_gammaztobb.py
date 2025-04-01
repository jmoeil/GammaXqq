import ROOT
import yaml
import json
import sys

from array import array
from math import floor, ceil
from scipy import stats
from scipy.signal import find_peaks
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
    
    #The subset of these jets which hold a specific flavour
    df = df.Define('Jet_TightID_Pt30_Central_partonFlavour', 'Jet_partonFlavour[Jet_TightID_Pt30_Central]')
    df = df.Define('Jet_TightID_Pt30_Central_PartonFlavour1', 'abs(Jet_TightID_Pt30_Central_partonFlavour)==1')
    df = df.Define('Jet_TightID_Pt30_Central_PartonFlavour2', 'abs(Jet_TightID_Pt30_Central_partonFlavour)==2')
    df = df.Define('Jet_TightID_Pt30_Central_PartonFlavour3', 'abs(Jet_TightID_Pt30_Central_partonFlavour)==3') 
    df = df.Define('Jet_TightID_Pt30_Central_PartonFlavour4', 'abs(Jet_TightID_Pt30_Central_partonFlavour)==4')
    df = df.Define('Jet_TightID_Pt30_Central_PartonFlavour5', 'abs(Jet_TightID_Pt30_Central_partonFlavour)==5')
    df = df.Define('Jet_TightID_Pt30_Central_PartonFlavour6', 'abs(Jet_TightID_Pt30_Central_partonFlavour)==6')
    #For now, consider only events with exactly 2 central jets and no other jet
    df = df.Filter('Sum(Jet_TightID_Pt30)==2&&Sum(Jet_TightID_Pt30_Central)==2','=2 central jets with pt>30 GeV, no additional jet')
    histos['photon_pt_2jselection'] = df.Histo1D(ROOT.RDF.TH1DModel('photon_pt_2jselection', '', 1000, 0, 1000), 'Photon_LooseID_Pt20_pt', 'LHEWeight_originalXWGTUP')

    #Compute the dijet invariant mass and remore invariant mass outside range of study
    df = df.Define('Mjj_Before_Deltas', 'InvariantMass(Jet_TightID_Pt30_Central_Pt[0], Jet_TightID_Pt30_Central_Eta[0], Jet_TightID_Pt30_Central_Phi[0], Jet_TightID_Pt30_Central_Mass[0], Jet_TightID_Pt30_Central_Pt[1], Jet_TightID_Pt30_Central_Eta[1], Jet_TightID_Pt30_Central_Phi[1], Jet_TightID_Pt30_Central_Mass[1])')
    df = df.Filter('Mjj_Before_Deltas < 200','Invariant mass clearly outside the range of this study')    

    # Delta eta
    df = df.Define('Jet_delta_eta','abs(Jet_TightID_Pt30_Central_Eta[0]-Jet_TightID_Pt30_Central_Eta[1])')
    histos['Jet_delta_eta'] = df.Histo1D(ROOT.RDF.TH1DModel('Jet_delta_eta','',100,0,5), 'Jet_delta_eta','LHEWeight_originalXWGTUP')
    # Delta phi
    df = df.Define('Jet_delta_phi','abs(acos(cos(Jet_TightID_Pt30_Central_Phi[0]-Jet_TightID_Pt30_Central_Phi[1])))')
    histos['Jet_delta_phi'] = df.Histo1D(ROOT.RDF.TH1DModel('Jet_delta_phi','',100,0,5), 'Jet_delta_phi','LHEWeight_originalXWGTUP')

    # Delta pT
    df = df.Define('Jet_pt1','Jet_TightID_Pt30_Central_Pt[0]')
    df = df.Define('Jet_pt2','Jet_TightID_Pt30_Central_Pt[1]')
    df = df.Define('Jet_delta_pT','abs(Jet_pt1-Jet_pt2)')
    histos['Jet_delta_pT'] = df.Histo1D(ROOT.RDF.TH1DModel('Jet_delta_pT','',1000,0,1000), 'Jet_delta_pT', 'LHEWeight_originalXWGTUP')

    # pT1/pT2
    df = df.Define('Jet_pT2pT1','Ratio_pt(Jet_pt1,Jet_pt2)')
    histos ['Jet_pT2pT1'] = df.Histo1D(ROOT.RDF.TH1DModel('Jet_pT2pT1','',100,0,1),'Jet_pT2pT1','LHEWeight_originalXWGTUP')

    # Angular distance R
    df = df.Define('Jet_delta_R','sqrt(pow(Jet_delta_eta,2) + pow(Jet_delta_phi,2))')
    histos['Jet_delta_R'] = df.Histo1D(ROOT.RDF.TH1DModel('Jet_delta_R','',100,0,5),'Jet_delta_R','LHEWeight_originalXWGTUP')

    # Veto Delta R(photon,j1) > 0.4
    df = df.Define('PJet_Delta_eta','abs(Photon_TightID_Pt100_eta[0]-Jet_TightID_Pt30_Central_Eta[0])') # Delta eta leading jet and photon
    df = df.Define('PJet_Delta_phi','abs(acos(cos(Photon_TightID_Pt100_phi[0]-Jet_TightID_Pt30_Central_Phi[0])))') # Delta phi leading jet and photon
    df = df.Define('PSubJet_Delta_eta','abs(Photon_TightID_Pt100_eta[0]-Jet_TightID_Pt30_Central_Eta[1])') # Delta eta subleading jet and photon
    df = df.Define('PSubJet_Delta_phi','abs(acos(cos(Photon_TightID_Pt100_phi[0]-Jet_TightID_Pt30_Central_Phi[1])))') # Delta phi subleading jet and photon
    
    df = df.Define('PJet_Delta_R','sqrt(pow(PJet_Delta_eta,2) + pow(PJet_Delta_phi,2))') # Delta R leading jet and photon
    df = df.Define('PSubJet_Delta_R','sqrt(pow(PSubJet_Delta_eta,2) + pow(PSubJet_Delta_phi,2))') # Delta R subleading jet and photon

    df = df.Filter('PJet_Delta_R > 0.4 && PSubJet_Delta_R > 0.4',f'Angular distance between the photon and both jets is > 0.4')
    
    # Plot these variables in histograms
    flavour_id = [k for k in range(1,6)]
    variables = {
    	'Jet_delta_eta': (100,0,5),
	'Jet_delta_phi': (100,0,5),
	'Jet_delta_pT': (1000,0,1000),
	'Jet_pT2pT1': (100,0,1),
	'Jet_delta_R': (100,0,5)
    }
    
    for fid in flavour_id:
    	# Create a filtered datafrom for each flavour
        flav_var = f'Jet_TightID_Pt30_Central_PartonFlavour{fid}'
        df_flav = df.Filter(f'Sum({flav_var})==2',f'Both jets have flavour {fid}')
	# Loop over each variable
        for varname, (nbins,xmin,xmax) in variables.items():
            hist_name = f'{varname}_PartonFlavour{fid}'
            histos[hist_name] = df_flav.Histo1D(ROOT.RDF.TH1DModel(hist_name,'',nbins,xmin,xmax),varname,'LHEWeight_originalXWGTUP')

    #Compute the dijet invariant mass 
    df = df.Define('Mjj', 'InvariantMass(Jet_TightID_Pt30_Central_Pt[0], Jet_TightID_Pt30_Central_Eta[0], Jet_TightID_Pt30_Central_Phi[0], Jet_TightID_Pt30_Central_Mass[0], Jet_TightID_Pt30_Central_Pt[1], Jet_TightID_Pt30_Central_Eta[1], Jet_TightID_Pt30_Central_Phi[1], Jet_TightID_Pt30_Central_Mass[1])')
    histos['mjj'] = df.Histo1D(ROOT.RDF.TH1DModel('mjj', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')    
    histos['mjj_partonflavour1'] = df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour1)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_partonflavour1', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')    
    histos['mjj_partonflavour2'] = df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour2)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_partonflavour2', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')    
    histos['mjj_partonflavour3'] = df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour3)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_partonflavour3', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')    
    histos['mjj_partonflavour4'] = df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour4)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_partonflavour4', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')    
    histos['mjj_partonflavour5'] = df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour5)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_partonflavour5', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')    
    histos['mjj_partonflavour6'] = df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour6)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_partonflavour6', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP') # Create an histogram for two top jets. Should be empty. This line is just a sanity check.
    histos['mjj_partonflavour1_bis'] = df.Filter('Jet_TightID_Pt30_Central_PartonFlavour1[0]&&Jet_TightID_Pt30_Central_PartonFlavour1[1]').Histo1D(ROOT.RDF.TH1DModel('mjj_partonflavour1_bis', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')

   # --- Mjj after vetoes: Delta R < 2 and pT2/pT1 > 0.4 ---
    vetoed_df = df.Filter('Jet_delta_R > 0.75 && Jet_delta_R < 2 && Jet_pT2pT1 > 0.4', 'DeltaR < 2 and pT2/pT1 > 0.4')

   # Flavour-filtered Mjj histograms AFTER vetoes
    histos['mjj_vetoed'] = vetoed_df.Histo1D(ROOT.RDF.TH1DModel('mjj_vetoed', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')
    histos['mjj_vetoed_flavour1'] = vetoed_df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour1)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_vetoed_flavour1', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')
    histos['mjj_vetoed_flavour2'] = vetoed_df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour2)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_vetoed_flavour2', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')
    histos['mjj_vetoed_flavour3'] = vetoed_df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour3)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_vetoed_flavour3', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')
    histos['mjj_vetoed_flavour4'] = vetoed_df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour4)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_vetoed_flavour4', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')
    histos['mjj_vetoed_flavour5'] = vetoed_df.Filter('Sum(Jet_TightID_Pt30_Central_PartonFlavour5)==2').Histo1D(ROOT.RDF.TH1DModel('mjj_vetoed_flavour5', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP') 
   
    # Define pT for each flavour

    df = df.Define('Jet_TightID_Pt30_Central_Pt_Flavour1','Jet_TightID_Pt30_Central_Pt[Jet_TightID_Pt30_Central_PartonFlavour1]')
    df = df.Define('Jet_TightID_Pt30_Central_Pt_Flavour2','Jet_TightID_Pt30_Central_Pt[Jet_TightID_Pt30_Central_PartonFlavour2]')
    df = df.Define('Jet_TightID_Pt30_Central_Pt_Flavour3','Jet_TightID_Pt30_Central_Pt[Jet_TightID_Pt30_Central_PartonFlavour3]')
    df = df.Define('Jet_TightID_Pt30_Central_Pt_Flavour4','Jet_TightID_Pt30_Central_Pt[Jet_TightID_Pt30_Central_PartonFlavour4]')
    df = df.Define('Jet_TightID_Pt30_Central_Pt_Flavour5','Jet_TightID_Pt30_Central_Pt[Jet_TightID_Pt30_Central_PartonFlavour5]') 

    # Plot pT for each flavour

    histos['jet_pt_flavour1'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_pt_partonflavour1', '', 100, 0, 500),'Jet_TightID_Pt30_Central_Pt_Flavour1','LHEWeight_originalXWGTUP')
    histos['jet_pt_flavour2'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_pt_partonflavour2', '', 100, 0, 500),'Jet_TightID_Pt30_Central_Pt_Flavour2','LHEWeight_originalXWGTUP')
    histos['jet_pt_flavour3'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_pt_partonflavour3', '', 100, 0, 500),'Jet_TightID_Pt30_Central_Pt_Flavour3','LHEWeight_originalXWGTUP')
    histos['jet_pt_flavour4'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_pt_partonflavour4', '', 100, 0, 500),'Jet_TightID_Pt30_Central_Pt_Flavour4','LHEWeight_originalXWGTUP')
    histos['jet_pt_flavour5'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_pt_partonflavour5', '', 100, 0, 500),'Jet_TightID_Pt30_Central_Pt_Flavour5','LHEWeight_originalXWGTUP')
 
    # Define eta for each flavour
    df = df.Define('Jet_TightID_Pt30_Central_Eta_Flavour1','Jet_TightID_Pt30_Central_Eta[Jet_TightID_Pt30_Central_PartonFlavour1]')
    df = df.Define('Jet_TightID_Pt30_Central_Eta_Flavour2','Jet_TightID_Pt30_Central_Eta[Jet_TightID_Pt30_Central_PartonFlavour2]')
    df = df.Define('Jet_TightID_Pt30_Central_Eta_Flavour3','Jet_TightID_Pt30_Central_Eta[Jet_TightID_Pt30_Central_PartonFlavour3]')
    df = df.Define('Jet_TightID_Pt30_Central_Eta_Flavour4','Jet_TightID_Pt30_Central_Eta[Jet_TightID_Pt30_Central_PartonFlavour4]')
    df = df.Define('Jet_TightID_Pt30_Central_Eta_Flavour5','Jet_TightID_Pt30_Central_Eta[Jet_TightID_Pt30_Central_PartonFlavour5]')

    # Plot eta for each flavour
    histos['jet_eta_flavour1'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_eta_partonflavour1', '', 50, -2.5, 2.5),'Jet_TightID_Pt30_Central_Eta_Flavour1','LHEWeight_originalXWGTUP')
    histos['jet_eta_flavour2'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_eta_partonflavour2', '', 50, -2.5, 2.5),'Jet_TightID_Pt30_Central_Eta_Flavour2','LHEWeight_originalXWGTUP')
    histos['jet_eta_flavour3'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_eta_partonflavour3', '', 50, -2.5, 2.5),'Jet_TightID_Pt30_Central_Eta_Flavour3','LHEWeight_originalXWGTUP')
    histos['jet_eta_flavour4'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_eta_partonflavour4', '', 50, -2.5, 2.5),'Jet_TightID_Pt30_Central_Eta_Flavour4','LHEWeight_originalXWGTUP')
    histos['jet_eta_flavour5'] = df.Histo1D(ROOT.RDF.TH1DModel('jet_eta_partonflavour5', '', 50, -2.5, 2.5),'Jet_TightID_Pt30_Central_Eta_Flavour5','LHEWeight_originalXWGTUP')

    # Count the number of events per flavour
    n1 = histos['mjj_partonflavour1'].Integral(60, 120)
    n2 = histos['mjj_partonflavour2'].Integral(60, 120)
    n3 = histos['mjj_partonflavour3'].Integral(60, 120)
    n4 = histos['mjj_partonflavour4'].Integral(60, 120)
    n5 = histos['mjj_partonflavour5'].Integral(60, 120)
    n = np.array([n1,n2,n3,n4,n5])
    total = sum(n)

    frac1 = n1/total # Valeur th√orique : 0.115
    frac2 = n2/total # Valeur th√orique : 0.156
    frac3 = n3/total # Valeur th√orique : 0.156
    frac4 = n4/total # Valeur th√orique : 0.115
    frac5 = n5/total # Valeur th√orique : 0.151
    BR_Obs = np.array([frac1,frac2,frac3,frac4,frac5])
    BR_Theo = np.array([0.156, 0.116, 0.156, 0.116, 0.156])  # d, u, s, c, b, t
    BR_Theo = BR_Theo * (1/np.sum(BR_Theo))
    sigma = np.sqrt(n)/total 
    ndof = len(BR_Obs) - 1
    difference = BR_Obs-BR_Theo
    chi2 = np.sum((BR_Obs-BR_Theo)**2/BR_Theo)
    p = 1 - stats.chi2.cdf(chi2,ndof)   
    
    # Get the number of bins

    #n_bins1 = histos['mjj_partonflavour1'].GetNbinsX()
    #n_bins2 = histos['mjj_partonflavour2'].GetNbinsX()
    #n_bins3 = histos['mjj_partonflavour3'].GetNbinsX()
    #n_bins4 = histos['mjj_partonflavour4'].GetNbinsX()
    #n_bins5 = histos['mjj_partonflavour5'].GetNbinsX()
    
    # Extract the bin centers and their contents (values)

    #bin_contents1 = np.array([histos['mjj_partonflavour1'].GetBinContent(i+1) for i in range(n_bins1)])  # Histogram values
    #bin_contents2 = np.array([histos['mjj_partonflavour2'].GetBinContent(i+1) for i in range(n_bins2)])  # Histogram values
    #bin_contents3 = np.array([histos['mjj_partonflavour3'].GetBinContent(i+1) for i in range(n_bins3)])  # Histogram values
    #bin_contents4 = np.array([histos['mjj_partonflavour4'].GetBinContent(i+1) for i in range(n_bins4)])  # Histogram values
    #bin_contents5 = np.array([histos['mjj_partonflavour5'].GetBinContent(i+1) for i in range(n_bins5)])  # Histogram values
    # Use find_peaks to extract the peaks

    #peaks1, properties1 = find_peaks(bin_contents1, height=300, prominence=100)  # Ajuste ces seuils selon tes donn√©e
    #peaks2, properties2 = find_peaks(bin_contents2, height=300, prominence=100)  # Ajuste ces seuils selon tes donn√©e
    #peaks3, properties3 = find_peaks(bin_contents3, height=300, prominence=100)  # Ajuste ces seuils selon tes donn√©e
    #peaks4, properties4 = find_peaks(bin_contents4, height=300, prominence=100)  # Ajuste ces seuils selon tes donn√©e
    #peaks5, properties5 = find_peaks(bin_contents5, height=300, prominence=100)  # Ajuste ces seuils selon tes donn√©
    #peaks = np.array([

 
    return df, histos, BR_Obs, n, total, ndof, chi2, difference, p
