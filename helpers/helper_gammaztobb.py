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

variables = {
        'Jet_delta_eta': (100,0,5),
        'Jet_delta_phi': (100,0,5),
        'Jet_delta_pT': (1000,0,1000),
        'Jet_pT2pT1': (100,0,1),
        'Jet_delta_R': (100,0,5)
    }
cut_conditions = {
	'pt2pt1': 'Jet_pT2pT1 > 0.02',
	'Delta_R': 'Jet_delta_R < 3.95'
}

def plot_jet_kinematics_by_flavour(df, histos, label_suffix=''):
    """
    Create pt and eta histograms of jets, separated by parton flavour (1 to 5).
    Optionally, a label_suffix can be passed to distinguish histos after cuts.
    """
    for k in range(1, 6):
        pt_col = f'Jet_TightID_Pt30_Central_Pt_Flavour{k}{label_suffix}'
        eta_col = f'Jet_TightID_Pt30_Central_Eta_Flavour{k}{label_suffix}'

        df = df.Define(pt_col, f'Jet_TightID_Pt30_Central_Pt[Jet_TightID_Pt30_Central_PartonFlavour{k}]')
        df = df.Define(eta_col, f'Jet_TightID_Pt30_Central_Eta[Jet_TightID_Pt30_Central_PartonFlavour{k}]')

        pt_hist_name = f'jet_pt_partonflavour{k}{label_suffix}'
        eta_hist_name = f'jet_eta_partonflavour{k}{label_suffix}'

        histos[pt_hist_name] = df.Histo1D(ROOT.RDF.TH1DModel(pt_hist_name, '', 100, 0, 500), pt_col, 'LHEWeight_originalXWGTUP')
        histos[eta_hist_name] = df.Histo1D(ROOT.RDF.TH1DModel(eta_hist_name, '', 50, -2.5, 2.5), eta_col, 'LHEWeight_originalXWGTUP')
    
    return df

def cut_fill_histos(df, condition_expr,label):
	histos_cut = {}
	df_cut = df.Filter(condition_expr,label)

	# Inclusive mjj
	histos_cut[f'mjj_cut_{label}'] = df_cut.Histo1D(ROOT.RDF.TH1DModel(f'mjj_cut_{label}', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')

	# Flavour-filtered mjj
	for k in range(1,6):
	   flav_label = f'PartonFlavour_{k}_cut_{label}'
	   histos_cut[f'mjj_{flav_label}'] = df_cut.Filter(f'Sum(Jet_TightID_Pt30_Central_PartonFlavour{k})==2').Histo1D(ROOT.RDF.TH1DModel(f'mjj_{flav_label}', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')

	# Inclusive Kinematic histograms
	for varname,(nbins,xmin,xmax) in variables.items():
	   histos_cut[f'{varname}_cut_{label}'] = df_cut.Histo1D(ROOT.RDF.TH1DModel(f'{varname}_cut_{label}','',nbins,xmin,xmax),varname,'LHEWeight_originalXWGTUP')

	# Flavour-filtered kinematic histograms
	for varname,(nbins,xmin,xmax) in variables.items():
	   for k in range(1,6):
	       histos_cut[f'{varname}_PartonFlavour{k}_cut_{label}'] = df_cut.Histo1D(ROOT.RDF.TH1DModel(f'{varname}_PartonFlavour{k}_cut_{label}','',nbins,xmin,xmax),varname,'LHEWeight_originalXWGTUP')

	return histos_cut

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
    for k in range(1,7):
        df = df.Define(f'Jet_TightID_Pt30_Central_PartonFlavour{k}',f'abs(Jet_TightID_Pt30_Central_partonFlavour)=={k}')
    #For now, consider only events with exactly 2 central jets and no other jet
    df = df.Filter('Sum(Jet_TightID_Pt30)==2&&Sum(Jet_TightID_Pt30_Central)==2','=2 central jets with pt>30 GeV, no additional jet')
    histos['photon_pt_2jselection'] = df.Histo1D(ROOT.RDF.TH1DModel('photon_pt_2jselection', '', 1000, 0, 1000), 'Photon_LooseID_Pt20_pt', 'LHEWeight_originalXWGTUP')

    #Compute the dijet invariant mass and remove invariant mass mjj outside range of study
    df = df.Define('Mjj_Before_Deltas', 'InvariantMass(Jet_TightID_Pt30_Central_Pt[0], Jet_TightID_Pt30_Central_Eta[0], Jet_TightID_Pt30_Central_Phi[0], Jet_TightID_Pt30_Central_Mass[0], Jet_TightID_Pt30_Central_Pt[1], Jet_TightID_Pt30_Central_Eta[1], Jet_TightID_Pt30_Central_Phi[1], Jet_TightID_Pt30_Central_Mass[1])')
    df = df.Filter('Mjj_Before_Deltas < 200','Invariant mass clearly outside the range of this study')    
    
    # --- Define key kinematic variables to study their behaviors ---
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

    # Veto Delta R(photon,j) > 0.4 with j = leading and subleading jet
    df = df.Define('PJet_Delta_eta','abs(Photon_TightID_Pt100_eta[0]-Jet_TightID_Pt30_Central_Eta[0])') # Delta eta leading jet and photon
    df = df.Define('PJet_Delta_phi','abs(acos(cos(Photon_TightID_Pt100_phi[0]-Jet_TightID_Pt30_Central_Phi[0])))') # Delta phi leading jet and photon
    df = df.Define('PSubJet_Delta_eta','abs(Photon_TightID_Pt100_eta[0]-Jet_TightID_Pt30_Central_Eta[1])') # Delta eta subleading jet and photon
    df = df.Define('PSubJet_Delta_phi','abs(acos(cos(Photon_TightID_Pt100_phi[0]-Jet_TightID_Pt30_Central_Phi[1])))') # Delta phi subleading jet and photon
    
    df = df.Define('PJet_Delta_R','sqrt(pow(PJet_Delta_eta,2) + pow(PJet_Delta_phi,2))') # Delta R leading jet and photon
    df = df.Define('PSubJet_Delta_R','sqrt(pow(PSubJet_Delta_eta,2) + pow(PSubJet_Delta_phi,2))') # Delta R subleading jet and photon

    df = df.Filter('PJet_Delta_R > 0.4 && PSubJet_Delta_R > 0.4',f'Angular distance between the photon and both jets is > 0.4')
    
    # Plot these variables in histograms
    
    for k in range(1,6):
    	# Create a filtered datafrom for each flavour
        flav_var = f'Jet_TightID_Pt30_Central_PartonFlavour{k}'
        df_flav = df.Filter(f'Sum({flav_var})==2',f'Both jets have flavour {k}')
	# Loop over each variable
        for varname, (nbins,xmin,xmax) in variables.items():
            hist_name = f'{varname}_PartonFlavour{k}'
            histos[hist_name] = df_flav.Histo1D(ROOT.RDF.TH1DModel(hist_name,'',nbins,xmin,xmax),varname,'LHEWeight_originalXWGTUP')

    #Compute the dijet invariant mass 
    df = df.Define('Mjj', 'InvariantMass(Jet_TightID_Pt30_Central_Pt[0], Jet_TightID_Pt30_Central_Eta[0], Jet_TightID_Pt30_Central_Phi[0], Jet_TightID_Pt30_Central_Mass[0], Jet_TightID_Pt30_Central_Pt[1], Jet_TightID_Pt30_Central_Eta[1], Jet_TightID_Pt30_Central_Phi[1], Jet_TightID_Pt30_Central_Mass[1])')
    
    histos['mjj'] = df.Histo1D(ROOT.RDF.TH1DModel('mjj', '', 1000, 0, 1000), 'Mjj', 'LHEWeight_originalXWGTUP')
    for k in range(1,7):
        histos[f'mjj_partonflavour{k}'] = df.Filter(f'Sum(Jet_TightID_Pt30_Central_PartonFlavour{k})==2').Histo1D(ROOT.RDF.TH1DModel(f'mjj_partonflavour{k}','',1000,0,1000),'Mjj','LHEWeight_originalXWGTUP')

    # Compute and add histograms for the different cuts 
    for label,condition in cut_conditions.items():
        extra_histos = cut_fill_histos(df, condition,label)
        histos.update(extra_histos)

    df = plot_jet_kinematics_by_flavour(df, histos)
    
    # Count the number of events per flavour
    n1 = histos['mjj_partonflavour1'].Integral(60, 120)
    n2 = histos['mjj_partonflavour2'].Integral(60, 120)
    n3 = histos['mjj_partonflavour3'].Integral(60, 120)
    n4 = histos['mjj_partonflavour4'].Integral(60, 120)
    n5 = histos['mjj_partonflavour5'].Integral(60, 120)
    n = np.array([n1,n2,n3,n4,n5])
    total = sum(n)

    frac1 = n1/total # Expected value : 0.115
    frac2 = n2/total # Expected value : 0.156
    frac3 = n3/total # Expected value : 0.156
    frac4 = n4/total # Expected value : 0.115
    frac5 = n5/total # Expected value : 0.151
    BR_Obs = np.array([frac1,frac2,frac3,frac4,frac5])
    BR_Theo = np.array([0.156, 0.116, 0.156, 0.116, 0.156])  # d, u, s, c, b, t
    BR_Theo = BR_Theo * (1/np.sum(BR_Theo))
    sigma = np.sqrt(n)/total 
    ndof = len(BR_Obs) - 1
    difference = BR_Obs-BR_Theo
    chi2 = np.sum((BR_Obs-BR_Theo)**2/BR_Theo)
    p = 1 - stats.chi2.cdf(chi2,ndof)       
 
    return df, histos, BR_Obs, n, total, ndof, chi2, difference, p
