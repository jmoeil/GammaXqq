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

'''
TO DO : uniformize mjj_partonflavour vs mjj_PartonFlavour
'''

# Global parameters

variables = {
        'Jet_delta_eta': (100,0,5),
        'Jet_delta_phi': (100,0,5),
        'Jet_delta_pT': (1000,0,1000),
        'Jet_pT2pT1': (100,0,1),
        'Jet_delta_R': (100,0,5)
    }
cut_conditions = {
        'pt2pt1': 'Jet_pT2pT1 > 0.02',
        'Delta_R': 'Jet_delta_R < 3.952'
}

# Key: Selected efficiency, correspoding cut values for (respectively) the leading and subleading jets and the cumulative histogram cut values.
efficiencies = {
        0.7: (0.733, 0.532, 0.629),
        0.8: (0.399, 0.234, 0.305),
        0.9: (0.11, 0.061, 0.077)
}

def defineWeight(df, isData):
    if isData:
        df = df.Define("unit_weight", "1.0")
        weight = "unit_weight"
    else:
        weight = "LHEWeight_originalXWGTUP"
    return df, weight

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

        histos[pt_hist_name] = df.Histo1D(ROOT.RDF.TH1DModel(pt_hist_name, '', 100, 0, 500), pt_col, weight)
        histos[eta_hist_name] = df.Histo1D(ROOT.RDF.TH1DModel(eta_hist_name, '', 50, -2.5, 2.5), eta_col, weight)
    
    return df

def cut_fill_histos(df,cut_expr,label, isData, weight,skipKinematics=0):
        '''
        Applies a given cut and fills histograms after each steps.

        Parameters:
- df             
- cut_expr
- label          

        returns:
- histos:               Updated with cumulative histograms

        '''
        histos_cut = {}
        df_cut = df.Filter(cut_expr,label)

        # Inclusive mjj
        histos_cut[f'mjj_cut_{label}'] = df_cut.Histo1D(
            ROOT.RDF.TH1DModel(f'mjj_cut_{label}', '', 1000, 0, 1000),
                'Mjj',
                weight
        )

        # Flavour-filtered mjj
        if not skipKinematics and not isData:
            for k in range(1,6):
                histos_cut[f'mjj_PartonFlavour_{k}_cut_{label}'] = df_cut.Filter(f'Sum(Jet_TightID_Pt30_Central_PartonFlavour{k})==2').Histo1D(ROOT.RDF.TH1DModel(f'mjj_PartonFlavour_{k}_cut_{label}', '', 1000, 0, 1000), 'Mjj', weight)

        # Inclusive Kinematic histograms
        if not skipKinematics:
            for varname,(nbins,xmin,xmax) in variables.items():
                histos_cut[f'{varname}_cut_{label}'] = df_cut.Histo1D(ROOT.RDF.TH1DModel(f'{varname}_cut_{label}','',nbins,xmin,xmax),varname,weight)

        # Flavour-filtered kinematic histograms
        if not skipKinematics and not isData:
            for varname,(nbins,xmin,xmax) in variables.items():
                for k in range(1,6):
                    flav_df_cut = df_cut.Filter(f'Sum(Jet_TightID_Pt30_Central_PartonFlavour{k})==2')
                    histos_cut[f'{varname}_PartonFlavour{k}_cut_{label}'] = flav_df_cut.Histo1D(ROOT.RDF.TH1DModel(f'{varname}_PartonFlavour{k}_cut_{label}','',nbins,xmin,xmax),varname,weight)

        return histos_cut

def apply_cumulative_cuts(df, cut_conditions, histos, isData,weight):
    """
    Applies a sequence of cuts cumulatively using cut_fill_histos.

    Parameters:
    - df: ROOT RDataFrame
    - cut_conditions: dict of {label: condition_expr} in the order of application
    - histos: histogram dictionary to be filled cumulatively

    Returns:
    - df: final filtered dataframe
    - histos: updated histogram dictionary
    """
    df_cut = df
    cumulative_labels = []

    for cut_label, cut_expr in cut_conditions.items():
        cumulative_labels.append(cut_label)
        label = "_".join(cumulative_labels)
        new_histos = cut_fill_histos(df_cut, cut_expr, label, isData,weight)
        histos.update(new_histos)
        df_cut = df_cut.Filter(cut_expr, cut_label)

    return df_cut, histos


def fill_btagPNetB(df, histos,isData, weight):
        """
        Fill inclusive and per-flavour histograms for Jet_btagPNetB of the two selected jets
        """

        # Define Jet_btagPNetB values for the two selected jets
        df = df.Define("Jet_TightID_Pt30_Central_btagPNetB",'Jet_btagPNetB[Jet_TightID_Pt30_Central]')
        df = df.Define("Jet_btagPNetB_1","Jet_TightID_Pt30_Central_btagPNetB[0]")
        df = df.Define("Jet_btagPNetB_2","Jet_TightID_Pt30_Central_btagPNetB[1]")
        df = df.Define("Jet_btagPNetB_mean","(Jet_btagPNetB_1 + Jet_btagPNetB_2)/2")

        histos["Jet_btagPNetB_1"] = df.Histo1D(ROOT.RDF.TH1DModel("Jet_btagPNetB_1",'',1000,0,1),"Jet_btagPNetB_1",weight)
        histos["Jet_btagPNetB_2"] = df.Histo1D(ROOT.RDF.TH1DModel("Jet_btagPNetB_2",'',1000,0,1),"Jet_btagPNetB_2",weight)
        histos["Jet_btagPNetB_mean"] = df.Histo1D(ROOT.RDF.TH1DModel("Jet_btagPNetB_mean", "", 1000, 0, 1),"Jet_btagPNetB_mean",weight)
        histos["Jet_btagPNetB"] = df.Histo2D(ROOT.RDF.TH2DModel("Jet_btagPNetB",'',1000,0,1,1000,0,1),"Jet_btagPNetB_1","Jet_btagPNetB_2",weight)
     
        if not isData:
            for k in range(1,6):
                filter_expression = f'Sum(Jet_TightID_Pt30_Central_PartonFlavour{k})==2'
                flav_df = df.Filter(filter_expression)
                histos[f"Jet_btagPNetB_PartonFlavour{k}_1"] = flav_df.Histo1D(ROOT.RDF.TH1DModel(f"Jet_btagPNetB_PartonFlavour{k}_1", "", 1000, 0, 1), "Jet_btagPNetB_1", weight)
                histos[f"Jet_btagPNetB_PartonFlavour{k}_2"] = flav_df.Histo1D(ROOT.RDF.TH1DModel(f"Jet_btagPNetB_PartonFlavour{k}_2","",1000,0,1),"Jet_btagPNetB_2",weight)
                histos[f"Jet_btagPNetB_PartonFlavour{k}_mean"] = flav_df.Histo1D(ROOT.RDF.TH1DModel(f"Jet_btagPNetB_PartonFlavour{k}_mean",'',1000,0,1),"Jet_TightID_Pt30_Central_btagPNetB",weight)
                histos[f"Jet_btagPNetB_PartonFlavour{k}"] = flav_df.Histo2D(ROOT.RDF.TH2DModel(f"Jet_btagPNetB_PartonFlavour{k}",'',1000,0,1,1000,0,1),f"Jet_btagPNetB_1",f"Jet_btagPNetB_2",weight)

        return df, histos

def GammaZSelection(df, year=2023, era='C', isData=False):
    '''
    Select events with = 1 photon with pT>100 GeV.
    The event must pass a single photon trigger. 
    Requires exactly 2 jets
    '''

    print("→ Number of events in file:", df.Count().GetValue())
    df, weight = defineWeight(df, isData)
    
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
    histos['photon_pt_aftertrigger'] = df.Histo1D(ROOT.RDF.TH1DModel('photon_pt_aftertrigger', '', 1000, 0, 1000), 'Photon_LooseID_Pt20_pt', weight)

    #Some simple event filters:
    df = df.Filter('Sum(Photon_LooseID_Pt20)==1','=1 loose photon with p_{T}>20 GeV')
    df = df.Filter('Sum(Photon_TightID_Pt100)==1','=1 tight photon with p_{T}>100 GeV')
    #Plot photon pt. Because of the two filters above, this distribution starts at 100 GeV.
    histos['photon_pt_afterptcut'] = df.Histo1D(ROOT.RDF.TH1DModel('photon_pt_afterptcut', '', 1000, 0, 1000), 'Photon_LooseID_Pt20_pt', weight)
    
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
    if not isData:
        df = df.Define('Jet_TightID_Pt30_Central_partonFlavour', 'Jet_partonFlavour[Jet_TightID_Pt30_Central]')
        for k in range(1,7):
            df = df.Define(f'Jet_TightID_Pt30_Central_PartonFlavour{k}',f'abs(Jet_TightID_Pt30_Central_partonFlavour)=={k}')
            
    #For now, consider only events with exactly 2 central jets and no other jet
    df = df.Filter('Sum(Jet_TightID_Pt30)==2&&Sum(Jet_TightID_Pt30_Central)==2','=2 central jets with pt>30 GeV, no additional jet')
    histos['photon_pt_2jselection'] = df.Histo1D(ROOT.RDF.TH1DModel('photon_pt_2jselection', '', 1000, 0, 1000), 'Photon_LooseID_Pt20_pt', weight)

    #Compute the dijet invariant mass and remove invariant mass mjj outside baseline
    df = df.Define('Mjj', 'InvariantMass(Jet_TightID_Pt30_Central_Pt[0], Jet_TightID_Pt30_Central_Eta[0], Jet_TightID_Pt30_Central_Phi[0], Jet_TightID_Pt30_Central_Mass[0], Jet_TightID_Pt30_Central_Pt[1], Jet_TightID_Pt30_Central_Eta[1], Jet_TightID_Pt30_Central_Phi[1], Jet_TightID_Pt30_Central_Mass[1])')
    df = df.Filter('Mjj < 200 && Mjj > 40','Invariant mass clearly outside the range of this study')    
    
    # --- Define key kinematic variables to study their behaviors ---
    # Delta eta
    df = df.Define('Jet_delta_eta','abs(Jet_TightID_Pt30_Central_Eta[0]-Jet_TightID_Pt30_Central_Eta[1])')

    # Delta phi
    df = df.Define('Jet_delta_phi','abs(acos(cos(Jet_TightID_Pt30_Central_Phi[0]-Jet_TightID_Pt30_Central_Phi[1])))')

    # Delta pT
    df = df.Define('Jet_pt1','Jet_TightID_Pt30_Central_Pt[0]')
    df = df.Define('Jet_pt2','Jet_TightID_Pt30_Central_Pt[1]')
    df = df.Define('Jet_delta_pT','abs(Jet_pt1-Jet_pt2)')

    # pT1/pT2
    df = df.Define('Jet_pT2pT1','Ratio_pt(Jet_pt1,Jet_pt2)')

    # Angular distance R
    df = df.Define('Jet_delta_R','sqrt(pow(Jet_delta_eta,2) + pow(Jet_delta_phi,2))')

    # Veto Delta R(photon,j) > 0.4 with j = leading and subleading jet
    df = df.Define('PJet_Delta_eta','abs(Photon_TightID_Pt100_eta[0]-Jet_TightID_Pt30_Central_Eta[0])') # Delta eta leading jet and photon
    df = df.Define('PJet_Delta_phi','abs(acos(cos(Photon_TightID_Pt100_phi[0]-Jet_TightID_Pt30_Central_Phi[0])))') # Delta phi leading jet and photon
    df = df.Define('PSubJet_Delta_eta','abs(Photon_TightID_Pt100_eta[0]-Jet_TightID_Pt30_Central_Eta[1])') # Delta eta subleading jet and photon
    df = df.Define('PSubJet_Delta_phi','abs(acos(cos(Photon_TightID_Pt100_phi[0]-Jet_TightID_Pt30_Central_Phi[1])))') # Delta phi subleading jet and photon
    
    df = df.Define('PJet_Delta_R','sqrt(pow(PJet_Delta_eta,2) + pow(PJet_Delta_phi,2))') # Delta R leading jet and photon
    df = df.Define('PSubJet_Delta_R','sqrt(pow(PSubJet_Delta_eta,2) + pow(PSubJet_Delta_phi,2))') # Delta R subleading jet and photon

    df = df.Filter('PJet_Delta_R > 0.4 && PSubJet_Delta_R > 0.4',f'Angular distance between the photon and both jets is > 0.4')

    # Selection cut to further select kinematic variables in range of study
    df_cut = df.Filter('Mjj > 40 && Mjj < 200', 'mjj in range of study')
    for varname, (nbins, xmin, xmax) in variables.items():
        hist_name = f'{varname}'
        histos[hist_name] = df_cut.Histo1D(ROOT.RDF.TH1DModel(hist_name,'',nbins,xmin,xmax),varname,weight)
    histos['Jet_delta_phi_vs_delta_eta'] = df_cut.Histo2D(ROOT.RDF.TH2DModel('Jet_delta_phi_vs_delta_eta','',100,0,5,100,0,5),'Jet_delta_phi','Jet_delta_eta',weight)

    df, histos = fill_btagPNetB(df_cut, histos,isData, weight)

    # Plot these variables in histograms
    if not isData:
        for k in range(1,6):
            # Create a filtered datafrom for each flavour
            flav_var = f'Jet_TightID_Pt30_Central_PartonFlavour{k}'
            df_flav = df_cut.Filter(f'Sum({flav_var})==2',f'Both jets have flavour {k}')
            # Loop over each variable
            for varname, (nbins,xmin,xmax) in variables.items():
                hist_name = f'{varname}_PartonFlavour{k}'
                histos[hist_name] = df_flav.Histo1D(ROOT.RDF.TH1DModel(hist_name,'',nbins,xmin,xmax),varname,weight)
    
    #Compute the dijet invariant mass  
    histos['Mjj'] = df.Histo1D(ROOT.RDF.TH1DModel('mjj', '', 1000, 0, 1000), 'Mjj', weight)
    if not isData:
        for k in range(1,7):
            histos[f'mjj_partonflavour{k}'] = df.Filter(f'Sum(Jet_TightID_Pt30_Central_PartonFlavour{k})==2').Histo1D(ROOT.RDF.TH1DModel(f'mjj_partonflavour{k}','',1000,0,1000),'Mjj',weight)

    # Compute and add histograms for the different cuts 
    df, histos = apply_cumulative_cuts(df, cut_conditions, histos, isData,weight)  
    for label,condition in cut_conditions.items():
        extra_histos = cut_fill_histos(df, condition,label,isData,weight)
        histos.update(extra_histos)

    if not isData:
        df = plot_jet_kinematics_by_flavour(df, histos)

    btag_cut = "Jet_btagPNetB_1 > 0.3 && Jet_btagPNetB_2 > 0.3"
    btag_histo = cut_fill_histos(df, btag_cut, 'final', isData,weight,skipKinematics=1)
    histos.update(btag_histo)


    # Count the number of events per flavour
    if not isData:
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
    else:
        BR_Obs = np.empty_like(np.array([0,0,0,0,0]))
        n = np.empty_like(BR_Obs)
        total = np.sum(n)
        ndof = total
        chi2 = ndof
        difference = ndof
        p = ndof
    return df, histos, BR_Obs, n, total, ndof, chi2, difference, p
