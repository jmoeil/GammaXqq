import ROOT
import numpy as np
import matplotlib.pyplot as plt

# --- Import ROOT files ---
file_sig = ROOT.TFile.Open('../../Source.root')
file_bkg = ROOT.TFile.Open('../../Background.root')

# If on new .root files, change "Jet_pT1pT2" to "Jet_pT2pT1"

# --- Base variable names ---
base_variables = {
    r'$\Delta\eta$': {'hist_name': 'Jet_delta_eta', 'xrange': (0, 5)},
    r'$\Delta p_T$': {'hist_name': 'Jet_delta_pT', 'xrange': (0, 250)},
    r'$\Delta\Phi$': {'hist_name': 'Jet_delta_phi', 'xrange': (-np.pi, np.pi)},
    r'$\frac{p_{T2}}{p_{T1}}$': {'hist_name': 'Jet_pT1pT2', 'xrange': (0, 1)},
    r'$\Delta R$': {'hist_name': 'Jet_delta_R', 'xrange': (0, 4)}
}

# --- Parton Flavours (1 to 5) ---
flavours = ['', 'PartonFlavour1', 'PartonFlavour2', 'PartonFlavour3', 'PartonFlavour4', 'PartonFlavour5']

# --- Process histograms ---
for varname, info in base_variables.items():
    for flavour in flavours:
        suffix = f'_{flavour}' if flavour else ''  # Empty string for main histograms
        
        sig_hist_name = f"{info['hist_name']}{suffix}_zg"
        bkg_hist_name = f"{info['hist_name']}{suffix}_gjets"
        
        # Get histograms
        sig_hist = file_sig.Get(sig_hist_name)
        bkg_hist = file_bkg.Get(bkg_hist_name)
        
        # Skip missing histograms
        if not sig_hist or not bkg_hist:
            print(f"Skipping {sig_hist_name} or {bkg_hist_name} (not found)")
            continue
        
        # Normalize histograms
        if sig_hist.Integral() > 0:
            sig_hist.Scale(1 / sig_hist.Integral())
        if bkg_hist.Integral() > 0:
            bkg_hist.Scale(1 / bkg_hist.Integral())
        
        # Define likelihood function
        def compute_likelihood(x):
            bin_sig = sig_hist.FindBin(x)
            bin_bkg = bkg_hist.FindBin(x)
            
            p_sig = sig_hist.GetBinContent(bin_sig)
            p_bkg = bkg_hist.GetBinContent(bin_bkg)
            
            denom = p_sig + p_bkg
            return p_sig / denom if denom > 0 else 0.5
        
        # Compute likelihood values
        nbins = sig_hist.GetNbinsX()
        x_min, x_max = info['xrange']
        x_vals = np.array([sig_hist.GetBinCenter(i) for i in range(1, nbins + 1)])
        x_vals = x_vals[(x_vals >= x_min) & (x_vals <= x_max)]
        likelihoods = np.array([compute_likelihood(x) for x in x_vals])
        
        # Plot
        plt.figure()
        plt.plot(x_vals, likelihoods, label='Likelihood ratio')
        plt.xlabel(varname)
        plt.ylabel('Likelihood Discriminant')
        plt.title(f'Likelihood Discriminant for {varname} {flavour}')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        
        # Save figure
        safename = f"{info['hist_name']}{suffix}"
        plt.savefig(f'Likelihood_{safename}.png', dpi=300, bbox_inches='tight')
        print(f'Successfully created Likelihood_{safename}.png')
        plt.close() # Close the figure to save memory
