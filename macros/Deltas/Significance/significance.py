import ROOT
import matplotlib.pyplot as plt
import numpy as np

# --- Import ROOT files ---
file_sg = ROOT.TFile('../../Source.root')
file_bg = ROOT.TFile('../../Background.root')

# --- Kinematic Variables ---
variables = {
    r'$\Delta\eta$': {'hist_name': 'Jet_delta_eta', 'xrange': (0, 5),'cut_dir': 'less'},
    r'$\Delta p_T$': {'hist_name': 'Jet_delta_pT', 'xrange': (0, 250),'cut_dir': 'less'},
    r'$\Delta\Phi$': {'hist_name': 'Jet_delta_phi', 'xrange': (0, np.pi),'cut_dir': 'less'},
    r'$\frac{p_{T2}}{p_{T1}}$': {'hist_name': 'Jet_pT2pT1', 'xrange': (0, 1),'cut_dir': 'greater'},
    r'$\Delta R$': {'hist_name': 'Jet_delta_R', 'xrange': (0, 4),'cut_dir': 'less'}
}

# --- Plot S/sqrt(S+B) ---

# --- Make a plot ---
def plot_significance(cut_points, significances, variable_name=None, save_path=None):
    """
    Plot signal significance S / sqrt(S + B) vs cut threshold.

    Parameters:
    - cut_points: list or np.array of threshold values
    - significances: list or np.array of significance values
    - variable_name: (optional) string, name of the variable being cut on
    - save_path: (optional) path to save the figure instead of displaying it
    """
    plt.figure(figsize=(8, 5))

    # Find the max significance
    max_idx = np.argmax(significances)
    best_cut = cut_points[max_idx]
    best_significance = significances[max_idx]

    # Plot curve
    plt.plot(cut_points, significances, label='Signal significance', lw=2)
    plt.xlabel(f'Cut value on {variable_name}' if variable_name else 'Cut value')
    plt.axvline(best_cut,color='red',linestyle='--',label=f'Highest signal significance is {best_significance:.3f} at {best_cut:.3f}')
    plt.scatter([best_cut],[best_significance],color='red',zorder=5)
    plt.ylabel(r'Significance $S / \sqrt{S + B}$')
    title = f'Signal Significance vs Cut on {variable_name}' if variable_name else 'Signal Significance vs Cut'
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()

# --- Plot significance curves ---

for label, info in variables.items():
    sg_hist = file_sg.Get(f'{info["hist_name"]}_zg')
    bg_hist = file_bg.Get(f'{info["hist_name"]}_gjets')

    if not sg_hist or not bg_hist or sg_hist.Integral() == 0 or bg_hist.Integral() == 0:
        print(f"Skipping {label} due to missing or empty histograms.")
        continue

    # Normalize
    sg_hist.Scale(1 / sg_hist.Integral())
    bg_hist.Scale(1 / bg_hist.Integral())

    x_min, x_max = info['xrange']
    x_vals = np.linspace(x_min, x_max, 500)
    significances,cut_points = [],[]
    for threshold in x_vals:
        bin_thr = sg_hist.GetXaxis().FindBin(threshold)
        if info['cut_dir'] == "less":
            S = sg_hist.Integral(1, bin_thr)
            B = bg_hist.Integral(1, bin_thr)
        else:
            S = sg_hist.Integral(bin_thr, sg_hist.GetNbinsX())
            B = bg_hist.Integral(bin_thr, bg_hist.GetNbinsX())
        
        if S + B > 0:
            significance = S/np.sqrt(S+B)
            significances.append(significance)
            cut_points.append(threshold)
    plot_significance(cut_points,significances,label,f'plots/Significance_{info["hist_name"]}.png')
    print(f'Significance curve for {label} plotted.')