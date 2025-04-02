import ROOT
import numpy as np
import matplotlib.pyplot as plt
import os

cut_dir = 'less'

# --- Import ROOT files ---
file_sg = ROOT.TFile('../../Source.root')
file_bg = ROOT.TFile('../../Background.root')

# --- Kinematic Variables ---
variables = {
    r'$\Delta\eta$': {'hist_name': 'Jet_delta_eta', 'xrange': (0, 5)},
    r'$\Delta p_T$': {'hist_name': 'Jet_delta_pT', 'xrange': (0, 250)},
    r'$\Delta\Phi$': {'hist_name': 'Jet_delta_phi', 'xrange': (-np.pi, np.pi)},
    r'$\frac{p_{T2}}{p_{T1}}$': {'hist_name': 'Jet_pT2pT1', 'xrange': (0, 1)},
    r'$\Delta R$': {'hist_name': 'Jet_delta_R', 'xrange': (0, 4)}
}

# --- Plot all ROC curves on the same figure ---
plt.figure(figsize=(8, 6))
print('*** Making global ROC curves ***\n')
for label, info in variables.items():
    sg_hist = file_sg.Get(f"{info['hist_name']}_zg")
    bg_hist = file_bg.Get(f"{info['hist_name']}_gjets")

    if not sg_hist or not bg_hist or sg_hist.Integral() == 0 or bg_hist.Integral() == 0:
        print(f"Skipping {label} due to missing or empty histograms.")
        continue

    # Normalize
    sg_hist.Scale(1 / sg_hist.Integral())
    bg_hist.Scale(1 / bg_hist.Integral())

    x_min, x_max = info['xrange']
    x_vals = np.linspace(x_min, x_max, 500)

    TPR, FPR = [], []
    for threshold in x_vals:
        bin_thr = sg_hist.GetXaxis().FindBin(threshold)
        if cut_dir == 'less':
            sig_pass = sg_hist.Integral(1, bin_thr)
            bg_pass = bg_hist.Integral(1, bin_thr)
        else:
            sig_pass = sg_hist.Integral(bin_thr, sg_hist.GetNbinsX())
            bg_pass = bg_hist.Integral(bin_thr, bg_hist.GetNbinsX())
	
        TPR.append(sig_pass)
        FPR.append(bg_pass)

    auc = np.trapz(TPR, FPR)
    plt.plot(FPR, TPR, label=f'{label} (AUC={auc:.3f})')
    print(f'ROC curve for {label} ready')


# --- Plot Styling ---
plt.plot([0,1], [0,1], 'k--', label='Random guess')
plt.xlabel(r'False Positive Rate (FPR) - Background efficiency $\epsilon_{background}$')
plt.ylabel(r'True Positive Rate (TPR) - Signal efficiency $\epsilon_{signal}$')
plt.title('ROC curves for signal discrimination')
plt.legend(loc='lower right')
plt.grid(True)
plt.tight_layout()
os.makedirs("plots", exist_ok=True)
plt.savefig(f"plots/roc_all_variables.png", dpi=300)
plt.close()

# --- Kinematic ROC per flavour ---

flavours = [f'PartonFlavour{i}' if i != 0 else '' for i in range(0,6)]
print('*** Making ROC curves per flavour ***\n')
for label, info in variables.items():
	print(f'\n=== Working on {label} ===\n')
	plt.figure(figsize=(8, 6))  # One figure per variable
	for flavour in flavours:
		suffix = f'_{flavour}' if flavour else ''
		
		sg_hist_name = f'{info["hist_name"]}{suffix}_zg'
		bg_hist_name = f'{info["hist_name"]}{suffix}_gjets'

		sg_hist = file_sg.Get(sg_hist_name)
		bg_hist = file_bg.Get(bg_hist_name)
		if not sg_hist or not bg_hist:
			print(f'Skipping {sg_hist_name} or {bg_hist_name} (not found)')
			continue

		if sg_hist.Integral() == 0 or bg_hist.Integral() == 0:
			print(f'Skipping empty histogram for {label} - {flavour or "all"}')
			continue
		sg_hist.Scale(1/sg_hist.Integral())
		bg_hist.Scale(1/bg_hist.Integral())

		x_min,x_max = info['xrange']
		x_vals = np.linspace(x_min,x_max,500)

		TPR,FPR = [],[]
		for threshold in x_vals:
			bin_thr = sg_hist.GetXaxis().FindBin(threshold)
			if cut_dir == 'less':
				sig_pass = sg_hist.Integral(1,bin_thr)
				bg_pass = bg_hist.Integral(1,bin_thr)
			else:
				sig_pass = sg_hist.Integral(bin_thr,sg_hist.GetNbinsX())
				bg_pass = bg_hist.Integral(bin_thr,bg_hist.GetNbinsX())

			TPR.append(sig_pass)
			FPR.append(bg_pass)
		auc = np.trapz(TPR,FPR)
		legend_label = 'all' if flavour=='' else flavour
		plt.plot(FPR,TPR,label=f'{legend_label} (AUC = {auc:.3f})')
		print(f'ROC curve for {label} {legend_label} ready')

	# --- Plot styling ---
	plt.plot([0,1],[0,1],'k--',label='Random guess')	
	plt.xlabel(r'False Positive Rate (FPR) - Background efficiency $\epsilon_{\mathrm{bkg}}$')
	plt.ylabel(r'True Positive Rate (TPR) - Signal efficiency $\epsilon_{\mathrm{sig}}$')
	plt.title(f'ROC curves for {label} by jet flavour')
	plt.legend(loc='lower right')
	plt.grid(True)
	plt.tight_layout()
	os.makedirs(f"plots/", exist_ok=True)
	plt.savefig(f"plots/ROC_{info['hist_name']}.png", dpi=300)
	plt.close()