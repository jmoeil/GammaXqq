import ROOT
import numpy as np
import matplotlib.pyplot as plt
import os

#cut_dir = 'less'

# --- Import ROOT files ---
file_sg = ROOT.TFile('../../Source.root')
file_bg = ROOT.TFile('../../Background.root')

# --- Kinematic Variables ---
variables = {
    r'$\Delta\eta$': {'hist_name': 'Jet_delta_eta', 'xrange': (0, 5),'cut_dir': 'less'},
    r'$\Delta p_T$': {'hist_name': 'Jet_delta_pT', 'xrange': (0, 250),'cut_dir': 'less'},
    r'$\Delta\Phi$': {'hist_name': 'Jet_delta_phi', 'xrange': (-np.pi, np.pi),'cut_dir': 'less'},
    r'$\frac{p_{T2}}{p_{T1}}$': {'hist_name': 'Jet_pT2pT1', 'xrange': (0, 1),'cut_dir': 'greater'},
    r'$\Delta R$': {'hist_name': 'Jet_delta_R', 'xrange': (0, 4),'cut_dir': 'less'}
}

# --- Make a plot ---
def plot_roc(FPR,TPR,labels,title,outpath,auc=None):
	'''
	Plot ROC curves with given FPR/TPR values.

	Parameters:
	- FPR		List of lists (one list per curve)
	- TPR		List of lists (one list per curve)
	- labels	List of legend labels
	- Title		Title for the plot
	- Outpath	Path to save the figure
	'''
	plt.figure(figsize=(8,6))
	for i in range(len(FPR)):
		auc_text = f' (AUC={auc[i]:.3f})' if auc else ''
		plt.plot(FPR[i],TPR[i],label=labels[i] + auc_text)

	plt.plot([0,1],[0,1],'k--',label='Random guess')
	plt.xlabel(r'False Positive Rate (FPR) - Background efficiency $\epsilon_{\mathrm{bkg}}$')
	plt.ylabel(r'True Positive Rate (TPR) - Signal efficiency $\epsilon_{\mathrm{sig}}$')
	plt.title(title)
	plt.legend(loc='lower right')
	plt.grid(True)
	plt.tight_layout()
	os.makedirs(os.path.dirname(outpath),exist_ok=True)
	plt.savefig(outpath,dpi=300)
	plt.close()

# --- Plot all ROC curves on the same figure ---
print('*** Making global ROC curves ***\n')

TPR_list,FPR_list,label_list,auc_list = [],[],[],[]
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
        if info['cut_dir'] == 'less':
            sig_pass = sg_hist.Integral(1, bin_thr)
            bg_pass = bg_hist.Integral(1, bin_thr)
        else:
            sig_pass = sg_hist.Integral(bin_thr, sg_hist.GetNbinsX())
            bg_pass = bg_hist.Integral(bin_thr, bg_hist.GetNbinsX())
	
        TPR.append(sig_pass)
        FPR.append(bg_pass)
    
    FPR,TPR = zip(*sorted(zip(FPR,TPR)))
    auc = np.trapz(TPR, FPR)
    TPR_list.append(TPR)
    FPR_list.append(FPR)
    label_list.append(f'{label} (AUC={auc:.3f})')
    auc_list.append(auc)
    print(f'ROC curve for {label} ready')

plot_roc(FPR_list,TPR_list,label_list,title='ROC curves for signal discrimination',outpath='plots/roc_all_variables.png')

# --- Kinematic ROC per flavour ---
flavours = [f'PartonFlavour{i}' if i != 0 else '' for i in range(0,6)]
print('*** Making ROC curves per flavour ***\n')
for label, info in variables.items():
	print(f'\n=== Working on {label} ===\n')
	plt.figure(figsize=(8, 6))  # One figure per variable
	TPR_list,FPR_list,label_list,auc_list = [],[],[],[]
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
			if info['cut_dir'] == 'less':
				sig_pass = sg_hist.Integral(1,bin_thr)
				bg_pass = bg_hist.Integral(1,bin_thr)
			else:
				sig_pass = sg_hist.Integral(bin_thr,sg_hist.GetNbinsX())
				bg_pass = bg_hist.Integral(bin_thr,bg_hist.GetNbinsX())

			TPR.append(sig_pass)
			FPR.append(bg_pass)
		FPR,TPR = zip(*sorted(zip(FPR,TPR)))
		auc = np.trapz(TPR,FPR)
		legend_label = 'all' if flavour=='' else flavour
		
		TPR_list.append(TPR)
		FPR_list.append(FPR)
		label_list.append(legend_label)
		auc_list.append(auc)
		print(f'ROC curve for {label} {legend_label} ready')

	# --- Plot styling ---
	plot_roc(TPR_list,FPR_list,label_list,title=f'ROC curves for {label} by jet flavour', outpath=f'plots/ROC_{info["hist_name"]}.png',auc=auc_list)