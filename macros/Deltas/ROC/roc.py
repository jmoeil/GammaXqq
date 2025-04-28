import ROOT
import numpy as np
import matplotlib
matplotlib.use("Agg")  # avoid GUI with matplotlib
import matplotlib.pyplot as plt
import os

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

# --- Different cuts ---

cuts = {
    'nocut': '',
    'pt2pt1': '_cut_pt2pt1',
    'Delta_R': '_cut_pt2pt1_Delta_R'
}

cut_labels = {
    'nocut': 'No cut',
    'pt2pt1': r'$p_{T2}/p_{T1} > 0.02$',
    'Delta_R': r'$\Delta R < 3.95$'
}

# --- Stock global curves ---

all_curves = {}

# --- Make a plot ---
def plot_roc(FPR,TPR,labels,title,outpath,auc=None):
        '''
        Plot ROC curves with given FPR/TPR values.

        Parameters:
        - FPR           List of lists (one list per curve)
        - TPR           List of lists (one list per curve)
        - labels        List of legend labels
        - Title         Title for the plot
        - Outpath       Path to save the figure
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
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(outpath,dpi=300)
        plt.close()

# --- Plot all ROC curves on the same figure ---

for cut_name,cut_suffix in cuts.items():
        print(f'\n*** Making global ROC curves after cut {cut_labels[cut_name]} ***')
        TPR_list,FPR_list,label_list,auc_list = [],[],[],[]

        for label, info in variables.items():
            hist_base = info['hist_name']
            sg_hist = file_sg.Get(f"{hist_base}{cut_suffix}_zg")
            bg_hist = file_bg.Get(f"{hist_base}{cut_suffix}_gjets")

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
            FPR = list(FPR)
            TPR = list(TPR)
            auc = np.trapz(TPR, FPR)
            all_curves[(label,cut_name)] = (FPR,TPR,auc)
            TPR_list.append(TPR)
            FPR_list.append(FPR)
            label_list.append(f'{label} [{cut_labels[cut_name]}]')
            auc_list.append(auc)
            print(f'ROC curve for {label} ready')
        folder_name = 'No_cut' if cut_name == 'nocut' else cut_name
        plot_roc(FPR_list,TPR_list,label_list,title=f'ROC curves for signal discrimination - {cut_labels[cut_name]}',outpath=f'plots/{folder_name}/ROC_all_variables_{cut_name}.png',auc=auc_list)

flavours = [f'PartonFlavour{k}' for k in range(1,6)]
for cut_name,cut_suffix in cuts.items():
        # --- Kinematic ROC per flavour ---
        print(f'\n*** Making ROC curves per flavour after cut {cut_name}***\n')
        for label, info in variables.items():
                print(f'\n=== Working on {label} ===\n')
                plt.figure(figsize=(8, 6))  # One figure per variable
                TPR_list,FPR_list,label_list,auc_list = [],[],[],[]

                if (label,cut_name) in all_curves:
                        FPR_all,TPR_all,auc_all = all_curves[label,cut_name]
                        FPR_list.append(list(FPR_all))
                        TPR_list.append(TPR_all)
                        auc_list.append(auc_all)
                        label_list.append("all")
                        print(f'ROC curve for {label} ready')
                        
                for flavour in flavours:
                        #suffix = f'_{flavour}' if flavour else ''
                        suffix = f'_{flavour}{cut_suffix}'
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
                        legend_label = flavour
                        
                        TPR_list.append(TPR)
                        FPR_list.append(FPR)
                        label_list.append(f'{legend_label} [{cut_labels[cut_name]}]')
                        auc_list.append(auc)
                        print(f'ROC curve for {label} {legend_label} ready')
                folder_name = 'No_cut' if cut_name == 'nocut' else cut_name
                cut_title = cut_labels[cut_name]
                plot_roc(FPR_list,TPR_list,label_list,title=f'ROC curves for {label} by jet flavour - {cut_title}', outpath=f'plots/{folder_name}/ROC_{info["hist_name"]}_{cut_name}.png',auc=auc_list)
