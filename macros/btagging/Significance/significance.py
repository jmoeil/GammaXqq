import ROOT
import numpy as np
import matplotlib
matplotlib.use("Agg")  # avoid GUI with matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

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

def plot_significance_surface(FPR,TPR,thresholds_1,thresholds_2,Significance,outpath,zlabel=r'Significance $S/\sqrt{S+B}$'):
        '''
        Plot a 3D ROC surface:
        - x-axis:   Threshold for Jet_btagPNetB_1
        - y-axis:   Threshold for Jet_btagPNetB_2
        - z-axis:   TPR (signal efficiency)
        - colorbar: FPR (background efficiency)

        Parameters:
        - FPR   2D array of signal efficiencies
        - TPR   2D array of background efficiencies
        - thresholds_1  1D array of thresholds for Jet_btagPNetB_1
        - thresholds_2  1D array of thresholds for Jet_btagPNetB_2
        - outpath   Path to save figure
        - zlabel    Label for Z-axis
        '''
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111,projection='3d')

        x,y= np.meshgrid(thresholds_1,thresholds_2)

        surf = ax.plot_surface(x,y,Significance.T,cmap='viridis',edgecolor='k',alpha=0.85)
        ax.set_xlabel('Threshold Jet_btagPNetB_1')
        ax.set_ylabel('Threshold Jet_btagPNetB_2')
        ax.set_zlabel(zlabel)
        ax.set_title("3D Significance Surface")

        fig.colorbar(surf,shrink=0.5,aspect=10,label=zlabel)
        plt.tight_layout()
        plt.savefig(outpath,dpi=300)
        plt.close()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# Global parameters

file_sg = ROOT.TFile('../../Source.root')
file_bg = ROOT.TFile('../../Background.root')

sg_1 = file_sg.Get('Jet_btagPNetB_1_zg')
sg_2 = file_sg.Get('Jet_btagPNetB_2_zg')
bg_1 = file_bg.Get('Jet_btagPNetB_1_gjets')
bg_2 = file_bg.Get('Jet_btagPNetB_2_gjets')

assert sg_1 and sg_2 and bg_1 and bg_2, "Missing histograms!"

# Normalize histograms

#for h in [sg_1,sg_2]:
#    h.Scale(1/h.Integral())
#for h in [bg_1,bg_2]:
#    h.Scale(1/h.Integral())

# Threshold scan ranges

n_thresholds = 1000
thresholds_1 = np.linspace(0,1,n_thresholds)
thresholds_2 = np.linspace(0,1,n_thresholds)

S = np.zeros((n_thresholds,n_thresholds))
B = np.zeros((n_thresholds,n_thresholds))

# Fill ROC grid
for i,t1 in enumerate(thresholds_1):
      for j,t2 in enumerate(thresholds_2):
            bin_t1_sg = sg_1.GetXaxis().FindBin(t1)
            bin_t2_sg = sg_2.GetXaxis().FindBin(t2)
            bin_t1_bg = bg_1.GetXaxis().FindBin(t1)
            bin_t2_bg = bg_2.GetXaxis().FindBin(t2)

            sg_nbins_1 = sg_1.GetNbinsX()
            sg_nbins_2 = sg_2.GetNbinsX()
            sg_pass = min(sg_1.Integral(bin_t1_sg,sg_nbins_1),sg_2.Integral(bin_t2_sg,sg_nbins_2))

            bg_nbins_1 = bg_1.GetNbinsX()
            bg_nbins_2 = bg_2.GetNbinsX()
            bg_pass = min(bg_1.Integral(bin_t1_bg,bg_nbins_1),bg_2.Integral(bin_t2_bg,bg_nbins_2))

            S[i,j] = sg_pass
            B[i,j] = bg_pass

Significance = np.zeros_like(S)
mask = (S+B) > 0
Significance[mask] = S[mask]/np.sqrt(S[mask]+B[mask])

plot_significance_surface(S,B,thresholds_1,thresholds_2,Significance,outpath="3D_significance_scan.png")

flat = Significance.flatten()
sorted_indices = flat.argsort()[::-1] # Sorted from max to min

for idx in sorted_indices:
     i,j = np.unravel_index(idx, Significance.shape)
     if thresholds_1[i] < 0.05 or thresholds_2[j] < 0.05:
          continue
     i_max, j_max = i,j
     break

best_t1,best_t2 = thresholds_1[i_max], thresholds_2[j_max]

print(f'{best_t1} with {i_max} and {best_t2} with {j_max}')
print(f"S({i_max},{j_max}) = {S[i_max,j_max]:.1f}, B({i_max},{j_max}) = {B[i_max,j_max]:.1f}, Significance = {Significance[i_max,j_max]:.3f}")

thresholds = np.linspace(0, 1, 200)
significances = []
for t in thresholds:
    bin_sg = sg_1.GetXaxis().FindBin(t)
    bin_bg = bg_1.GetXaxis().FindBin(t)
    s = sg_1.Integral(bin_sg, sg_1.GetNbinsX())
    b = bg_1.Integral(bin_bg, bg_1.GetNbinsX())
    if s + b > 0:
        sig = s / np.sqrt(s + b)
    else:
        sig = 0
    significances.append(sig)

plt.plot(thresholds, significances)
plt.xlabel("Threshold Jet_btagPNetB_1")
plt.ylabel("Significance S/√(S+B)")
plt.grid(True)
plt.savefig("significance_vs_cut_1.png")

thresholds = np.linspace(0, 1, 200)
significances = []
for t in thresholds:
    bin_sg = sg_2.GetXaxis().FindBin(t)
    bin_bg = bg_2.GetXaxis().FindBin(t)
    s = sg_2.Integral(bin_sg, sg_2.GetNbinsX())
    b = bg_2.Integral(bin_bg, bg_2.GetNbinsX())
    if s + b > 0:
        sig = s / np.sqrt(s + b)
    else:
        sig = 0
    significances.append(sig)

plt.plot(thresholds, significances)
plt.xlabel("Threshold Jet_btagPNetB_1")
plt.ylabel("Significance S/√(S+B)")
plt.grid(True)
plt.savefig("significance_vs_cut_2.png")