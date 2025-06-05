import ROOT
import numpy as np
from math import sqrt
import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend (no display needed)
import matplotlib.pyplot as plt

def extract_array_2D(hist):
    nbins_x = hist.GetNbinsX()
    nbins_y = hist.GetNbinsY()
    x_edges = np.array([hist.GetXaxis().GetBinLowEdge(i+1) for i in range(nbins_x + 1)])
    y_edges = np.array([hist.GetYaxis().GetBinLowEdge(i+1) for i in range(nbins_y + 1)])
    z_vals = np.array([[hist.GetBinContent(i+1,j+1) for i in range(nbins_x)] for j in range(nbins_y)])


   
    return x_edges,y_edges, z_vals

def plot_heatmap_2D(file,hist_name,outname):
    '''
    Plot a 2D histogram as a heatmap with contours.
    '''
    f = ROOT.TFile.Open(file)
    h = f.Get(hist_name)
    assert h, f"Histogram {hist_name} not found in {file}"

    # Binning and values
    x_edges, y_edges, z_vals = extract_array_2D(h)

    # Plot
    plt.figure(figsize=(8,6))
    im = plt.pcolormesh(x_edges,y_edges,z_vals,shading='auto',cmap='viridis')
    plt.contour(x_edges[:-1],y_edges[:-1], z_vals,levels=5,colors='white', linewidths=0.8)
    plt.xlabel("Jet_btagPNetB_1")
    plt.ylabel("Jet_btagPNetB_2")
    plt.colorbar(im,label="Entries")
    plt.title(f'{hist_name} heatmap')
    plt.tight_layout()
    plt.savefig(outname, dpi=300)
    plt.close()
    print(f'Histogram {hist_name} sucessfully saved as {outname}!')

def plot_significance_heatmap(hist_sg,hist_bg, outname):
    '''
    Plot a 2D significance heatmap from signal and background TH2F histograms.
    '''
    hist_sg.RebinX(10)
    hist_sg.RebinY(10)
    hist_bg.RebinX(10)
    hist_bg.RebinY(10)
    hist_significance = hist_sg.Clone()
    nbins_x = hist_sg.GetNbinsX()
    nbins_y = hist_sg.GetNbinsY()

    for i in range(nbins_x):
        print(i)
        for j in range(nbins_y):
            s = hist_sg.Integral(i+1, nbins_x+1, j+1, nbins_y+1) 
            b = hist_bg.Integral(i+1, nbins_x+1, j+1, nbins_y+1) 
            significance = 0
            if s > 0 and b >0 : 
                significance = s/sqrt(s+b)
            hist_significance.SetBinContent(i+1, j+1, significance)
    hist_significance.SaveAs("hist_significance.root")


# Call for signal and background
#plot_heatmap_2D("../Source.root","Jet_btagPNetB_PartonFlavour5_zg", "heatmap_signal.png")
#plot_heatmap_2D("../Background.root", "Jet_btagPNetB_gjets", "heatmap_background.png")

# Rebin to x,y to 5
# Set min to 0
# Set log z for background

file_sg = ROOT.TFile.Open("../Source.root")
file_bg = ROOT.TFile.Open('../Background.root')
hist_sg = file_sg.Get(("Jet_btagPNetB_PartonFlavour5_zg"))
hist_bg = file_bg.Get("Jet_btagPNetB_gjets")
assert hist_sg and hist_bg, "Missing histograms!"

plot_significance_heatmap(hist_sg, hist_bg, "significance_heatmap.png")