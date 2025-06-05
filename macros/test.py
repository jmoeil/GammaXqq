import ROOT

# Enable batch mode (no GUI pop-up)
ROOT.gROOT.SetBatch(True)

# Open ROOT files
signal_file = ROOT.TFile.Open("Source.root")
background_file = ROOT.TFile.Open("Background.root")

# Retrieve histograms
signal_hist = signal_file.Get("mjj_cut_final_zg")
background_hist = background_file.Get("mjj_cut_final_gjets")

# Check if histograms were loaded correctly
if not signal_hist or not background_hist:
    print("Error: Could not load one or both histograms.")
    exit(1)

# Normalize histograms if needed (optional)
signal_hist.Scale(1.0 / signal_hist.Integral())
background_hist.Scale(1.0 / background_hist.Integral())

# Set styles
signal_hist.SetLineColor(ROOT.kRed)
signal_hist.SetLineWidth(2)
signal_hist.SetTitle("Invariant Mass Distribution; m_{jj} [GeV]; Normalized Entries")

background_hist.SetLineColor(ROOT.kBlue)
background_hist.SetLineWidth(2)

# Create canvas
canvas = ROOT.TCanvas("canvas", "mjj comparison", 800, 600)
signal_hist.Draw("HIST")
background_hist.Draw("HIST SAME")

# Add legend
legend = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
legend.AddEntry(signal_hist, "Signal (Z#gamma)", "l")
legend.AddEntry(background_hist, "Background (#gamma+jets)", "l")
legend.Draw()

# Save as JPG
canvas.SaveAs("mjj_comparison.jpg")

# Clean up
signal_file.Close()
background_file.Close()