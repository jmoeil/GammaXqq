import ROOT
import numpy as np
import matplotlib.pyplot as plt
import csv

# ---- Load Delta R from CSV ---
def load_DeltaR_CSV(filename):
	deltaR_vals = []
	entries = []
	with open(filename,newline='') as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			deltaR_vals.append(float(row['Delta R']))
			entries.append(float(row['Entries']))
	return np.array(deltaR_vals),np.array(entries)


# Load signal and background data

deltaR_sig,entries_sig = load_DeltaR_CSV('../Approx_DeltaR_zg.csv')
deltaR_bkg,entries_bkg = load_DeltaR_CSV('../Approx_DeltaR_gjets.csv')

# Normalize to get PDF

pdf_sig = entries_sig / np.sum(entries_sig)
pdf_bkg = entries_bkg / np.sum(entries_bkg)

# Compute likelihood

likelihood_vals = []
for p_sig,p_bkg in zip(pdf_sig,pdf_bkg):
	denom = p_sig + p_bkg
	likelihood = p_sig/denom if denom>0 else 0.5
	likelihood_vals.append(likelihood)
likelihood_vals = np.array(likelihood_vals)

# Plot
plt.figure(figsize=(8,6))
plt.plot(deltaR_sig, likelihood_vals, label=r'Likelihood $\mathcal{L}(\Delta R)$', color='black')
plt.xlabel(r'$\Delta R$')
plt.ylabel('Likelihood Discriminant')
plt.title(r'Likelihood Discriminant from Approx. $\Delta R$')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("Likelihood_deltaR.png", dpi=300)
plt.show()
print("Successfully created Likelihood_deltaR.png")


# --- Import ROOT files ---
file_sig = ROOT.TFile.Open('../../Source.root')
file_bkg = ROOT.TFile.Open('../../Background.root')

# --- Define variables to process ---

variables = {
	r'$\Delta\eta$' : {
		'sig_hist': 'Jet_delta_eta_zg',
		'bkg_hist': 'Jet_delta_eta_gjets',
		'xlabel': r'$\Delta\eta$',
		'xrange': (0,5)
	},
	r'$\Delta p_T$': {
		'sig_hist': 'Jet_delta_pT_zg',
		'bkg_hist': 'Jet_delta_pT_gjets',
		'xlabel': r'$\Delta p_T$',
		'xrange': (0,250)
	},
	r'$\Delta\Phi$': {   
		'sig_hist': 'Jet_delta_phi_zg',
		'bkg_hist': 'Jet_delta_phi_gjets',
		'xlabel': r'$\Delta\Phi$',
		'xrange': (0,np.pi)
	},
	r'$\frac{p_{T2}}{p_{T1}}$': {
		'sig_hist': 'Jet_pT1pT2_zg',
		'bkg_hist': 'Jet_pT1pT2_gjets',
		'xlabel': r'$\frac{p_{T2}}{p_{T1}}$',
		'xrange': (0,1)
	}
}

# --- Loop over variables ---

for varname,info in variables.items():
	# --- Reset plots ---
	plt.figure()

	# --- Get histograms ---
	sig_hist = file_sig.Get(info['sig_hist'])
	bkg_hist = file_bkg.Get(info['bkg_hist'])

	# --- Normalize to make proper PDFs ---
	sig_hist.Scale(1 / sig_hist.Integral())
	bkg_hist.Scale(1 / bkg_hist.Integral())

	# --- Define likelihood function ---
	def compute_likelihood(x):
		bin_sig = sig_hist.FindBin(x)
		bin_bkg = bkg_hist.FindBin(x)

		p_sig = sig_hist.GetBinContent(bin_sig)
		p_bkg = bkg_hist.GetBinContent(bin_bkg)

		denom = p_sig + p_bkg
		return p_sig/denom if denom > 0 else 0.5

	# --- Scan across x-axis and compute likelihood values ---
	nbins = sig_hist.GetNbinsX()
	x_min,x_max = info['xrange']
	x_vals = np.array([sig_hist.GetBinCenter(i) for i in range(1,nbins+1)])
	x_vals = x_vals[(x_vals >= x_min) & (x_vals <= x_max)]
	likelihoods = np.array([compute_likelihood(x) for x in x_vals])
	# --- Name of variables ---
	safename = info['sig_hist'].split('_zg')[0]
	# --- Plotting ---
	plt.plot(x_vals,likelihoods,label='Likelihood ratio')
	plt.xlabel(r'$\Delta\eta$')
	plt.ylabel('Likelihood Discriminant')
	plt.title(f'Likelihood Discriminant for {varname}')
	plt.grid(True)
	plt.legend()
	plt.tight_layout()
	plt.savefig(f'Likelihood_{safename}.png',dpi=300,bbox_inches='tight')
	print(f'Sucessfully created Likelihood_{safename}.png')
