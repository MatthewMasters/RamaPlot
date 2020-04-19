import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

def read_pdb(input_file):
	# Parse PDB file into pandas dataframe
	data = []
	fi = open(input_file)
	for l in fi:
		if l[:4] == 'ATOM': data.append([l[12:16].strip(), l[21:22].strip(), l[17:20].strip(),
				                         int(l[22:26]), float(l[30:38]), float(l[38:46]), float(l[46:54])])
	fi.close()
	df = pd.DataFrame(data, columns=['atomName', 'chain', 'resName', 'resId', 'x', 'y', 'z'])
	return df

def compute_angle(p):
	# Compute dihedral angle provided a 4x3 matrix of coordinates
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array([v - (v.dot(b[1]) / b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    x = np.dot(v[0], v[1])
    y = np.dot(np.cross(v[0], b[1] / np.linalg.norm(b[1])), v[1])
    return np.degrees(np.arctan2(y, x))

def get_dihedrals(df):
	# Compute dihedral angles for each residue in PDB
	angles = []
	for chain in pd.unique(df.chain):
		chain_df = df[df['chain'] == chain]
		for resId in pd.unique(chain_df.resId):
			# Get residues i-1, i, i+1
			prevRes = chain_df[(chain_df['resId'] == resId - 1)]
			currRes = chain_df[chain_df['resId'] == resId]
			nextRes = chain_df[(chain_df['resId'] == resId + 1)]
			
			# Get position of atoms required to compute dihedrals
			positions = []
			for res, name in [[prevRes, 'C'], [currRes, 'N'], [currRes, 'CA'], [currRes, 'C'], [nextRes, 'N']]:
				q = res[res['atomName'] == name]
				if len(q) == 0: break
				positions.append(q[['x','y','z']].values)
			if len(positions) != 5: continue
			positions = np.vstack(positions)

			# Compute dihedrals
			phi = compute_angle(positions[:-1])
			psi = compute_angle(positions[1:])
			angles.append([phi, psi])
	return np.array(angles)

def get_ax():
	# Returns formatted axis
	fig, ax = plt.subplots(figsize=(6,6))
	ax.set_xlim([-180, 180])
	ax.set_ylim([-180, 180])
	ax.set_xlabel('Phi')
	ax.set_ylabel('Psi')
	return ax

def get_heatmap(dihedrals):
	# Computes heatmap by binning points
	bins = (dihedrals).astype(np.int) + 180
	heatmap = np.zeros((360, 360))
	for i,j in bins: heatmap[i,j] += 1
	heatmap = gaussian_filter(heatmap, sigma=8)
	return heatmap

def plot_points(dihedrals):
	# Plots Ramachandran plot with points representing each residue
	phis, psis = dihedrals[:,0], dihedrals[:,1]
	ax = get_ax()
	ax.plot([-180,180],[0,0], c='gray')
	ax.plot([0,0],[-180,180], c='gray')
	ax.scatter(phis, psis, zorder=10)
	plt.show()

def plot_heatmap(dihedrals):
	# Plots Ramachandran plot with heatmap of likelyhood
	heatmap = get_heatmap(dihedrals)[:-1, :-1]
	vmax = np.abs(heatmap).max()/8
	y, x = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-180, 180, 360))
	ax = get_ax()
	ax.pcolormesh(x, y, heatmap, cmap='jet', vmin=0, vmax=vmax)
	plt.show()

def plot_contour(dihedrals):
	# Plots Ramachandran plot with contour-plot of likelyhood
	heatmap = get_heatmap(dihedrals)
	vmax = np.abs(heatmap).max()
	y, x = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-180, 180, 360))
	ax = get_ax()
	ax.contour(x, y, heatmap, cmap='jet', vmin=0, vmax=vmax)
	plt.show()

if __name__ == '__main__':
	input_file = '1atp.pdb'
	pdb_df = read_pdb(input_file)
	dihedrals = get_dihedrals(pdb_df)
	plot_points(dihedrals)
	plot_heatmap(dihedrals)
	plot_contour(dihedrals)