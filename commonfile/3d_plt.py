import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_function(data, clust, system):
  # Extract data dimensions
  x = data[:, 0]
  y = data[:, 1]
  z = data[:, 2]

  # Define plot ranges
  xmin, xmax = -6, 13
  ymin, ymax = -6, 13
  zmin, zmax = -6, 13

  # Calculate density
  hist, edges = np.histogramdd(r, bins = (5, 8, 4))
  hist, xedges, yedges, zedges = np.histogramdd(x, y, z, bins=(100,100,100), range=[[xmin, xmax], [ymin, ymax], [zmin, zmax]], density=True)

  # Apply log scale and normalize
  hist = -0.6 * np.log(hist)
  hist -= np.amin(hist)

  # Create 3D figure
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  # Plot density surface
  ax.plot_trisurf(xedges[1:], yedges[1:], zedges[1:], hist.T, cmap='viridis', alpha=0.8)

  # Set axis limits and labels
  ax.set_xlim3d(xmin + 2, xmax - 2)
  ax.set_ylim3d(ymin + 2, ymax - 2)
  ax.set_zlim3d(zmin + 2, zmax - 2)
  ax.set_xlabel('LD1', fontsize=16)
  ax.set_ylabel('LD2', fontsize=16)
  ax.set_zlabel('LD3', fontsize=16)

  # Add grid lines
  ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')

  # Plot cluster centroids and labels
  offset = ((-2.11, -2.43, 0), (1.18, 1.14, 0), (1.55, -1.31, 0), (-1.24, 0.94, 0), (0.58, -2.8, 0))
  for i in range(len(clust)):
    ax.scatter(clust[i, 0], clust[i, 1], clust[i, 2], color='red', s=100)
    ax.text(clust[i, 0] + offset[i][0], clust[i, 1] + offset[i][1], clust[i, 2] + offset[i][2], str(i+1), color='black', fontsize=16)

  # Create and customize colorbar
  cbar = fig.colorbar(hist, ax=ax)
  cbar.set_label('PMF (kcal/mol)', rotation=270, labelpad=18, fontsize=16)
  cbar.ax.tick_params(labelsize=16)

  # Save the figure
  fig.savefig('%s.pmf.png' % system)
  plt.close()

# Load data and cluster centers
clust = np.loadtxt("centers.dat")
data = np.loadtxt("lda.dat")

# Run the plot function
plot_function(data, clust, "tq_pbd_3d")

