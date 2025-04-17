# base
import math
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
from sklearn.metrics.pairwise import haversine_distances
from scipy.spatial.distance import pdist, squareform
import statsmodels.api as sm
import pickle

# viz
import matplotlib.pyplot as plt
from matplotlib import gridspec
import cartopy.crs as ccrs

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz, Objective, query_node_attributes
from feems.cross_validation import run_cv, comp_mats

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"

def cov_to_dist(S):
    s2 = np.diag(S).reshape(-1, 1)
    ones = np.ones((s2.shape[0], 1))
    D = s2 @ ones.T + ones @ s2.T - 2 * S
    return(D)

data_path = "/mendel-nas1/dhooper/feems"

# read the genotype data and mean impute missing data
(bim, fam, G) = read_plink("{}/220418.BTFs.miss65.mm80.feems".format(data_path))
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)
n, p = genotypes.shape

print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))

# setup graph
coord = np.loadtxt("{}/feems_BTFs.miss65.mm80.samples.coord".format(data_path))  # sample coordinates
outer = np.loadtxt("{}/feems_BTFs.range.outer".format(data_path))  # outer coordinates
grid_path = "{}/world_triangle_res8.shp".format(data_path)  # path to discrete global grid

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=False, 
                                             buffer=0,
                                             outer=outer)

sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

projection = ccrs.EquidistantConic(central_longitude=143.5, central_latitude=-18.3)

fig = plt.figure(dpi=750)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=10, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=5)
v.draw_map()
v.draw_samples()
v.draw_edges(use_weights=False)
v.draw_obs_nodes(use_ids=False)

plt.savefig('feems_220418.BTFs.pops.res8.png', bbox_inches='tight')

#Fit feems with cross-validation
LoadCVFromDisk = False

# define grid
# Publication grid for this dataset
lamb_grid = np.geomspace(1e-6, 1e2, 20)[::-1]
# Exploratory grid: 
#lamb_grid = np.geomspace(1e-4, 1e1, 10)[::-1]

# run cross-validation
if not LoadCVFromDisk:
    cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)
    pickle.dump(cv_err,open("cv_err.pkl","wb"))

#Plot results of cross-validation
LoadCVFromDisk = True
if LoadCVFromDisk: 
    cv_err = pickle.load(open("cv_err.pkl","rb"))
    
# average over folds
mean_cv_err = np.mean(cv_err, axis=0)

# argmin of cv error
lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)])

fig, ax = plt.subplots(dpi=500)
ax.plot(np.log10(lamb_grid), mean_cv_err, ".");
ax.set_xlabel("log10(lambda)");
ax.set_ylabel("L2 CV Error");
ax.axvline(np.log10(lamb_cv), color = "orange")
lamb_cv

plt.savefig('feems_220418.BTFs.cross-validation.png', bbox_inches='tight')

#Re-fit feems using best lambda
# re-fit
sp_graph.fit(lamb_cv)

# Plot the FEEMS result
fig = plt.figure(dpi=1000)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=5)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
v.draw_edge_colorbar()

plt.savefig('feems_220418.BTFs.cross-validation_bestfit.png', bbox_inches='tight')

# Initialize figure
fig = plt.figure(constrained_layout=True, dpi=750, figsize=(8, 6))
spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)

# (A) Genetic distance vs geographic distance
D_geno = squareform(pdist(genotypes, metric="sqeuclidean")) / p
coord_rad = coord[:,::-1] * math.pi / 180.0
D_geo = haversine_distances(coord_rad) * 6371000/1000
tril_idx = np.tril_indices(n, k=-1)
x = D_geo[tril_idx]
y = D_geno[tril_idx]
X = sm.add_constant(x)
mod = sm.OLS(y, X)
res = mod.fit()
muhat, betahat = res.params

ax_00 = fig.add_subplot(spec[0, 0])
ax_00.set_title("A", loc='left')
ax_00.scatter(x, 
              y, 
              marker=".", 
              alpha=1, 
              zorder=0, 
              color="black",
              s=3)

x_ = np.linspace(np.min(x), np.max(x), 20)
ax_00.plot(x_, muhat + betahat * x_, zorder=2, color="orange", linestyle='--', linewidth=1)
ax_00.text(3500, .6, "R²={:.4f}".format(res.rsquared))
ax_00.set_xlabel("great circle distance (km)")
ax_00.set_ylabel("genetic distance")

# (B) Genetic distance vs fitted distance for constant w model 
tril_idx = np.tril_indices(sp_graph.n_observed_nodes, k=-1)
ax_01 = fig.add_subplot(spec[0, 1])
ax_01.set_title("B", loc='left')
sp_graph.fit_null_model()
sp_graph.comp_graph_laplacian(sp_graph.w)

obj = Objective(sp_graph)
fit_cov, _, emp_cov = comp_mats(obj)
fit_dist = cov_to_dist(fit_cov)[tril_idx]
emp_dist = cov_to_dist(emp_cov)[tril_idx]
X = sm.add_constant(fit_dist)
mod = sm.OLS(emp_dist, X)
res = mod.fit()
muhat, betahat = res.params
ax_01.scatter(fit_dist, 
              emp_dist, 
              marker=".", 
              alpha=1, 
              zorder=0, 
              color="black",
              s=3)

x_ = np.linspace(np.min(fit_dist), np.max(fit_dist), 20)
ax_01.plot(x_, muhat + betahat * x_, zorder=2, color="orange", linestyle='--', linewidth=1)
ax_01.text(4.2, 2, "R²={:.4f}".format(res.rsquared))
ax_01.set_xlabel("fitted distance (constant w)")

# (C) Genetic distance vs fitted distance for lambda = lambda_cv
tril_idx = np.tril_indices(sp_graph.n_observed_nodes, k=-1)
ax_10 = fig.add_subplot(spec[1, 0])
ax_10.set_title("C", loc='left')
lamb = lamb_cv
sp_graph.fit(lamb=lamb,
             lb=math.log(1e-6), 
             ub=math.log(1e+6))
sp_graph.comp_graph_laplacian(sp_graph.w)

obj = Objective(sp_graph)
fit_cov, _, emp_cov = comp_mats(obj)
fit_dist = cov_to_dist(fit_cov)[tril_idx]
emp_dist = cov_to_dist(emp_cov)[tril_idx]
X = sm.add_constant(fit_dist)
mod = sm.OLS(emp_dist, X)
res = mod.fit()
muhat, betahat = res.params
ax_10.scatter(fit_dist,
              emp_dist,
              marker=".", 
              alpha=1, 
              zorder=0, 
              color="black",
              s=3)

x_ = np.linspace(np.min(fit_dist), np.max(fit_dist), 20)
ax_10.plot(x_, muhat + betahat * x_, zorder=2, color="orange", linestyle='--', linewidth=1)
ax_10.text(4.2, 2, "R²={:.4f}".format(res.rsquared))
ax_10.set_xlabel("fitted distance ($\lambda = \lambda_{CV}$)")
ax_10.set_ylabel("genetic distance")

# (D) Genetic distance vs fitted distance for lambda = 1e-3*lambda_cv
tril_idx = np.tril_indices(sp_graph.n_observed_nodes, k=-1)
ax_11 = fig.add_subplot(spec[1, 1])
ax_11.set_title("D", loc='left')
lamb = 0.001*lamb_cv
sp_graph.fit(lamb=lamb,
             lb=math.log(1e-6), 
             ub=math.log(1e+6))
sp_graph.comp_graph_laplacian(sp_graph.w)

obj = Objective(sp_graph)
fit_cov, _, emp_cov = comp_mats(obj)
fit_dist = cov_to_dist(fit_cov)[tril_idx]
emp_dist = cov_to_dist(emp_cov)[tril_idx]
X = sm.add_constant(fit_dist)
mod = sm.OLS(emp_dist, X)
res = mod.fit()
muhat, betahat = res.params
ax_11.scatter(fit_dist,
              emp_dist,
              marker=".", 
              alpha=1, 
              zorder=0, 
              color="black",
              s=3)

x_ = np.linspace(np.min(fit_dist), np.max(fit_dist), 20)
ax_11.plot(x_, muhat + betahat * x_, zorder=2, color="orange", linestyle='--', linewidth=1)
ax_11.text(4.2, 2, "R²={:.4f}".format(res.rsquared))
ax_11.set_xlabel("fitted distance ($\lambda = \lambda_{CV}\cdot 10^{-3}$)")
ax_11.set_ylabel("genetic distance")

plt.savefig('feems_kelsie.BTFs.model_comparison.res8.png', bbox_inches='tight')

#INSPECT THE OBSERVED VS FITTED COVARIANCES ON SAMPLED DEMES
K = 20

# Plot Observed vs Fitted Covariances and residuals
obj = Objective(sp_graph)
fit_cov, _, emp_cov = comp_mats(obj)
# Set diagonals to zero to focus on covariances
np.fill_diagonal(fit_cov,0)
np.fill_diagonal(emp_cov,0)

# Plot matrices
fig = plt.figure(figsize = [15,5], dpi=1000)
ax1 = fig.add_subplot(141)
ax1.imshow(emp_cov, interpolation='nearest')
ax1.set_xlabel("empirical covariance")
ax2 = fig.add_subplot(142)
ax2.imshow(fit_cov, interpolation='nearest')
ax2.set_xlabel("fitted covariance")
ax3 = fig.add_subplot(143)
ax3.imshow(emp_cov-fit_cov, interpolation='nearest')
ax3.set_xlabel("residuals")
ax4 = fig.add_subplot(144)
ax4.hist(emp_cov.flatten())
ax4.set_xlabel("empirical covariance")
print("Min empirical covariance {:.2f}".format(np.min(emp_cov)))
print("25th percentile empirical covariance {:.2f}".format(np.percentile(emp_cov,0.25)))
print("75th percentile empirical covariance {:.2f}".format(np.percentile(emp_cov,0.75)))
print("Max empirical covariance {:.2f}".format(np.max(emp_cov)))

plt.savefig('feems_220418.BTFs.fit-evaluation1.res8.png', bbox_inches='tight')

# Plot the map with top K poorly fit covariances highlighted
fig = plt.figure(dpi=1000)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=5)
v.draw_map()
v.draw_obs_nodes(use_ids=False) 
v.draw_edges(use_weights=True)
#v.draw_edge_colorbar()

# Compute matrix of residuals and order them
res_mat = emp_cov-fit_cov
iu1 = np.triu_indices(res_mat.shape[0])
res_mat[iu1] = 0
abs_res_mat = abs(res_mat)

res_ord=np.argsort(abs_res_mat,axis=None)

permuted_idx = query_node_attributes(v.sp_graph, "permuted_idx")
obs_perm_ids = permuted_idx[: v.sp_graph.n_observed_nodes]
obs_grid = v.grid[obs_perm_ids, :]
for i in range(1,K):
    pair = np.unravel_index(res_ord[-i], res_mat.shape)

    v.ax.plot(
        [obs_grid[pair[0], 0],obs_grid[pair[1],0]],
        [obs_grid[pair[0], 1],obs_grid[pair[1],1]],
        color = 'gray', linewidth = 0.5
    )
    v.ax.text(
        obs_grid[pair[0], 0],obs_grid[pair[0], 1], str(pair[0]),
        horizontalalignment="center", verticalalignment="bottom",
        size=v.obs_node_textsize*0.8, zorder=v.obs_node_zorder,
    )
    v.ax.text(
        obs_grid[pair[1], 0],obs_grid[pair[1], 1], str(pair[1]),
        horizontalalignment="center", verticalalignment="bottom",
        size=v.obs_node_textsize*0.8, zorder=v.obs_node_zorder,
    )
    v.ax.text(
        np.mean([obs_grid[pair[0], 0],obs_grid[pair[1], 0]]),
        np.mean([obs_grid[pair[0], 1],obs_grid[pair[1], 1]]),
        str("{:.1f}".format(res_mat[pair[0],pair[1]])),
        horizontalalignment="center",
        verticalalignment="center",
        size=v.obs_node_textsize*0.7
    )

plt.savefig('feems_220418.BTFs.fit-evaluation2.res8.png', bbox_inches='tight')
