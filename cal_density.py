def cal_density(adata, mapping='X_tsne'):
    import numpy as np
    from scipy.stats import gaussian_kde

    x = adata.obsm[mapping][:,0]
    y = adata.obsm[mapping][:,1]


    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    dens_name = mapping[2:] + '_density'
    adata.obs[dens_name]=z

    return(adata)


def plot_density(adata,adata_full, mapping='X_tsne', max_scale=0.15, save=False):
    import matplotlib.pyplot as pl

    from scipy.stats import gaussian_kde
    dens_name = mapping[2:] + '_density'
    z = adata.obs[dens_name]
    x = adata.obsm[mapping][:,0]
    y = adata.obsm[mapping][:,1]
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    fig, ax = pl.subplots(frameon=False)
    im = ax.scatter(adata_full.obsm[mapping][:,0],
               adata_full.obsm[mapping][:,1], c='grey', edgecolor='', rasterized=True)
    im = ax.scatter(x, y, c=z, s=10, edgecolor='',cmap='Spectral_r', rasterized=True)
    cbar = fig.colorbar(im)
    #pl.imshow(ax)
    cbar.set_clim(0.,max_scale)
    pl.grid(b=False)
  
    if (save!=False):
        pl.savefig(save, dpi=300)
    pl.show()
    return()
