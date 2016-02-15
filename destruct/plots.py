from __future__ import absolute_import
import itertools
import io
import tarfile
import scipy
import scipy.stats
import numpy as np


class gaussian_kde_set_covariance(scipy.stats.gaussian_kde):
    def __init__(self, dataset, covariance):
        self.covariance = covariance
        scipy.stats.gaussian_kde.__init__(self, dataset)
    def _compute_covariance(self):
        self.inv_cov = 1.0 / self.covariance
        self._norm_factor = np.sqrt(2*np.pi*self.covariance) * self.n


def filled_density(ax, data, c, a, xmin, xmax, cov):
    density = gaussian_kde_set_covariance(data, cov)
    xs = [xmin] + list(np.linspace(xmin, xmax, 2000)) + [xmax]
    ys = density(xs)
    ys[0] = 0.0
    ys[-1] = 0.0
    ax.plot(xs, ys, color=c, alpha=a)
    ax.fill(xs, ys, color=c, alpha=a)


def filled_density_weighted(ax, data, weights, c, a, xmim, xmax, cov):
    samples = np.array(list(itertools.chain(*[itertools.repeat(idx, cnt) for idx, cnt in enumerate(np.random.multinomial(10000, weights / weights.sum()))])))
    filled_density(ax, data[samples], c, a, xmim, xmax, cov)


def savefig_tar(tar, fig, filename):
    plot_buffer = io.BytesIO()
    fig.savefig(plot_buffer, format='pdf')
    info = tarfile.TarInfo(name=filename)
    info.size = plot_buffer.tell()
    plot_buffer.seek(0)
    tar.addfile(tarinfo=info, fileobj=plot_buffer)


