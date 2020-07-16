import emcee
import sys
import os
import numpy as np
import emcee
import h5py

#This might not track the fake convergence where a walker is stuck out of prior range.
#This computes the autocorrelation time slightly differently from the sampler. They are similar.

def auto_window(taus, c):
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1


def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i

def autocorr_func_1d(x, norm=True):
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two(len(x))

    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
    acf /= 4 * n

    # Optionally normalize
    if norm:
        acf /= acf[0]

    return acf

def autocorr_new(y, c=5.0):
    f = np.zeros(y.shape[1])
    for yy in y:
        f += autocorr_func_1d(yy)
    f /= len(y)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]

if __name__ == '__main__':
	data_path = sys.argv[1]
	sim = sys.argv[2]
	spec_or_photo = sys.argv[3]
	noise = sys.argv[4]
	model = sys.argv[5]
	factor = float(sys.argv[6])
	output = sys.argv[7] == 'True'

	if output:
		output = np.zeros(97)

	f_base = os.path.join(data_path,f'{sim}.{spec_or_photo}.noise_{noise}.{model}.')
	for i in range(97):
		f_name = f_base + f'{i}.mcmc.hdf5'
		if not os.path.exists(f_name):
			print(f'igal{i} not found')
			continue
		else:
			file = h5py.File(f_name,'r')
			i = -1
			for k in list(file.keys()):
				if 'mcmc_chain' in k:
					i+=1
			mcmc_chains = []
			for idx in range(i):
				mcmc_chains.append(file[f'mcmc_chain{idx}'][...])

			ndim = np.array(mcmc_chains).shape[3]
			mcmc_chains = np.array(mcmc_chains).reshape(-1,40,ndim)
			niter = mcmc_chains.shape[0]
			taus = []
			for dim in range(ndim):
				each_chain = mcmc_chains[:,:,dim]
				tau = autocorr_new(each_chain.T)
				taus.append(tau)

			converged = (np.all(taus > (np.ones_like(taus) * niter * factor)))
			if converged:
				output[i] = 2
			else:
				output[i] = 1

	np.save(f'{sim}.{spec_or_photo}.noise_{noise}.{model}.convergence.npy', output)

