import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob as glob
from matplotlib import cm
from matplotlib import colors as c
from gqp_mc import data as Data 
from gqp_mc import fitters as Fitters
import os
import sys

path = os.environ.get('GQPMC_DIR')
path = os.path.join(path,'mini_mocha','ispeculator','james')
print(path)

class plotter():
    def __init__(self,walkers,data_dir, igal):
        self.num_walkers = walkers
        self.igal = igal
        self.data_dir = data_dir
        self.f = h5py.File(self.data_dir,'r')
        self.prior = self.f['prior_range'][...]
        self.mcmc_data = None
        idx1 = data_dir.index('lgal.')
        idx2 = data_dir.index('hdf5')
        self.obj_name = data_dir[idx1:idx2+4]
        
    def plot_total_median(self,ax,param_idx, latest = True, inc = 1000, thin = 1):
        keys = list(self.f.keys())

        num = -1
        for k in keys:
            if 'mcmc_chain' in k:
                num += 1
        if not latest:
            self.mcmc_data = np.array([])
            for idx in range(num+1):
                self.mcmc_data = np.append(self.mcmc_data,self.f[f'mcmc_chain{idx}'][...][:,:,param_idx])
        else:
            self.mcmc_data = self.f[f'mcmc_chain{num}'][...][:,:,param_idx]
        
        walker_mcmc = []
        for i in range(self.num_walkers):
            walker_mcmc.append(self.mcmc_data[i::self.num_walkers])
            
        
        for i, walker in enumerate(walker_mcmc):
            running_median = []
            walker = walker[::thin]
            length = len(walker)
            for ii in range(length//inc):
                running_median.append(np.median(walker[:(ii+1)*inc:thin]))
            
            walker_mcmc[i] = running_median
        
        
        total_med = np.median(walker_mcmc,axis = 0)
        ax.plot(np.arange(len(total_med))*inc*self.num_walkers, total_med, c = 'steelblue', lw = 0.5)
        ax.set_ylim(self.prior[param_idx])
        label = np.array(self.f['theta_names'][...]).astype(str)[param_idx]
        ax.set_ylabel(f'Prior Range\n{label})')
        ax.set_title(self.obj_name+f'\n{label}')
        ax.set_xlabel('Iteration')

    def close(self):
        if self.f:
            self.f.close()
        else:
            pass

        return None

if __name__ == '__main__':
    data_path = sys.argv[1]
    num_walker = int(sys.argv[2])
    sim = sys.argv[3]
    spec_or_photo = sys.argv[4]
    noise = sys.argv[5]
    model = sys.argv[6]
    latest = sys.argv[7] == 'True'

    save_path = '/global/homes/k/kgb0255/packages/gqp_mc/doc/data_list/param_plots/'

    sample_f_name = f'{sim}.{spec_or_photo}.noise_{noise}.{model}.0.mcmc.hdf5'
    sample_path = os.path.join(data_path,sample_f_name)
    f_sample = h5py.File(sample_path,'r')
    param_list = np.array(f_sample['theta_names'][...]).astype(str)
    f_sample.close()
    n_param = len(param_list)

#    fig, axs = plt.subplot(97,n_param, figsize = (

    for i in range(97):
        try:
            fig, axs = plt.subplots(1,n_param,figsize=(5*n_param,5))
            obj = plotter(num_walker,f'{data_path}/{sim}.{spec_or_photo}.noise_{noise}.{model}.{i}.mcmc.hdf5',i)
            for n, param in enumerate(param_list):
                obj.plot_total_median(axs[n],n,latest)
                fig.savefig(f'{save_path}{sim}.{spec_or_photo}.noise_{noise}.{model}.{n}.param_plots.pdf', format = 'pdf')
                fig.savefig(f'{save_path}{sim}.{spec_or_photo}.noise_{noise}.{model}.{n}.param_plots.png', format = 'png')
            #   obj.close()
        except:
            pass

#    fig.savefig(f'{save_path}{sim}.{spec_or_photo}.noise_{noise}.{model}.param_plots.pdf', format = 'pdf')
#    fig.savefig(f'{save_path}{sim}.{spec_or_photo}.noise_{noise}.{model}.param_plots.png', format = 'png')


