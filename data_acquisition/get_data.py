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
        self.mcmc_data = None
        specs,meta = Data.Spectra(sim='lgal',noise = 'none', lib = 'bc03', sample = 'mini_mocha')
        photo, _ = Data.Photometry(sim='lgal', noise= 'none', lib='bc03', sample = 'mini_mocha') 
        self.meta_data = meta
        self.prior = self.f['prior_range']
        
        
    def plot_raw(self,param_idx,all_data = False):
        keys = list(self.f.keys())
        # Get the number of mcmc_chain columns
        num = -1
        for k in keys:
            if 'mcmc_chain' in k:
                num += 1
        
        # If all_data == False, take the last mcmc_chain column
        if not all_data:
            self.mcmc_data = self.f[f'mcmc_chain{num}'][...][:,:,param_idx]

        # Otherwise, take all 
        else:
            self.mcmc_data = np.array([])
            for idx in range(num+1):
                self.mcmc_data = np.append(self.mcmc_data,self.f[f'mcmc_chain{idx}'][...][:,:,param_idx])
        
        fig, ax = plt.subplots(1,1,figsize = (10,10))
        length = len(self.mcmc_data)//self.num_walkers
        for i in range(self.num_walkers):
            plt.plot(np.arange(length),self.mcmc_data[i::self.num_walkers], c = 'steelblue', lw = 0.5)
        plt.ylim(self.prior[param_idx])
        label = np.array(self.f['theta_names'][...]).astype(str)[param_idx]
        plt.ylabel(f'Prior Range\n{label})')
        plt.title(f'igal{self.igal} Walker Distribution {label}')
        plt.xlabel('Iteration')
        
    def plot_total_median(self,param_idx, inc = 1000, thin = 0):
        keys = list(self.f.keys())
        
        num = -1
        for k in keys:
            if 'mcmc_chain' in k:
                num += 1
                
        self.mcmc_data = np.array([])
        for idx in range(num+1):
            self.mcmc_data = np.append(self.mcmc_data,self.f[f'mcmc_chain{idx}'][...][:,:,param_idx])
        
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
        fig, ax = plt.subplots(1,1,figsize = (10,10))
        plt.plot(np.arange(len(total_med))*inc*self.num_walkers, total_med, c = 'steelblue', lw = 0.5)
        plt.ylim(self.prior[param_idx])
        label = np.array(self.f['theta_names'][...]).astype(str)[param_idx]
        plt.ylabel(f'Prior Range\n{label})')
        plt.title(f'igal{self.igal} Median {label}')
        plt.xlabel('Iteration')
        
    def plot_walker_median(self,param_idx,inc = 1000, thin = 0):
        keys = list(self.f.keys())
        
        num = -1
        for k in keys:
            if 'mcmc_chain' in k:
                num += 1
                
        self.mcmc_data = np.array([])
        for idx in range(num+1):
            self.mcmc_data = np.append(self.mcmc_data,self.f[f'mcmc_chain{idx}'][...][:,:,param_idx])
        
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
        
        
        fig, ax = plt.subplots(1,1,figsize = (10,10))
        for walker in walker_mcmc:
            ax.plot(np.arange(len(walker))*inc, walker, c = 'steelblue', lw = 0.5)
        plt.ylim(self.prior[param_idx])
        label = np.array(self.f['theta_names'][...]).astype(str)[param_idx]
        plt.ylabel(f'Prior Range\n{label})')
        plt.title(f'igal{self.igal} Walker Median {label}')
        plt.xlabel('Iteration')
            
    def report_median(self,param_idx,all_data = False):
        keys = list(self.f.keys())
        
        num = -1
        for k in keys:
            if 'mcmc_chain' in k:
                num += 1
        if not all_data:
            self.mcmc_data = self.f[f'mcmc_chain{num}'][...][:,:,param_idx]
        else:
            self.mcmc_data = np.array([])
            for idx in range(num+1):
                self.mcmc_data = np.append(self.mcmc_data,self.f[f'mcmc_chain{idx}'][...][:,:,param_idx])
        med = np.median(self.mcmc_data)
        return med

if __name__ == '__main__':
	data_path = sys.argv[1]
	num_walker = int(sys.argv[2])
	sim = sys.argv[3]
	spec_or_photo = sys.argv[4]
	noise = sys.argv[5]
	model = sys.argv[6]
	
	param_medians = {}
	
	sample_f_name = f'{sim}.{spec_or_photo}.noise_{noise}.{model}.0.mcmc.hdf5'
	sample_path = os.path.join(data_path,sample_f_name)
	f_sample = h5py.File(sample_path,'r')
	param_list = np.array(f_sample['theta_names'][...]).astype(str)
    f_sample.close()

	for n, param in enumerate(param_list):
		param_med = []
		for i in range(97):
			try:
				obj = plotter(num_walker,f'{data_path}/{sim}.{spec_or_photo}.noise_{noise}.{model}.{i}.mcmc.hdf5',i)
				med = obj.report_median(n,True)
				param_med.append(med)
			except:
				print(f'igal{i} file not found')
				param_med.append('N/A')
		param_medians[param] = param_med

	sfr_param_medians = {}

	sample_sfr_f_name = f'{sim}.{spec_or_photo}.noise_{noise}.{model}.0.postproc.hdf5'
	sfr_path = os.path.join(data_path,sample_sfr_f_name)
	f_sample_sfr = h5py.File(sfr_path,'r')
	sfr_param_list = np.array(f_sample_sfr['theta_names'][...]).astype(str)
    f_sample_sfr.close()

	for n, param in enumerate(sfr_param_list):
		sfr_param_med = []
		for i in range(97):
			try:
				obj = plotter(num_walker,f'{data_path}/{sim}.{spec_or_photo}.noise_{noise}.{model}.{i}.postproc.hdf5',i)
				med = obj.report_median(n,True)
				sfr_param_med.append(med)
			except:
				print(f'igal{i} file not found')
				sfr_param_med.append('N/A')
		sfr_param_medians[param] = sfr_param_med

	save_path = '/global/homes/k/kgb0255/packages/gqp_mc/doc/data_list'
	f_mcmc_name = f'{sim}.{spec_or_photo}.noise_{noise}.{model}.mcmc.npy'
	f_sfr_name = f'{sim}.{spec_or_photo}.noise_{noise}.{model}.sfr.npy'

	np.save(os.path.join(save_path,f_mcmc_name),param_medians)
	np.save(os.path.join(save_path,f_sfr_name),sfr_param_medians)