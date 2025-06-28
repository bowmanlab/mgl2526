# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 07:09:02 2025

@author: jeff bowman, jsbowman@ucsd.edu
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colormaps
import glob
import seawater
import itertools

bl_names = ['bottle', 'bottle_2', 'date.time', 'index_start', 'index_stop']

# cnv_names = ['depSM',
#              'CStarTr0',
#              'c0S/m',
#              'c1S/m',
#              'flECO-AFL',
#              'sbox0Mm/Kg',
#              'sbox1Mm/Kg',
#              'sal00',
#              'sal11',
#              't090C',
#              't190C',
#              'turbWETntu0',
#              'bpos',
#              'nbf',
#              'flag']

cnv_input_dir = "Z://public//CTD//"
bl_input_dir = "Z://public//CTD//raw_CTD//"

all_T = []
all_S = []
all_rho = []

bl_out = pd.DataFrame()
cast_out = pd.DataFrame()

for filename in glob.glob(cnv_input_dir + '*clean.csv'):
    os.remove(filename) 
    
def get_header_length(f):
    i = 0
    with open(cnv_input_dir + f, 'r') as f_file:
        for line in f_file.readlines():
            i = i + 1
            if "END" in line:
                return(i)
            
def get_var_names(f):
    var_names = []
    with open(cnv_input_dir + f, 'r') as f_file:
        for line in f_file.readlines():
            if line.startswith('# name'):
                var_name = line.strip().rstrip().split(' = ')[1]
                var_names.append(var_name)
        return(var_names)

files = os.listdir(cnv_input_dir)
files.sort()

for f in files:
    if f.endswith('.cnv'):
        print(f)
        skip = get_header_length(f)
        var_names = get_var_names(f)
        cnv_in = pd.read_csv(cnv_input_dir + f, header = None, sep = r'\s+', names = var_names, index_col = False, skiprows = skip)
        name = f.split('.cnv')[0]
        bl_in = pd.read_csv(bl_input_dir + name + '.bl', header = None, skiprows = 2, index_col = 0, names = bl_names)
        
        ## scrub the upcast by eliminating all points after max depth
        
        max_depth_i = cnv_in['depSM: Depth [salt water, m]'].idxmax()
        cnv_downcast = cnv_in.drop(index = range(max_depth_i, len(cnv_in.index)))
        
        ## remove any atmospheric measurements, including from PAR sensor at top of cage
        
        cnv_downcast = cnv_downcast[cnv_downcast['depSM: Depth [salt water, m]'] > 2]
        
        ## calculate density
        
        s = cnv_downcast['sal00: Salinity, Practical [PSU]']
        t = cnv_downcast['t090C: Temperature [ITS-90, deg C]']
        p = cnv_downcast['depSM: Depth [salt water, m]']
        cnv_downcast['rho [kg/m^3]'] = seawater.eos80.dens(s, t, p)
        
        ## % irradiance - this is commented out because PAR is PAR/Irradiance 
        
        #cnv_downcast['irradiance%'] = (cnv_downcast['par: PAR/Irradiance, Biospherical/Licor']/cnv_in['spar: SPAR, Biospherical/Licor']) * 100
        
        ## bin data by 2bdar
        
        cnv_downcast['bin'] = pd.cut(cnv_downcast['depSM: Depth [salt water, m]'], bins=300)
        cnv_binned = cnv_downcast.groupby('bin', observed = True).mean()
        
        s_binned = cnv_binned['sal00: Salinity, Practical [PSU]']
        t_binned = cnv_binned['t090C: Temperature [ITS-90, deg C]']
        p_binned = cnv_binned['depSM: Depth [salt water, m]']
            
        ## buoyancy frequency
                
        all_T = all_T + list(t)
        all_S = all_S + list(s)
        all_rho = all_rho + list(cnv_downcast['rho [kg/m^3]'])
        
        n2, q, p_ave = seawater.geostrophic.bfrq(s_binned, t_binned, p_binned)
        n2_flat = list(itertools.chain.from_iterable(n2))
        q_flat = list(itertools.chain.from_iterable(q))
        bvf = pd.DataFrame(zip(n2_flat, q_flat), index = [round(item, 4) for item in list(p_binned[1:])], columns = ['n2', 'q'])
        
        cast_out.loc[name, 'ml_depth'] = bvf.n2.idxmax()
        
        ## get bottle depth
        
        for bl in bl_in.index:
            start = bl_in.loc[bl, 'index_start']
            stop = bl_in.loc[bl, 'index_stop']
            z = cnv_in.loc[start:stop, 'depSM: Depth [salt water, m]'].mean().round(4) # depth
            z_i = abs(cnv_downcast['depSM: Depth [salt water, m]'] - z).idxmin()
            bl_out.loc[bl, 'depSM: Depth [salt water, m]'] = z
            bl_out.loc[bl, 'downcast_index'] = z_i
            
            ## populate bottle file
            
            bl_out.loc[bl, 'downcast_index'] = z_i
            
            for param in cnv_downcast.columns:
                bl_out.loc[bl, param] = cnv_downcast.loc[z_i, param]

        ## write out files
            
        cnv_downcast.to_csv(cnv_input_dir + name + '_clean.csv')
        bvf.to_csv(cnv_input_dir + name + '_bvf.csv')
            
        with PdfPages('Z://public//CTD//' + name + '_plots.pdf') as pdf:
            for param in cnv_downcast.columns:
                
                ## try clause is necessary because some parameters like "bin" can't be plotted
                
                try:
                    fig, ax = plt.subplots(figsize = (8.5, 11))
                    ax.plot(cnv_downcast[param], cnv_downcast['depSM: Depth [salt water, m]'])
                    ax.invert_yaxis()
                    ax.set_xlabel(param)
                    ax.set_ylabel('depSM: Depth [salt water, m]')
                    plt.title(name)
                    
                    if param == 'par: PAR/Irradiance, Biospherical/Licor':
                        temp_delta = abs(cnv_downcast[param] - 0.01)
                        n_i = temp_delta.idxmin()
                        comp_depth = cnv_downcast['depSM: Depth [salt water, m]'][n_i]
                        max_par = cnv_downcast[param].max()
                        plt.plot([0, max_par], [comp_depth, comp_depth])
                        cast_out.loc[name, 'comp_depth'] = comp_depth
                    
                    try:
                        plt.plot(bl_out[param], bl_out['depSM: Depth [salt water, m]'], 'ro', markersize = 8)
                    except KeyError:
                        plt.plot(0, bl_out['depSM: Depth [salt water, m]'], 'ro', markersize = 8)                   
                    
                    pdf.savefig()
                    plt.close()
                except TypeError:
                    continue
                    
            for param in bvf.columns:
                
                fig, ax = plt.subplots(figsize = (8.5, 11))
                ax.plot(bvf[param], bvf.index)
                ax.invert_yaxis()
                ax.set_xlabel(param)
                ax.set_ylabel('depSM: Depth [salt water, m]')
                plt.title(name)
                
                plt.plot([0] * 24,
                         bl_out['depSM: Depth [salt water, m]'], 'ro', markersize = 8)
                pdf.savefig()
                plt.close()
                
bl_out.to_csv('Z://public//CTD//mgl2506_bottles.csv')
cast_out.to_csv('Z://public//CTD//mgl2506_castmetadata.csv')

## master T-S plot

with PdfPages('Z://public//CTD//all_TS_plots.pdf') as pdf:

    viridis = colormaps['viridis'].resampled(100)
    fig, ax = plt.subplots(figsize = (8.5, 8.5))
    ts_scatter = plt.scatter(all_S, all_T, c = all_rho, cmap = 'viridis')
    ax.set_xlabel('S')
    ax.set_ylabel('T')
    cbar = plt.colorbar(ts_scatter)
    cbar.set_label(rotation = 270, label = 'rho')
    #plt.show()
    pdf.savefig()
    plt.close()
                    
            