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
from matplotlib import lines
import glob
import seawater
import itertools
from datetime import datetime
from datetime import timezone
from zoneinfo import ZoneInfo

bl_names = ['bottle', 'bottle_2', 'date.time', 'index_start', 'index_stop']

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
    
def parse_header_file(directory, name):
    with open(directory + name + '.hdr', 'r') as header_file:
        for line in header_file:
            if 'NMEA Latitude' in line:
                line = line.split()
                deg_lat = line[4]
                min_lat = line[5]
                ddeg_lat = round(float(deg_lat) + float(min_lat)/60, 5)
            elif 'NMEA Longitude' in line:
                line = line.split()
                deg_lon = line[4]
                min_lon = line[5]
                ddeg_lon = round(float(deg_lon) + float(min_lon)/60, 5)     
            elif 'NMEA UTC (Time)' in line:
                line = line.split(' = ')[1]
                line = line + ' +0000'
                time_utc = datetime.strptime(line, '%b %d %Y %H:%M:%S %z')
                time_local = time_utc.astimezone(ZoneInfo("America/Los_Angeles"))
                time_utc_out = datetime.strftime(time_utc, '%Y-%m-%d %H:%M:%S %z')
                time_local_out = datetime.strftime(time_local, '%Y-%m-%d %H:%M:%S %z')
                
    header_out = {'lat':ddeg_lat,
                  'lon':ddeg_lon,
                  'time_utc': time_utc_out,
                  'time_local_out': time_local_out}
            
    return(header_out)

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
        
        ## parse the header file 
        
        header = parse_header_file(bl_input_dir, name)
        
        for key in header.keys():
            cast_out.loc[name, key] = header[key]
        
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
        
        ## % irradiance
        
        cnv_downcast['irradiance%'] = (cnv_downcast['par: PAR/Irradiance, Biospherical/Licor']/cnv_in['spar: SPAR, Biospherical/Licor']) * 100
        
        ## bin data by 2bdar
        
        cnv_downcast['bin'] = pd.cut(cnv_downcast['depSM: Depth [salt water, m]'], bins=300)
        cnv_binned = cnv_downcast.groupby('bin', observed = True).mean()
        
        s_binned = cnv_binned['sal00: Salinity, Practical [PSU]']
        t_binned = cnv_binned['t090C: Temperature [ITS-90, deg C]']
        p_binned = cnv_binned['depSM: Depth [salt water, m]']
            
        ## calculate buoyancy frequency and find MLD
                
        all_T = all_T + list(t)
        all_S = all_S + list(s)
        all_rho = all_rho + list(cnv_downcast['rho [kg/m^3]'])
        
        n2, q, p_ave = seawater.geostrophic.bfrq(s_binned, t_binned, p_binned)
        n2_flat = list(itertools.chain.from_iterable(n2))
        q_flat = list(itertools.chain.from_iterable(q))
        p_flat = list(itertools.chain.from_iterable(p_ave))
        bvf = pd.DataFrame(zip(p_flat, n2_flat, q_flat), columns = ['p', 'n2', 'q'])
        ml_depth = bvf['p'][bvf.n2.idxmax()]
        
        ## find cmax and max c
        
        cmax = cnv_downcast['depSM: Depth [salt water, m]'][cnv_downcast['flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]'].idxmax()]
        maxc = cnv_downcast['flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]'].max()
        
        ## find compensation depth (1 % light level)
        
        temp_delta = abs(cnv_downcast['irradiance%'] - 1)
        temp_i = temp_delta.idxmin()
        comp_depth = cnv_downcast['depSM: Depth [salt water, m]'][temp_i]

        ## add cmax, mld, maxc, and compensation depth to cast_out
        
        cast_out.loc[name, 'ml_depth_[m]'] = ml_depth
        cast_out.loc[name, 'cmax_[m]'] = cmax
        cast_out.loc[name, 'maxc_[mg/m^3]'] = maxc
        cast_out.loc[name, 'comp_depth_[m]'] = comp_depth
        
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
        
        ## make some depth profiles
            
        with PdfPages('Z://public//CTD//' + name + '_plots.pdf') as pdf:
            for param in ['CStarTr0: Beam Transmission, WET Labs C-Star [%]',
                   'c0S/m: Conductivity [S/m]', 'c1S/m: Conductivity, 2 [S/m]',
                   'flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]',
                   'sbox0Mm/Kg: Oxygen, SBE 43 [umol/kg]',
                   'sbox1Mm/Kg: Oxygen, SBE 43, 2 [umol/kg]',
                   'sal00: Salinity, Practical [PSU]',
                   'sal11: Salinity, Practical, 2 [PSU]',
                   't090C: Temperature [ITS-90, deg C]',
                   't190C: Temperature, 2 [ITS-90, deg C]',
                   'turbWETntu0: Turbidity, WET Labs ECO [NTU]',
                   'par: PAR/Irradiance, Biospherical/Licor',
                   'rho [kg/m^3]', 'irradiance%']:
                
                fig, ax = plt.subplots(figsize = (8.5, 11))
                ax.plot(cnv_downcast[param], cnv_downcast['depSM: Depth [salt water, m]'])
                ax.invert_yaxis()
                ax.set_xlabel(param)
                ax.set_ylabel('depSM: Depth [salt water, m]')
                plt.title(name)
                
                line_colors = {'ml_depth_[m]':'b' ,
                               'cmax_[m]':'g',
                               'comp_depth_[m]':'k'}
                
                try:
                    plt.plot(bl_out[param], bl_out['depSM: Depth [salt water, m]'], 'ro', markersize = 8)
                except KeyError:
                    plt.plot(0, bl_out['depSM: Depth [salt water, m]'], 'ro', markersize = 8)
                    
                for special_param in line_colors.keys():
                    special_param_z = cast_out.loc[name, special_param]
                    plt.plot([0, max(cnv_downcast[param])], [special_param_z, special_param_z],
                             color = line_colors[special_param],
                             linestyle = 'dashed')
                    
                legend_lines = [lines.Line2D([0], [0], color = 'b', linestyle = 'dashed'),
                                lines.Line2D([0], [0], color = 'g', linestyle = 'dashed'),
                                lines.Line2D([0], [0], color = 'k', linestyle = 'dashed')]
                    
                ax.legend(legend_lines,
                          ['ml_depth_[m]', 'cmax_[m]', 'comp_depth_[m]'],
                          loc = 'lower right')

                pdf.savefig()
                plt.close()
                    
            for param in ['n2']:
                
                fig, ax = plt.subplots(figsize = (8.5, 11))
                ax.plot(bvf[param], bvf['p'])
                ax.invert_yaxis()
                ax.set_xlabel(param)
                ax.set_ylabel('depSM: Depth [salt water, m]')
                plt.title(name)
                
                plt.plot([0] * 24,
                         bl_out['depSM: Depth [salt water, m]'], 'ro', markersize = 8)
                
                for special_param in line_colors.keys():
                    special_param_z = cast_out.loc[name, special_param]
                    plt.plot([0, max(bvf[param])], [special_param_z, special_param_z],
                             color = line_colors[special_param],
                             linestyle = 'dashed')
                    
                    ax.legend(legend_lines,
                              ['ml_depth_[m]', 'cmax_[m]', 'comp_depth_[m]'],
                              loc = 'lower right')
                    
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
                    
            