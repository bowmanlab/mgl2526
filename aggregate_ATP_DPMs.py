# -*- coding: utf-8 -*-
"""
Created on Sat Jul  5 05:25:03 2025

@author: jeff
"""

import pandas as pd
import os

all_data = []

for f in os.listdir():
    if f.endswith('ATP.xls'):
        data = pd.read_excel('2025062804_18_ATP.xls', engine='xlrd', skiprows = 33)
        data['file'] = f
        data.drop(columns = ['SID', 'TDCR', 'LUMI', 'LLD', 'MDA', 'TRH', 'Unnamed: 9'], inplace = True)
        column_order = ['file', 'POS', 'CPM', 'DPM']
        data = data[column_order]
        all_data.append(data)
        
data_out = pd.concat(all_data)

data_out.to_csv('combined_ATP.csv')
        
