# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 12:52:11 2015

Stereology Analysis:
Reads through a series of csvs containing exported StereoInvestigator cell
count data and produces a DataFrame containing information on ratios between
cell tyoes

Naming convention: 
function = snake_case
variable = CamelCase
"""

import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

pd.set_option('display.precision', 5)

def get_frame(filename, region):
    '''
    Loops through all folders in 'To Analyze' folder, combining average Right &
    Left hemisphere data from all mice into one DataFrame
    '''
    os.chdir(filename)
    framelist = []
    
    for i in os.listdir():
        if region in i:
            frame = pd.read_csv(i, index_col='Marker', engine='python')
            framelist.append(frame)
            
    return framelist
    
def create_ratioframe(framelist):
    '''
    loops through framelist, creating a dataframe with ratio data for each mouse
    '''
    ratioframe = pd.DataFrame(np.zeros((0,3)),columns = ['Arc/Neun','Arc/DAPI','NeuN/DAPI'])
    ratios = pd.DataFrame(np.zeros((1,3)),columns = ['Arc/Neun','Arc/DAPI','NeuN/DAPI'])
    for i in framelist:
        if 'colabel' in i.index.values.tolist():
            arcneun = i['Total Markers Counted']['brdu'] / (i['Total Markers Counted']['colabel']*.8)
            arcdapi = np.nan
            neundapi = np.nan
            if 'cfos' in i.index.values.tolist():
                arcdapi = i['Total Markers Counted']['brdu'] / i['Total Markers Counted']['cfos']
                neundapi = i['Total Markers Counted']['colabel'] / i['Total Markers Counted']['cfos']
        else:
            arcdapi = i['Total Markers Counted']['brdu'] / i['Total Markers Counted']['cfos']
            arcneun = 'na'
            neundapi = 'na'
        
        ratios.iloc[0] = [arcneun,arcdapi,neundapi]
        ratioframe = ratioframe.append(ratios)
    
    return ratioframe
  


if __name__ == "__main__":
    framelist = get_frame('Memory Strength Counting Data','LA')
    ratioframe = create_ratioframe(framelist)

