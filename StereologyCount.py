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
    Creates a list of dataframes containing all counting data for particular region
    '''
    os.chdir(filename)
    framelist = []
    
    for i in os.listdir():
        if region in i:
            frame = pd.read_csv(i, index_col='Marker', engine='python')
            framelist.append(frame)
    os.chdir(os.pardir)
            
    return framelist
    
def create_ratioframe(framelist):
    '''
    loops through framelist, creating a dataframe with ratio data for each mouse
    '''
    ratioframe = pd.DataFrame(np.zeros((0,5)),columns = ['Arc/Neun','Arc/DAPI','NeuN/DAPI','Estimated Population using Mean Section Thickness','Coefficient of Error (Gundersen), m=1'])
    ratios = pd.DataFrame(np.zeros((1,5)),columns = ['Arc/Neun','Arc/DAPI','NeuN/DAPI','Estimated Population using Mean Section Thickness','Coefficient of Error (Gundersen), m=1'])
    indexframe = pd.DataFrame(np.zeros((0,1)),columns = ['Number'])
    index = pd.DataFrame(np.zeros((1,1)),columns = ['Number'])
    for i in framelist:
        if 'colabel' in i.index.values.tolist():
            arcneun = (i['Total Markers Counted']['brdu'] / (i['Total Markers Counted']['colabel']))*100
            arcdapi = np.nan
            neundapi = np.nan
            estpop = i['Estimated Population using Mean Section Thickness']['brdu']
            ce = i['Coefficient of Error (Gundersen), m=1']['brdu']
            if 'cfos' in i.index.values.tolist():
                arcdapi = (i['Total Markers Counted']['brdu'] / i['Total Markers Counted']['cfos'])*100
                neundapi = (i['Total Markers Counted']['colabel'] / i['Total Markers Counted']['cfos'])*100
                
        ratios.iloc[0] = [arcneun,arcdapi,neundapi,estpop, ce]
        ratioframe = ratioframe.append(ratios)
        data_file = i['Data File']['brdu']
        index.iloc[0] = [int(data_file[0:(len(data_file)-4)])]
        indexframe = indexframe.append(index)
    
    return ratioframe, indexframe
  
def make_tables(ratioframe, indexframe):
    '''
    '''
    Key = pd.read_csv('Stereology Key.csv',index_col='Number')
    ratioframe.index=indexframe['Number']
    ratioframe['Name'] = Key['Name']
    ratioframe['Group'] = Key['Group']
    
    ArcNeunRatio = ratioframe.pivot('Name','Group', values='Arc/Neun').sort_index(axis=1,
        ascending=False)
        
    NeunDAPIRatio = ratioframe.pivot('Name','Group', values='NeuN/DAPI').sort_index(axis=1,
        ascending=False)
        
    EstPop = ratioframe.pivot('Name','Group', values='Estimated Population using Mean Section Thickness').sort_index(axis=1,
        ascending=False)
    
    
    return ArcNeunRatio, NeunDAPIRatio, EstPop

if __name__ == "__main__":
    framelist = get_frame('Memory Strength Counting Data','LA')
    ratioframe, indexframe = create_ratioframe(framelist)
    ArcNeunRatio, NeunDAPIRatio, EstPop = make_tables(ratioframe, indexframe)