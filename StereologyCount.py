# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 12:52:11 2015

Stereology Analysis:
Reads through a series of csvs containing exported StereoInvestigator cell
count data and produces a set of DataFrames containing information on ratios between available
cell types (Arc/NeuN, Arc/DAPI, NeuN/DAPI)
"""

import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


def get_frame(filename, region):
    '''
    Creates a list of dataframes containing all counting data for desired region
    '''
    os.chdir(filename)
    framelist = []
    indexframe = pd.DataFrame(np.zeros((0,1)),columns = ['Number'])
    index = pd.DataFrame(np.zeros((1,1)),columns = ['Number'])
    
    for i in np.arange(1,24):
        filename = str(i) + region + '.csv'
        print(filename)
        if os.path.isfile(filename):
            frame = pd.read_csv(str(i)+region+'.csv', index_col='Marker', engine='python')
            framelist.append(frame)
            index.iloc[0] = i
            indexframe = indexframe.append(index)
    os.chdir(os.pardir)
            
    return framelist, indexframe
    
def create_ratioframe(framelist):
    '''
    loops through framelist, creating a dataframe with ratio data for each mouse
    '''
    ratioframe = pd.DataFrame(np.zeros((0,5)),columns = ['Arc/Neun','Arc/DAPI','NeuN/DAPI',
    'Estimated Population using Mean Section Thickness','Coefficient of Error (Gundersen), m=1'])
    ratios = pd.DataFrame(np.zeros((1,5)),columns = ['Arc/Neun','Arc/DAPI','NeuN/DAPI',
    'Estimated Population using Mean Section Thickness','Coefficient of Error (Gundersen), m=1'])
    
    
    for i in framelist:
        
        if 'colabel' in i.index.values.tolist():
            arcneun = (i['Total Markers Counted']['brdu'] / (i['Total Markers Counted']['colabel']))*100
            arcdapi = np.nan
            ce = np.nan
            neundapi = np.nan
            estpop = i['Estimated Population using Mean Section Thickness']['brdu']
            ce = i['Coefficient of Error (Gundersen), m=1']['brdu']
            if 'cfos' in i.index.values.tolist():
                arcdapi = (i['Total Markers Counted']['brdu'] / i['Total Markers Counted']['cfos'])*100
                neundapi = (i['Total Markers Counted']['colabel'] / i['Total Markers Counted']['cfos'])*100
        elif 'cfos' in i.index.values.tolist():
            arcdapi = (i['Total Markers Counted']['brdu'] / (i['Total Markers Counted']['cfos']))*100
            estpop = i['Estimated Population using Mean Section Thickness']['brdu']
            arcneun = np.nan
            neundapi = np.nan
            print(i.index.values.tolist())
            ce = 'placehold'
                
        ratios.iloc[0] = [arcneun,arcdapi,neundapi,estpop, ce]
        ratioframe = ratioframe.append(ratios)
        
        #Stores blinded index number of each mouse by removing the final 4chars (ie. '.DAT') from filename
        
        
    
    return ratioframe

def assign_groups(key, ratioframe, index):
    '''
    Assigns real mouse names and groups to ratioframe based on key
    '''
    Key = pd.read_csv(key,index_col='Number')
    ratioframe.index=index['Number']
    ratioframe['Name'] = Key['Name']
    ratioframe['Group'] = Key['Group']

    return ratioframe
    
def make_tables(ratioframe):
    '''
    Constructs DataFrame tables of ratio data based on group ID
    '''
    
    ArcNeunRatio = ratioframe.pivot('Name','Group', values='Arc/Neun').sort_index(axis=1,
        ascending=False)
        
    NeunDAPIRatio = ratioframe.pivot('Name','Group', values='NeuN/DAPI').sort_index(axis=1,
        ascending=False)
        
    EstPop = ratioframe.pivot('Name','Group', values='Estimated Population using Mean Section Thickness').sort_index(axis=1,
        ascending=False)
        
    ArcDapi = ratioframe.pivot('Name','Group', values='Arc/DAPI').sort_index(axis=1,
        ascending=False)
    
    return ArcNeunRatio, NeunDAPIRatio, EstPop, ArcDapi

if __name__ == "__main__":
    framelist, index = get_frame('PV Total Proportion','BA')
    ratioframe  = create_ratioframe(framelist)
    ratioframe = assign_groups('PV Stereology Key.csv', ratioframe, index)
    ArcNeunRatio, NeunDAPIRatio, EstPop, ArcDapi = make_tables(ratioframe)
    
    ax = sns.barplot(x='Group', y='Arc/DAPI', data=ratioframe)