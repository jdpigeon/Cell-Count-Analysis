# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 16:17:54 2016

Subnuclei Analysis:
Calculates the Arc+ proportion of cells in amygdala slice counts, broken down by
subnuclei. Implements two different forms of subnuclei assignment depending on
whether subnuclei were delineated with contours during counting or post-hoc.
"""

import csv
import pandas as pd
import numpy as np
import os

def get_frames(filename):
    '''
    Extracts arc and dapi counts from csvs containing stereology count data.
    Excises sections containing 'cfos' (DAPI) and 'brdu' (arc) counts and returns
    arc and dapi DataFrames.
    '''
    search_terms = ['cfos Counts By Site ','cfos CE Gundersen ','brdu Counts By Site ', 'brdu CE Gundersen ']
    stop_terms = ['Total']
    slice_list = []

    with open(filename) as inf:
        reader=csv.reader(inf)
        i = 1
        for row in reader:
            if row[0] in search_terms:
                slice_list.append(i)
                
            i = i+1
                        
    arc = pd.read_csv(filename,skiprows=slice_list[0],nrows=slice_list[1]-(slice_list[0]+6)).dropna(1,'all')
    
    dapi = pd.read_csv(filename,skiprows=slice_list[2],nrows=slice_list[3]-(slice_list[2]+6)).dropna(1,'all')
    
    ''' Slice and restructure dataframes so they only contain numeric values'''
    #arc = arc[0:(arc['Run 1'].count()-1)]
    del arc['Marker']
    arc = arc.set_index('Site')
    del dapi['Marker']
    dapi = dapi.set_index('Site')
    dapi = pd.to_numeric(dapi, errors='coerce')
    return dapi, arc
    
def get_key(filename):
    '''
    Checks for post-hoc or contour subnuclei assignment key and
    passes off to appropriate function
    '''
    key_filename = filename[0:len(filename)-6]+'Key.csv'
    if os.path.isfile(key_filename):
        return posthoc_key(key_filename)
    else:
        return contour_key(filename)
    
def posthoc_key(filename):
    '''
    For post-hoc subnuclei assignment, read key csv
    '''
    key = pd.read_csv(filename)
     
    key = key.set_index('Site')
    return key

def contour_key(filename):
    '''
    For contour subnuclei assignment, read key csv
    '''
    search_terms = ['Section Order ']
    s = 0
    
    with open(filename) as inf:
        reader=csv.reader(inf)
        i = 1
        for row in reader:
            if row[0] in search_terms:
                s = i
                
            i = i+1
    ''' Manually assigning column names here because blank space in csv prevents\
    utf-8 codec from decoding it correctly'''
    key = pd.read_csv(filename,skiprows=s+1,header=None).dropna(1,'all')
    key.columns=['Blank','Section','Contour Name','Z of Section (µm)','Average Z of Contour (µm)']	

    key = key['Contour Name']
    
    return key
    
def compute_subnuclei_sums(df,key):
    '''
    Checks whether data is contour or post-hoc and passes to appropriate function
    '''
    # (len(key) == len(arc.columns)), "data and key do not match"
    if isinstance(key, pd.DataFrame):
        return posthoc_subnuclei_sums(df,key)
    else:
        return contour_subnuclei_sums(df,key)
 

def posthoc_subnuclei_sums(df, key):
    '''
    Loops through both columns and rows of key dataframe, adding corresponding
    value in arc or dapi dataframe based on value of key
    '''
    
    key.columns = df.columns.values.tolist()
    
    ladsum = 0
    lamsum = 0
    lavsum = 0
    
    
    for col,row in key.iteritems():
        #row = pd.to_numeric(row,errors='coerce')
        df[col] = pd.to_numeric(df[col],errors='coerce')
        i = 1
        for y in row:
            print(col)
            print (df[col][i])
            print(y)
            if df[col][i] != df[col][i]:
                break
            
            if 'd' in y:
                ladsum = ladsum + int(df[col][i])
            elif 'm' in y:
                lamsum = lamsum + int(df[col][i])
            elif 'v' in y:
                lavsum = lavsum + int(df[col][i])
            
            i = i+1
        
            
    return ladsum, lamsum, lavsum
'''    
def ant_post(df, key):
    
    loops through arc or dapi dataframe, creating sums based on position before
    or after ant/post split point indicated in key
    
    antsum = 0
    postsum = 0
    
    for col, row in df.iteritems():
        i = 1
        row = pd.to_numeric(row,errors='coerce')
        for y in row:
            if len(df.columns)-i >= key:
                postsum = postsum + y
            else:
                antsum = antsum + y
    
    return antsum, postsum
'''
       
def contour_subnuclei_sums(df, key):
    '''
    Loops through arc or dapi dataframe, adding cell counts together based on
    value of key at same position in dataframe.
    '''
    
    ladsum = 0
    lamsum = 0
    lavsum = 0
    i = 0
    for col,row in df.iteritems():
        row = pd.to_numeric(row,errors='coerce')
        
        if 'd' in key[i]:
            ladsum = ladsum + np.sum(row)
        
        elif 'm' in key[i]:
            lamsum = lamsum + np.sum(row)
        elif 'v' in key[i]:
            lavsum = lavsum + np.sum(row)
    
        i = i+1
    return ladsum, lamsum, lavsum
    
    
if __name__ == "__main__":
    os.chdir('PV Subnuclei')
    ratioframe = pd.DataFrame(np.zeros((0,3)),columns = ['dorsal','medial','ventral'])
    for i in np.arange(1,24):
        filename = str(i) + 'LA.csv'
        print(filename)
        if os.path.isfile(filename):
            arc, dapi = get_frames(filename)
            key = get_key(filename)
            dapi_d,dapi_m,dapi_v = compute_subnuclei_sums(dapi, key)
            arc_d,arc_m,arc_v = compute_subnuclei_sums(arc, key)
            nucleiratios = pd.DataFrame([(arc_d/dapi_d)*100, (arc_m/dapi_m)*100, (arc_v/dapi_v)*100]).transpose()
            nucleiratios.columns = ['dorsal','medial','ventral']
            nucleiratios.index = [i]
            ratioframe = ratioframe.append(nucleiratios)
            ratioframe['all']=ratioframe.mean(axis=1)
            
    
            
    
    
    