# -*- coding: utf-8 -*-
'''
Created on Mon Jul 27 12:52:11 2015

Cell Count Analysis:
Reads through a series of folders (1-n) containing cell count data and produces
tables that can be copied into GraphPad for visualization. Data is structured
as lists of ImageJ ROIs, with each file representing a section of the lateral
amygdala (LA). First entries are the area of the entire counted LA region and
the last 5 entries are samples of background staining intensity.
Each folder stores 3-8 sections for each LA from both the L and R hemisphere.

Naming convention: 
function = snake_case
variable = CamelCase
'''
import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

pd.set_option('display.precision', 5)

def threshold(filename):
    '''
    Calculates the threshold of a set of cells by averaging the last 5 (control) entries
    '''
    threshold = filename['Mean'][(len(filename.index) - 5):len(filename.index)].mean(axis=0)
    
    return threshold

def neun_count(filename, timesbaseline):
    SliceList = pd.DataFrame()
    for n in np.arange(1, 8):
            if not os.path.isfile(filename % n):
                break
            Slice = pd.DataFrame(np.random.randn(1,3),
                columns=['Number', 'Area', 'Density'])
            NeuNCount = pd.read_csv(filename % n, index_col=" ")
            NeuNCells = NeuNCount['Mean'][(NeuNCount['Mean'] > (threshold(NeuNCount)*timesbaseline))]
            Slice['Area'][0] = NeuNCount['Area'][1]
            Slice['Number'][0] = len(NeuNCells.index) + 1
            # 44696.932 transforms to area/100um^2 based on .473um pixel size
            # 178034.139970446 transforms to area/100um based on .237um pixel size
            Slice['Density'][0] = Slice['Number'][0] / (Slice['Area'][0] / 44696.932)
            SliceList = SliceList.append(Slice)
            
    return SliceList
    
    
def amygdala_loop(filename, timesbaseline, AllCells):
    '''
    Loops through amygdala sections,  constructing a DataFrame for each section
    containing Total Cell Density, Area, Number of Cells, Mean Pixel Intensity,
    Arc+ proportion values. At the end of each iteration,
    these DataFrames are appended, producing one that represents all sections
    taken form one hemisphere of one animal.
    '''
    
    SliceList = pd.DataFrame()
    for n in np.arange(1, 8):
            if not os.path.isfile(filename % n):
                break
            Slice = pd.DataFrame(np.random.randn(1,5),
                columns=['Number', 'Area', 'Intensity', 'Density', 'Proportion'])
            ArcCount = pd.read_csv(filename % n, index_col=" ")
            ArcCells = ArcCount['Mean'][(ArcCount['Mean'] > (threshold(ArcCount)*timesbaseline))]
            Slice['Area'][0] = ArcCount['Area'][1]
            Slice['Number'][0] = len(ArcCells.index) + 1
            Slice['Intensity'][0] = ArcCells.mean(axis=0) - threshold(ArcCount)
            # 44696.932 transforms to area/100um^2 based on .473um pixel size
            # 178034.139970446 transforms to area/100um based on .237um pixel size
            Slice['Density'][0] = Slice['Number'][0] / (Slice['Area'][0] / 44696.932)
            Slice['Proportion'][0] = Slice['Density'][0] / 16
            SliceList = SliceList.append(Slice)
            AllCells = np.append(AllCells,ArcCells.tolist()/threshold(ArcCount))
            
    return SliceList, AllCells


def make_tables(Key, CombinedMice, MiceLeft, MiceRight, CombinedAnt, CombinedPost):
    '''
    Constructs data tables for average Arc Cell Density, Area, Arc+ proportion,
    Pixel Intensity, Left & Right cell Density, and Anterior
    & Posterior Cell Density for export to GraphPad
    
    Still clunky. Can't figure out how to make nice tables with different keys
    '''
    CombinedMice.index += 1
    CombinedMice['Name'] = Key['Name']
    CombinedMice['Group'] = Key['Group']
    
    AverageDensity = CombinedMice.pivot('Name','Group', values='Density').sort_index(axis=1,
        ascending=False)
    
    AverageArea = CombinedMice.pivot('Name','Group', values='Area').sort_index(axis=1,
        ascending=False).sort_index(axis=0, ascending=False)
        
    AverageProportion = CombinedMice.pivot('Name','Group', values='Proportion').sort_index(axis=1,
        ascending=False).sort_index(axis=0, ascending=False)
    
    AverageIntensity = CombinedMice.pivot('Name','Group', values='Intensity').sort_index(axis=1,
        ascending=False).sort_index(axis=0, ascending=False)
    
    MiceLeft.index += 1
    MiceLeft['Name'] = Key['Name']
    MiceLeft['Group'] = Key['Group']
    
    LeftDensity = MiceLeft.pivot('Name','Group', values='Density').sort_index(axis=1,
        ascending=False).sort_index(axis=0, ascending=False)
            
    MiceRight.index += 1
    MiceRight['Name'] = Key['Name']
    MiceRight['Group'] = Key['Group']
    
    RightDensity = MiceRight.pivot('Name','Group', values='Density').sort_index(axis=1,
        ascending=False).sort_index(axis=0, ascending=False)
        
    CombinedAnt.index += 1
    CombinedAnt['Name'] = Key['Name']
    CombinedAnt['Group'] = Key['Group']
    
    AntDensity = CombinedAnt.pivot('Name','Group', values='Density').sort_index(axis=1,
        ascending=False).sort_index(axis=0, ascending=False)
        
    CombinedPost.index += 1
    CombinedPost['Name'] = Key['Name']
    CombinedPost['Group'] = Key['Group']
    
    PostDensity = CombinedPost.pivot('Name','Group', values='Density').sort_index(axis=1,
        ascending=False).sort_index(axis=0, ascending=False)
            
    return  AverageDensity, AverageArea, AverageProportion, AverageIntensity, LeftDensity, RightDensity, AntDensity, PostDensity, CombinedMice
      
def make_histogram(x):
    # the histogram of the data
    n, bins, patches = plt.hist(x, 100, normed=1, facecolor='green', alpha=0.75)

    plt.xlabel('Intensity / threshold')
    plt.ylabel('Frequency')
    plt.axis([1.1, 3,0,2.5])
    plt.show()
    
    
def cell_count(folder, region, timesbaseline):
    '''
    Loops through all folders in 'To Analyze' folder, combining average Right &
    Left hemisphere data from all mice into one DataFrame
    '''
    
    # Instantiate all of the global Dataframes
    
    AllCells = np.array([])
    MiceLeft = pd.DataFrame()
    MiceRight = pd.DataFrame()
    AntLeft = pd.DataFrame()
    PostLeft = pd.DataFrame()
    AntRight = pd.DataFrame()
    PostRight = pd.DataFrame()
    
    if folder == 'LA':
        os.chdir('LA Count')
    if folder == 'BA':
        os.chdir('BA Count')
    if folder == 'PV':
        os.chdir('PV Count')
        
    if region == 'LA':
        region = ''
    
    
    # Left
    i = 1
    # Loop through folders, building MiceLeft, AntLeft, and Post,Left
    while (i <= len(os.listdir())):
        os.chdir('%d' % i)
        Amygdala = pd.DataFrame()

        # Loop through slices, building up Amygdala
        Amygdala, AllCells = amygdala_loop('%dl' + region + 'Count.csv', timesbaseline, AllCells)
        MiceLeft = MiceLeft.append(Amygdala.mean(axis=0), ignore_index=True)
        
        # Seperate anterior and posterior slices
        AntLeft = AntLeft.append(
            Amygdala[0:(int(len(Amygdala.index) / 2))].mean(axis=0), ignore_index=True)
        PostLeft = PostLeft.append(
            Amygdala[(int(len(Amygdala.index) / 2)):len(Amygdala.index)].mean(axis=0), ignore_index=True)
        os.chdir(os.pardir)
        i = i + 1

    # Right
    i = 1
    # Loop through folders, building up MiceRight, AntRight, and PostRight
    while (i <= len(os.listdir())):
        os.chdir('%d' % i)
        Amygdala = pd.DataFrame()

        # Loop through slices
        Amygdala, AllCells = amygdala_loop('%dr' + region +'Count.csv', timesbaseline, AllCells)
        MiceRight = MiceRight.append(Amygdala.mean(axis=0), ignore_index=True)
        
        # Seperate anterior and posterior slices
        AntRight = AntRight.append(
            Amygdala[0:(int(len(Amygdala.index) / 2))].mean(axis=0), ignore_index=True)
        PostRight = PostRight.append(Amygdala[(int(len(Amygdala.index) / 2)):len(Amygdala.index)].mean(axis=0), ignore_index=True)
        os.chdir(os.pardir)
        i = i + 1


    
    CombinedMice = pd.concat([MiceLeft, MiceRight]).groupby(level=0).mean()
    CombinedAnt = pd.concat([AntLeft, AntRight]).groupby(level=0).mean()
    CombinedPost= pd.concat([PostLeft, PostRight]).groupby(level=0).mean()
    os.chdir(os.pardir)
    
    if folder == 'PV':
       Key = pd.read_csv('PVKey.csv',index_col='Number')
       
    if folder == 'BA':
       Key = pd.read_csv('BAKey.csv',index_col='Number')
        
    if folder == 'LA':
       Key = pd.read_csv('LAKey.csv',index_col='Number')
       
    # Make data tables that can be copied into GraphPad   
    return make_tables(Key, CombinedMice, MiceLeft, MiceRight, CombinedAnt, CombinedPost) 

if __name__ == "__main__":
    # LA Threshold = 1.5
    #BA Threshold = 1.2
    #PV Threshold = 1.15
    AverageDensity, AverageArea, AverageProportion, AverageIntensity, LeftDensity, RightDensity, AntDensity, PostDensity, CombinedMice = cell_count('PV','LA', 1.15)
    
    