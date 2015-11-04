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
    Arc+ cell Density and Arc+ proportion values. At the end of each iteration,
    these DataFrames are appended, producing one that represents all sections
    taken form one hemisphere of one animal.
    '''
    
    SliceList = pd.DataFrame()
    for n in np.arange(1, 8):
            if not os.path.isfile(filename % n):
                break
            Slice = pd.DataFrame(np.random.randn(1,6),
                columns=['Number', 'Area', 'Intensity', 'Density', 'Proportion', 'Cell Density'])
            ArcCount = pd.read_csv(filename % n, index_col=" ")
            ArcCells = ArcCount['Mean'][(ArcCount['Mean'] > (threshold(ArcCount)*timesbaseline))]
            Slice['Cell Density'][0] = 16
            Slice['Area'][0] = ArcCount['Area'][1]
            Slice['Number'][0] = len(ArcCells.index) + 1
            Slice['Intensity'][0] = ArcCells.mean(axis=0) - threshold(ArcCount)
            # 44696.932 transforms to area/100um^2 based on .473um pixel size
            # 178034.139970446 transforms to area/100um based on .237um pixel size
            Slice['Density'][0] = Slice['Number'][0] / (Slice['Area'][0] / 44696.932)
            Slice['Proportion'][0] = Slice['Density'][0] / Slice['Cell Density'][0]
            SliceList = SliceList.append(Slice)
            AllCells = AllCells + ArcCells.tolist()
            
    return SliceList, AllCells

def make_tables_BA(CombinedMice, MiceLeft, MiceRight, CombinedAnt, CombinedPost):
    
    hc = [20,21,22,23,24,25,26,27]
    low = [13,14,15,16,17,18,19]
    med = [5,6,7,8,9,10,11,12]
    high = [0,1,2,3,4]
    
    AverageDensity = pd.DataFrame(        [
    
            CombinedMice['Density'][hc].tolist(),
            CombinedMice['Density'][low].tolist(),
            CombinedMice['Density'][med].tolist(),
            CombinedMice['Density'][high].tolist()],
            index=['hc','low','med','high']).transpose()
            
    AverageArea = pd.DataFrame(        [
            CombinedMice['Area'][hc].tolist(),
            CombinedMice['Area'][low].tolist(),
            CombinedMice['Area'][med].tolist(),
            CombinedMice['Area'][high].tolist()],
            index=['hc','low','med','high']).transpose()
            
    AverageProportion = pd.DataFrame(        [
            CombinedMice['Proportion'][hc].tolist(),
            CombinedMice['Proportion'][low].tolist(),
            CombinedMice['Proportion'][med].tolist(),
            CombinedMice['Proportion'][high].tolist()],
            index=['hc','low','med','high']).transpose()
            
    AverageIntensity = pd.DataFrame(        [
            CombinedMice['Intensity'][hc].tolist(),
            CombinedMice['Intensity'][low].tolist(),
            CombinedMice['Intensity'][med].tolist(),
            CombinedMice['Intensity'][high].tolist()],
            index=['hc','low','med','high']).transpose()
                    
    
    LeftDensity = pd.DataFrame(        [
            MiceLeft['Density'][hc].tolist(),
            MiceLeft['Density'][low].tolist(),
            MiceLeft['Density'][med].tolist(),
            MiceLeft['Density'][high].tolist()],
            index=['hc','low','med','high']).transpose()
            
            
    RightDensity = pd.DataFrame(        [
            MiceRight['Density'][hc].tolist(),
            MiceRight['Density'][low].tolist(),
            MiceRight['Density'][med].tolist(),
            MiceRight['Density'][high].tolist()],
            index=['hc','low','med','high']).transpose()
            
    AntDensity = pd.DataFrame(        [
            CombinedAnt['Density'][hc].tolist(),
            CombinedAnt['Density'][low].tolist(),
            CombinedAnt['Density'][med].tolist(),
            CombinedAnt['Density'][high].tolist()],
            index=['hc','low','med','high']).transpose()
            
    PostDensity = pd.DataFrame(        [
            CombinedPost['Density'][hc].tolist(),
            CombinedPost['Density'][low].tolist(),
            CombinedPost['Density'][med].tolist(),
            CombinedPost['Density'][high].tolist()],
            index=['hc','low','med','high']).transpose()
            
    return  AverageDensity, AverageProportion, AverageIntensity, LeftDensity, RightDensity, AverageArea, AntDensity, PostDensity, CombinedMice   
    
    
def make_tables_PV(CombinedMice, MiceLeft, MiceRight, CombinedAnt, CombinedPost):
    saline = [2, 5, 6, 7]
    cno = [0, 1, 3, 4, 8] 
    
    AverageDensity = pd.DataFrame(        [
            CombinedMice['Density'][saline].tolist(),
            CombinedMice['Density'][cno].tolist()],
            index=['saline','cno']).transpose()
            
    AverageArea = pd.DataFrame(        [
            CombinedMice['Area'][saline].tolist(),
            CombinedMice['Area'][cno].tolist()],
            index=['saline','cno']).transpose()
    
    AverageProportion = pd.DataFrame(        [
            CombinedMice['Proportion'][saline].tolist(),
            CombinedMice['Proportion'][cno].tolist()],
            index=['saline','cno']).transpose()

    AverageIntensity = pd.DataFrame(        [
            CombinedMice['Intensity'][saline].tolist(),
            CombinedMice['Intensity'][cno].tolist()],
            index=['saline','cno']).transpose()
            
    LeftDensity = pd.DataFrame(        [
            MiceLeft['Density'][saline].tolist(),
            MiceLeft['Density'][cno].tolist()],
            index=['saline','cno']).transpose()
            
            
    RightDensity = pd.DataFrame(        [
            MiceRight['Density'][saline].tolist(),
            MiceRight['Density'][cno].tolist()],
            index=['saline','cno']).transpose()
            
    AntDensity = pd.DataFrame(        [
            CombinedAnt['Density'][saline].tolist(),
            CombinedAnt['Density'][cno].tolist()],
            index=['saline','cno']).transpose()
            
    PostDensity = pd.DataFrame(        [
            
            CombinedPost['Density'][saline].tolist(),
            CombinedPost['Density'][cno].tolist()],
            index=['saline','cno']).transpose()
            
                
    return  AverageDensity, AverageProportion, AverageIntensity, LeftDensity, RightDensity, AverageArea, AntDensity, PostDensity, CombinedMice
    
def make_tables_LA(CombinedMice, MiceLeft, MiceRight, CombinedAnt, CombinedPost):
    '''
    Constructs data tables for average Arc Cell Density, LA Area, Arc+ proportion,
    Pixel Intensity, Total Cell Density, Left & Right cell Density, and Anterior
    & Posterior Cell Density for export to GraphPad
    
    Clunky implementation currently
    '''
    hc = 4
    low = 20
    med = 18
    high = 10
    
    AverageDensity = pd.DataFrame(        [
            CombinedMice['Density'][0:hc].tolist(),
            CombinedMice['Density'][low:].tolist(),
            CombinedMice['Density'][high:med].tolist(),
            CombinedMice['Density'][hc:high].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    AverageArea = pd.DataFrame(        [
            CombinedMice['Area'][0:hc].tolist(),
            CombinedMice['Area'][low:].tolist(),
            CombinedMice['Area'][high:med].tolist(),
            CombinedMice['Area'][hc:high].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    AverageProportion = pd.DataFrame(        [
            CombinedMice['Proportion'][0:hc].tolist(),
            CombinedMice['Proportion'][low:].tolist(),
            CombinedMice['Proportion'][high:med].tolist(),
            CombinedMice['Proportion'][hc:high].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    AverageIntensity = pd.DataFrame(        [
            CombinedMice['Intensity'][0:hc].tolist(),
            CombinedMice['Intensity'][low:].tolist(),
            CombinedMice['Intensity'][high:med].tolist(),
            CombinedMice['Intensity'][hc:high].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
    
    LeftDensity = pd.DataFrame(        [
            MiceLeft['Density'][0:hc].tolist(),
            MiceLeft['Density'][low:].tolist(),
            MiceLeft['Density'][high:med].tolist(),
            MiceLeft['Density'][hc:high].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    CellDensity = pd.DataFrame(        [
            CombinedMice['Cell Density'][0:hc].tolist(),
            CombinedMice['Cell Density'][low:].tolist(),
            CombinedMice['Cell Density'][high:med].tolist(),
            CombinedMice['Cell Density'][hc:high].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    RightDensity = pd.DataFrame(        [
            MiceRight['Density'][0:hc].tolist(),
            MiceRight['Density'][low:].tolist(),
            MiceRight['Density'][high:med].tolist(),
            MiceRight['Density'][hc:high].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    AntDensity = pd.DataFrame(        [
            CombinedAnt['Density'][0:hc].tolist(),
            CombinedAnt['Density'][low:].tolist(),
            CombinedAnt['Density'][high:med].tolist(),
            CombinedAnt['Density'][hc:high].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    PostDensity = pd.DataFrame(        [
            CombinedPost['Density'][0:hc].tolist(),
            CombinedPost['Density'][low:].tolist(),
            CombinedPost['Density'][high:med].tolist(),
            CombinedPost['Density'][hc:high].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    return  AverageDensity, AverageProportion, AverageIntensity, LeftDensity, RightDensity, AverageArea, AntDensity, PostDensity, CombinedMice
      
    
def cell_count(folder, region, timesbaseline):
    '''
    Loops through all folders in 'To Analyze' folder, combining average Right &
    Left hemisphere data from all mice into one DataFrame
    '''
    
    # Instantiate all of the global Dataframes
    
    AllCells = []
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
    

    # Make data tables that can be copied into GraphPad
    if folder == 'PV':
       return make_tables_PV(CombinedMice, MiceLeft, MiceRight, CombinedAnt, CombinedPost)
    if folder == 'BA':
       return make_tables_BA(CombinedMice, MiceLeft, MiceRight, CombinedAnt, CombinedPost) 
    if folder == 'LA':
       return make_tables_LA(CombinedMice, MiceLeft, MiceRight, CombinedAnt, CombinedPost) 

if __name__ == "__main__":

    Density, Proportion, Intensity, LeftDensity, RightDensity, Area, AntDensity, PostDensity, CombinedMice = cell_count('LA','LA', 1.5)
    