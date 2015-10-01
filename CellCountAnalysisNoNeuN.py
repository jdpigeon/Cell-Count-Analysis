# -*- coding: utf-8 -*-
"""
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
"""
import pandas as pd
import numpy as np
import seaborn as sns
import os

pd.set_option('display.precision', 5)

# Calculates the threshold of a set of cells by averaging
# the last 5 (control) entries


def threshold(filename):

    threshold = filename['Mean'][(len(filename.index) - 5):len(filename.index)].mean(axis=0)
    
    return threshold
    
def amygdala_loop(filename, timesbaseline):
    
    SliceList = pd.DataFrame()
    for n in np.arange(1, 8):
            if not os.path.isfile(filename % n):
                break
            Slice = pd.DataFrame(np.random.randn(1,4),
                columns=['Number', 'Area', 'Intensity', 'Density'])
            ArcCount = pd.read_csv(filename % n, index_col=" ")
            ArcCells = ArcCount['Mean'][(ArcCount['Mean'] > (threshold(ArcCount)*timesbaseline))]
            Slice['Area'][0] = ArcCount['Area'][1]
            Slice['Number'][0] = len(ArcCells.index) + 1
            Slice['Intensity'][0] = ArcCells.mean(axis=0) - threshold(ArcCount)
            # 44696.932 transforms to area/100um^2 based on .473um pixel size
            Slice['Density'][0] = Slice['Number'][0] / (Slice['Area'][0] / 44696.932)
            SliceList = SliceList.append(Slice)
            
    return SliceList
    
def make_tables(Groups, CombinedMice, MiceLeft, MiceRight, AntLeft, CombinedAnt, CombinedPost):
    GroupList = list(Groups)
    hc = [20,21,22,23,24,25,26,27]
    low = [13,14,15,16,17,18,19]
    med = [5,6,7,8,9,10,11,12]
    high = [0,1,2,3,4]
    
    AverageDensity = pd.DataFrame(        [
    
            CombinedMice['Density'][hc].tolist(),
            CombinedMice['Density'][low].tolist(),
            CombinedMice['Density'][med].tolist(),
            CombinedMice['Density'][high].tolist()],
            index=GroupList).transpose()
                    
    
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
            
    return  AverageDensity, LeftDensity, RightDensity, AntDensity, PostDensity
    """

def make_tables(CombinedMice, MiceLeft, MiceRight, AntLeft, CombinedAnt, CombinedPost):
    saline = [2, 5, 6, 7]
    cno = [0, 1, 3, 4, 8] 
    
    AverageDensity = pd.DataFrame(        [
            CombinedMice['Density'][saline].tolist(),
            CombinedMice['Density'][cno].tolist()],
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
            
    return  AverageDensity, LeftDensity, RightDensity, AntDensity, PostDensity
          
    """
def cell_count():
    
        # Instantiate all of the global Dataframes
    timesbaseline = 1.5
    MiceLeft = pd.DataFrame()
    MiceRight = pd.DataFrame()
    AntLeft = pd.DataFrame()
    PostLeft = pd.DataFrame()
    AntRight = pd.DataFrame()
    PostRight = pd.DataFrame()
    os.chdir('To Analyze')
    
    # Left
    i = 1
    # Loop through folders, building MiceLeft, AntLeft, and Post,Left
    while (i <= len(os.listdir())):
        os.chdir('%d' % i)
        Amygdala = pd.DataFrame()
        

        # Loop through slices, building up Amygdala
        Amygdala = amygdala_loop('%dlBACount.csv', timesbaseline)
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
        Amygdala = amygdala_loop('%drBACount.csv', timesbaseline)
        MiceRight = MiceRight.append(Amygdala.mean(axis=0), ignore_index=True)
        
        # Seperate anterior and posterior slices
        AntRight = AntRight.append(
            Amygdala[0:(int(len(Amygdala.index) / 2))].mean(axis=0), ignore_index=True)
        PostRight = PostRight.append(Amygdala[(int(len(Amygdala.index) / 2)):len(Amygdala.index)].mean(axis=0), ignore_index=True)
        os.chdir(os.pardir)
        i = i + 1

    CombinedMice = (MiceRight + MiceLeft) / 2
    CombinedAnt = (AntRight + AntLeft) / 2
    CombinedPost = (PostRight + PostLeft) / 2

    # Make data tables that can be copied into GraphPad
    return make_tables(CombinedMice, MiceLeft, MiceRight, AntLeft, CombinedAnt, CombinedPost)    

if __name__ == "__main__":

    Density, LeftDensity, RightDensity, AntDensity, PostDensity = cell_count()
    
    #Plan for group organization
    #Insert a new column into mice dataframe that includes treatment group, then
    #simply draw on values corresponding to treatment condition when making tables