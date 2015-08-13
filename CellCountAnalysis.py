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
    
def amygdala_loop(filename, NeuNCells, NeuNCount, timesbaseline):
    
    SliceList = pd.DataFrame()
    for n in np.arange(1, 8):
            if not os.path.isfile(filename % n):
                break
            Slice = pd.DataFrame(np.random.randn(1,6),
                columns=['Number', 'Area', 'Intensity', 'Density', 'Proportion', 'Cell Density'])
            ArcCount = pd.read_csv(filename % n, index_col=" ")
            ArcCells = ArcCount['Mean'][(ArcCount['Mean'] > (threshold(ArcCount)*timesbaseline))]
            Slice['Cell Density'][0] = (len(NeuNCells.index) + 1) / (NeuNCount['Area'][1])
            Slice['Area'][0] = ArcCount['Area'][1]
            Slice['Number'][0] = len(ArcCells.index) + 1
            Slice['Intensity'][0] = ArcCells.mean(axis=0) - threshold(ArcCount)
            # 44696.932 transforms to area/100um^2 based on .473um pixel size
            Slice['Density'][0] = Slice['Number'][0] / (Slice['Area'][0] / 44696.932)
            Slice['Proportion'][0] = Slice['Density'][0] * Slice['Cell Density'][0]
            SliceList = SliceList.append(Slice)
            
    return SliceList

def make_tables(CombinedMice, MiceLeft, MiceRight, AntLeft, CombinedAnt, CombinedPost):
    
    AverageDensity = pd.DataFrame(        [
            CombinedMice['Density'][0:4].tolist(),
            CombinedMice['Density'][19:25].tolist(),
            CombinedMice['Density'][10:19].tolist(),
            CombinedMice['Density'][4:10].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    AverageArea = pd.DataFrame(        [
            CombinedMice['Area'][0:4].tolist(),
            CombinedMice['Area'][19:25].tolist(),
            CombinedMice['Area'][10:19].tolist(),
            CombinedMice['Area'][4:10].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    AverageProportion = pd.DataFrame(        [
            CombinedMice['Proportion'][0:4].tolist(),
            CombinedMice['Proportion'][19:25].tolist(),
            CombinedMice['Proportion'][10:19].tolist(),
            CombinedMice['Proportion'][4:10].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
    
    LeftDensity = pd.DataFrame(        [
            MiceLeft['Density'][0:4].tolist(),
            MiceLeft['Density'][19:25].tolist(),
            MiceLeft['Density'][10:19].tolist(),
            MiceLeft['Density'][4:10].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    RightDensity = pd.DataFrame(        [
            MiceRight['Density'][0:4].tolist(),
            MiceRight['Density'][19:25].tolist(),
            MiceRight['Density'][10:19].tolist(),
            MiceRight['Density'][4:10].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    AntDensity = pd.DataFrame(        [
            CombinedAnt['Density'][0:4].tolist(),
            CombinedAnt['Density'][19:25].tolist(),
            CombinedAnt['Density'][10:19].tolist(),
            CombinedAnt['Density'][4:10].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    PostDensity = pd.DataFrame(        [
            CombinedPost['Density'][0:4].tolist(),
            CombinedPost['Density'][19:25].tolist(),
            CombinedPost['Density'][10:19].tolist(),
            CombinedPost['Density'][4:10].tolist()],
            index=['HC', '.3 mA', '.5 mA', '.75 mA']).transpose()
            
    return  AverageDensity, AverageProportion, LeftDensity, RightDensity, AverageArea, AntDensity, PostDensity
            
    
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
        NeuNCount = pd.read_csv('NeuNlCount.csv', index_col=" ")
        NeuNCells = NeuNCount['Mean'][NeuNCount['Mean'] > (threshold(NeuNCount) * timesbaseline)]

        # Loop through slices, building up Amygdala
        Amygdala = amygdala_loop('%dlCount.csv', NeuNCells, NeuNCount, timesbaseline)
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
        NeuNCount = pd.read_csv('NeuNrCount.csv', index_col=" ")
        NeuNCells = NeuNCount['Mean'][(NeuNCount['Mean'] > (threshold(NeuNCount) * timesbaseline))]

        # Loop through slices
        Amygdala = amygdala_loop('%drCount.csv', NeuNCells, NeuNCount, timesbaseline)
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

    Density, Proportion, LeftDensity, RightDensity, Area, AntDensity, PostDensity = cell_count()
    