# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 14:18:51 2022

@author: Yang
"""

"""
This is some tools for analyzing the simulation result from AtomECS (SimECS).
"""

import scipy.constants as cts
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.rc('text', usetex = True)
#plt.rcParams['font.family'] = "Times New Roman"
import CONST


def get_T_Rho(FileName, cap_r, r_max, dr):
    #define useful constants 
    m_Rb = CONST.M_RB
    Kb = CONST.KB
    
    # load the file and add labels to each columns of them
    trj = pd.read_table(FileName, skiprows=9, sep=' ', skipinitialspace=True)
    trj.columns = ['id','atom','x','y','z','vx','vy','vz','speed','vxy','t'] 

    # define the atom cloud by selecting atom under at certain threshold 
    trj_cloud = trj[(trj.x)**2 + (trj.y)**2 + (trj.z)**2 <= cap_r**2]

    # getting coordinates 
    vxs, vys, vzs = get_velocities(trj_cloud)

    # substracting the COM velocity
    vxs = vxs - np.mean(vxs)
    vys = vys - np.mean(vys)
    vzs = vzs - np.mean(vzs)

    # calculating the temperature of the atom cloud
    speeds2 = vxs**2 + vys**2 + vzs**2   
    T = ((0.5*m_Rb*speeds2)/(1.5*Kb)).mean()

    # calculate rhos 
    R, rho_shell, rho_e, psd_e, psd_mean, psd_max = get_Rho(trj, cap_r, r_max, dr, T)

    return (len(trj_cloud),T, R, np.array(rho_shell), rho_e, psd_e, psd_mean, psd_max)

def get_Rho (trj, cap_r, r_max, dr, T):
    m_Rb = CONST.M_RB
    Kb = CONST.KB
    h = CONST.H 

    # define the atom cloud by selecting atom under at certain threshold 
    trj_core = trj[(trj.x)**2 + (trj.y)**2 + (trj.z)**2 <= cap_r**2]

    # get the coordinates of each direction
    xs, ys, zs = get_coord(trj, use_com = True)

    # calculate the distance of each atoms to the origin
    r2 = xs**2 + ys**2 + zs**2
    # the r value for each atom for calculating the histogram
    r = np.sqrt(r2)
    rg = np.sqrt(np.mean(r**2))
    r_sphere = np.sqrt(5/3) * rg

    # calculate the number of bins for the RDF
    nbins = int(r_max/dr)
    Natom_shell, R = np.histogram(r, bins = nbins, range = (0, r_max), density = False)
    R = np.delete(R,0)

    #calculate the volume of each shell, or, in another words, the normalization factor
    norm_factors = (4/3)*np.pi*(R**3 - (R-dr)**3)
    rho_shell = (Natom_shell/norm_factors)

    # calculate different rhos
    e_radius = get_eradius(rho_shell, dr)
    rho_e = len(trj_core) / calc_volume(e_radius)
    rho_mean = len(trj_core) / calc_volume(r_sphere)
    rho_max = np.max(rho_shell)

    # calculate lambda 
    lamb_da = h/np.sqrt(2*np.pi*m_Rb*Kb*T)
    # calculate the PSD
    psd_e = rho_e * lamb_da**3
    psd_mean = rho_mean * lamb_da**3
    psd_max = rho_max * lamb_da**3

    return (np.array(R), np.array(rho_shell), rho_mean, psd_e, psd_mean, psd_max)


def get_coord(trj, use_com):
    if use_com == True:
        comx = np.mean(np.array(trj.iloc[:, 2]))
        comy = np.mean(np.array(trj.iloc[:, 3]))
        comz = np.mean(np.array(trj.iloc[:, 4]))
    else:
        comx = 0
        comy = 0
        comz = 0
    return (np.array(trj.iloc[:, 2]) - comx , np.array(trj.iloc[:, 3]) - comy, np.array(trj.iloc[:, 4]) - comz) 

def get_velocities(trj):
    return (trj.vx, trj.vy, trj.vz)

def calc_volume (r):
    return ((4/3) * np.pi * r**3)

def get_eradius (rho_shell, dr):
    e_rho_max = np.max(rho_shell)/np.e
    e_r_ndx = min(range(len(rho_shell)), key = lambda i: abs(rho_shell[i] - e_rho_max))
    return e_r_ndx*dr

def get_instant_laser_intersection(timestep, frequency, lx0 = 0.0, ly0 = 0.0, lz0 = 0.0):
    lx = 0.0002 * np.sin(frequency*2*np.pi * 1e-6 * timestep)
    ly = 0.0002 * np.sin(frequency*2*np.pi * 1e-6 * timestep)
    lz = 0.0002 * np.sin(frequency*2*np.pi * 1e-6 * timestep)
    return (lx, ly, lz)

def get_Ti(FileName):
    m_Rb = 86.909*cts.value('atomic mass constant')
    Kb = cts.value('Boltzmann constant')

    trj = pd.read_table(FileName, skiprows=9, sep=' ', skipinitialspace=True)
    trj.columns = ['id','atom','x','y','z','vx','vy','vz','speed','vxy'] 
    vxi = np.array(trj.iloc[:, 6])
    vyi = np.array(trj.iloc[:, 7])
    vzi = np.array(trj.iloc[:, 8])
    speedsi2 = vxi**2 + vyi**2 + vzi**2   
    Ti = ((0.5*m_Rb*speedsi2)/(1.5*Kb)).mean()
    return (round(Ti,3))

def trj_analysis (features, pre_directory, tot_steps, d_step, cap_r, r_max, dr, output_dir, dt, skipfirstframe):

  
    # Initiate lists for storing calculated data
    Tini = [] # initial temperature
    Nini = [] # initial number
    TFinals = [] # final temperature
    NFinals = [] # final number
    PSDFinalsE = [] # final psd with e_radius determined density
    PSDFinalsMean = [] # final psd with mean density
    PSDFinalsMax = [] # final psd wihm max density

    for feature, tot_step in zip(features, tot_steps):
        
        if feature == '':
            feature = '.'

        # defining the path
        directory = pre_directory + feature +'/trjs/' 
        
        STEP = []
        T = []
        Natom = []
        RHO = []
        RHOE = []
        PSDE = []
        PSDMEAN = []
        PSDMAX = []
        
        # initialize the step number based on skipping the first step or not
        if skipfirstframe == True:
            step = d_step
        else:
            step = 0

        
        if step == 0:
            filename = '1.trj' # because we don't have "0.trj"
        else:
            filename = str(step) + '.trj'
        counter = 0
        #print(step, tot_step)
        
        while step <= tot_step:
            # print(counter)
            num, temp, R, rho, rho_mean, psd_e, psd_mean, psd_max = get_T_Rho(directory + filename, cap_r, r_max, dr)  
            
            # accumulating rho values
            if len(RHO) == 0:
                RHO = rho
            else:
                RHO += rho
               
            # append values to the storing lists    
            STEP.append(step) 
            T.append(temp)
            Natom.append(num)
            RHOE.append(rho_mean)
            PSDE.append(psd_e)
            PSDMEAN.append(psd_mean)
            PSDMAX.append(psd_max)
            step += d_step
            filename = str(step) + '.trj'
            counter += 1
            
        RHO = RHO/counter
        NFinals.append(num)
        TFinals.append(temp)
        PSDFinalsE.append(psd_e)
        PSDFinalsMean.append(psd_mean)
        PSDFinalsMax.append(psd_max)

        
        pd.DataFrame(RHO, R*1000).to_csv(output_dir + '/Rho_' + feature + '.csv')
        pd.DataFrame(T, np.array(STEP)*dt).to_csv(output_dir + '/T_' + feature + '.csv')
        pd.DataFrame(Natom, np.array(STEP)*dt).to_csv(output_dir + '/N_' + feature + '.csv')
        pd.DataFrame(np.array(RHOE), np.array(STEP)*dt).to_csv(output_dir + '/rho_mean_' + feature + '.csv')
        pd.DataFrame(np.array(PSDE), np.array(STEP)*dt).to_csv(output_dir + '/psd_eradius_' + feature + '.csv')
        pd.DataFrame(np.array(PSDMEAN), np.array(STEP)*dt).to_csv(output_dir + '/psd_mean_' + feature + '.csv')
        pd.DataFrame(np.array(PSDMAX), np.array(STEP)*dt).to_csv(output_dir + '/psd_max_' + feature + '.csv')
        
        print(str(feature) + " : " + str(Natom[-1]) + " atoms left." + " Final T: " + str(T[-1]*1e6) + " uK."
              + "eradius rho: " + str(rho_mean) + ". eradius psd: " + str(psd_e), ". mean psd: " + str(psd_mean), ". max psd: " + str(psd_max) + ".")



