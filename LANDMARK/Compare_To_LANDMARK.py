# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:15:32 2024

@author: brickd
"""

import numpy as np 
import pandas as pd 
import os
import matplotlib.pyplot as plt


"""
Load the Landmark data
"""
Landmark_base_path = "C:/Users/brickd/Documents/Grad_classes/Scientific_Computing/Project/HypiC/LANDMARK/Case_2"
Base_E = pd.read_csv(os.path.join(Landmark_base_path, "electric_field_hybrid.csv"), header=0).to_numpy(np.float64)
Base_energy = pd.read_csv(os.path.join(Landmark_base_path, "energy_hybrid.csv"), header=0).to_numpy(np.float64)
Base_ionization = pd.read_csv(os.path.join(Landmark_base_path, "ionization_hybrid.csv"), header=0).to_numpy(np.float64)
Base_nn = pd.read_csv(os.path.join(Landmark_base_path, "neutral_density_hybrid.csv"), header=0).to_numpy(np.float64)
Base_ne = pd.read_csv(os.path.join(Landmark_base_path, "plasma_density_hybrid.csv"), header=0).to_numpy(np.float64)
Base_phi = pd.read_csv(os.path.join(Landmark_base_path, "potential_hybrid.csv"), header=0).to_numpy(np.float64)



"""
Load the simulation 
"""
Simulation_path = "C:/Users/brickd/Documents/Grad_classes/Scientific_Computing/Project/Build/test/Test_Output.csv"
Model_Results = pd.read_csv(Simulation_path, index_col=False).to_numpy(np.float64)


"""
Make the plots
"""
save_path = "C:/Users/brickd/Documents/Grad_classes/Scientific_Computing/Project/Build/test/"
dpi=300

#Electric Field
plt.figure()
plt.plot(Base_E[:,0], Base_E[:,1], 'k')
plt.plot(Model_Results[:,0], Model_Results[:,10], 'r')
plt.xlabel("Axial Position, meters")
plt.ylabel("Electric Field, V/m")
plt.savefig(os.path.join(save_path, "Electric_Field.png"), dpi=dpi, bbox_inches='tight') 

#Energy
plt.figure()
plt.plot(Base_energy[:,0], Base_energy[:,1], 'k')
plt.plot(Model_Results[:,0], 1.5 * Model_Results[:,8] * Model_Results[:,4], 'r')
plt.xlabel("Axial Position, meters")
plt.ylabel("Electron Energy Density, eV/m^3")
plt.savefig(os.path.join(save_path, "Energy.png"), dpi=dpi, bbox_inches='tight') 

#Ionization Rate
plt.figure()
plt.plot(Base_ionization[:,0], Base_ionization[:,1], 'k')
plt.plot(Model_Results[:,0], Model_Results[:,12], 'r')
plt.xlabel("Axial Position, meters")
plt.ylabel("Ionization Rate, 1/m^3s")
plt.savefig(os.path.join(save_path, "Ionization_Rate.png"), dpi=dpi, bbox_inches='tight') 

#Neutral Number Density
plt.figure()
plt.plot(Base_nn[:,0], Base_nn[:,1], 'k')
plt.plot(Model_Results[:,0], Model_Results[:,1], 'r')
plt.xlabel("Axial Position, meters")
plt.ylabel("Neutral Number Density, 1/m^3")
plt.savefig(os.path.join(save_path, "nn.png"), dpi=dpi, bbox_inches='tight') 

#Electron number density
plt.figure()
plt.plot(Base_ne[:,0], Base_ne[:,1], 'k')
plt.plot(Model_Results[:,0], Model_Results[:,4], 'r')
plt.xlabel("Axial Position, meters")
plt.ylabel("Electron Number Density, 1/m^3")
plt.savefig(os.path.join(save_path, "ne.png"), dpi=dpi, bbox_inches='tight') 

#Potential
plt.figure()
plt.plot(Base_phi[:,0], Base_phi[:,1], 'k')
plt.plot(Model_Results[:,0], Model_Results[:,13], 'r')
plt.xlabel("Axial Position, meters")
plt.ylabel("Potential, V")
plt.savefig(os.path.join(save_path, "Potential.png"), dpi=dpi, bbox_inches='tight') 