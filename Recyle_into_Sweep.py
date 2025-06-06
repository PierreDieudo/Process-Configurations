

import numpy as np
import pandas as pd
import os
from PIL import Image
from Hub import Hub_Connector


''' General information here: 

    This configuration is for a simple one stage membrane, with Retentate recycling into the sweep.

    Please refer to me (s1854031@ed.ac.uk) for any questions or issues
    
    Pierre
 '''

def show_image():  
    img = Image.open('Recycle_into_sweep.png') #image of the configuration - only considfering membranes 1,2 and 3 so far
    img.show()

#show_image() #show the image of the configuration - comment if unecessary

#-----------------------------------------#
#--------- User input parameters ---------#
#-----------------------------------------#

directory = 'C:\\Users\\s1854031\\Desktop\\' #input file path here.

User_Sweep = {
    "Sweep_Composition": [0, 1, 0], # molar fraction
    "Sweep_Flow": 50,                           # mol/s (PS: 1 mol/s = 3.6 kmol/h)
}

Membrane_1 = {
    "Name": 'Membrane_1',
    "Solving_Method": 'CC',                     # 'CC' or 'CO' - CC is for counter-current, CO is for co-current
    "Temperature": 25+273.15,                   # Kelvin
    "Pressure_Feed": 2,                         # bar
    "Pressure_Permeate": 0.2,                   # bar
    "Area": 2000,                               # m2
    "Permeance": [500, 20, 15],                 # GPU
    "Pressure_Drop": True,
    "Recycling_Ratio": 0.3,                     # Ratio of Retentate flow recycled into the Sweep flow"
    
    "Feed_Composition": [0.3, 0.6, 0.1],        # molar fraction
    "Feed_Flow": 250,                           # mol/s (PS: 1 mol/s = 3.6 kmol/h)                   
    }
 
Component_properties = {
    "Viscosity_param": ([0.0479,0.6112],[0.0466,3.8874],[0.0558,3.8970]), # Viscosity parameters for each component: slope and intercept for the viscosity correlation wiht temperature (in K) - from NIST
    "Molar_mass": [44.009, 28.0134, 31.999], # Molar mass of each component in kg/kmol"        
    }

Fibre_Dimensions = {
    "D_in" : 150 * 1e-6, # Inner diameter in m (from µm)
    "D_out" : 300 * 1e-6, # Outer diameter in m (from µm)
    }

J = len(Membrane_1["Permeance"]) #number of components

def Run(Membrane):

    Export_to_mass_balance = Membrane, Component_properties, Fibre_Dimensions

    results, profile = Hub_Connector(Export_to_mass_balance)
    Membrane["Retentate_Composition"],Membrane["Permeate_Composition"],Membrane["Retentate_Flow"],Membrane["Permeate_Flow"] = results

    print(f"Overall mass balance error of membrane {Membrane["Name"]}: Feed + Sweep  - Retentate - Permeate = {abs(Membrane["Feed_Flow"] + Membrane["Sweep_Flow"] - Membrane["Retentate_Flow"] - Membrane["Permeate_Flow"]):.3e}")
    
    if np.any(profile<-1e-6):
       print(profile)
       raise ValueError("Negative values in the membrane profile") #check for negative values in the profile

    #Reformat Permeance and Pressure values to the initial units
    Membrane["Permeance"] = [p / ( 3.348 * 1e-10 ) for p in Membrane["Permeance"]]  # convert from mol/m2.s.Pa to GPU
    Membrane["Pressure_Feed"] *= 1e-5  #convert to bar
    Membrane["Pressure_Permeate"] *= 1e-5  

    #print(profile)
    return results, profile

### Run iterations for process recycling loops - Specific to this configuration!

max_iter = 300
tolerance = 1e-6
J = len(Membrane_1["Permeance"]) #number of components
Placeholder_1={ #Intermediade data storage for the recycling loop entering the second membrane used to check for convergence
    "Sweep_Composition": [0] * J,
    "Sweep_Flow": 0,              
    } 


for i in range(max_iter):

    print(f"Iteration {i+1}")

    if i ==0:
        Membrane_1["Sweep_Composition"] = User_Sweep["Sweep_Composition"]
        Membrane_1["Sweep_Flow"] = User_Sweep["Sweep_Flow"]
        
    results_1, profile_1 = Run(Membrane_1) # results "[0]: x", "[1]: y", "[2]: Q_ret", "[3]: Q_perm"


    Placeholder_1["Sweep_Flow"] = results_1[2] * Membrane_1["Recycling_Ratio"] + User_Sweep["Sweep_Flow"] # Fraction of the Retentate Recycled + Initial Sweep
    for j in range(J):
        Placeholder_1["Sweep_Composition"][j] = ( results_1[0][j] * results_1[2] * Membrane_1["Recycling_Ratio"] + User_Sweep["Sweep_Composition"][j] * User_Sweep["Sweep_Flow"] ) / (Placeholder_1["Sweep_Flow"]) #Component Mole Fractions

    Convergence_Composition = sum(abs(np.array(Placeholder_1["Sweep_Composition"]) - np.array(Membrane_1["Sweep_Composition"])))
    Convergence_Flowrate = abs( ( (Placeholder_1["Sweep_Flow"]) - (Membrane_1["Sweep_Flow"] ) ) / (Membrane_1["Sweep_Flow"] ) / 100 )
    print(f'Convergence Composition: {Convergence_Composition:.3e}, Convergence Flowrate: {Convergence_Flowrate*100:.3e}')
    print()

    #check for convergence
    if i > 0 and Convergence_Composition < tolerance and Convergence_Flowrate < tolerance:  
        print(f"Converged after {i+1} iterations")
        print()
        print("checking convergence for debugging:")
        print(f'Membrane 1 Sweep Composition after iteration is {[f"{comp:.4g}" for comp in Placeholder_1["Sweep_Composition"]]} with flow {Placeholder_1["Sweep_Flow"]:.4g} mol/s')
        print(f'Membrane 1 Sweep Composition before iteration is {[f"{comp:.4g}" for comp in Membrane_1["Sweep_Composition"]]} with flow {Membrane_1["Sweep_Flow"]:.4g} mol/s')
        print()
        print(f'Component 1 final Recovery is {results_1[1][0]*results_1[3]/(Membrane_1["Feed_Composition"][0]*Membrane_1["Feed_Flow"])*100:.4f}%')
        print(f'Component 1 final Purity is {results_1[1][0]*100:.4f}%')
        print()
        print(f'Membrane 1 Profile')
        print(profile_1)
        print()

        break
    #Ready for next iteration
    Membrane_1["Sweep_Composition"] = Placeholder_1["Sweep_Composition"]
    Membrane_1["Sweep_Flow"] = Placeholder_1["Sweep_Flow"]

else: print("Max iterations reached")


print("Done - probably")
