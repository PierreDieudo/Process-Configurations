
import numpy as np
import pandas as pd
import os
from PIL import Image
from Hub import Hub_Connector


''' General information here: 

    This configuration is part of a PFD Kamran provided with recycling loops in between membranes.

    Please refer to me (s1854031@ed.ac.uk) for any questions or issues
    
    Pierre
 '''
def show_image():  
    img = Image.open('Kamran_Config_Test.png') #image of the configuration - only considfering membranes 1,2 and 3 so far
    img.show()

#show_image() #show the image of the configuration - comment if unecessary

#-----------------------------------------#
#--------- User input parameters ---------#
#-----------------------------------------#

directory = 'C:\\Users\\s1854031\\Desktop\\' #input file path here.


Feed = {
    "Feed_Composition": [0.2, 0.6, 0.1, 0.1], # molar fraction
    "Feed_Flow": 1e4,                           # mol/s (PS: 1 mol/s = 3.6 kmol/h)
}


Membrane_1 = {
    "Name": 'Membrane_1',
    "Solving_Method": 'CO',                     # 'CC' or 'CO' - CC is for counter-current, CO is for co-current
    "Temperature": 25+273.15,                   # Kelvin
    "Pressure_Feed": 2,                         # bar
    "Pressure_Permeate": 0.2,                   # bar
    "Area": 1e3,                                # m2
    "Permeance": [10000, 500, 150, 100],        # GPU
    "Pressure_Drop": True,
    }

Membrane_2 = {
    "Name": 'Membrane_2',
    "Solving_Method": 'CO',                     # 'CC' or 'CO' - CC is for counter-current, CO is for co-current
    "Temperature": 25+273.15,                   # Kelvin    
    "Pressure_Feed": 2,                         # bar
    "Pressure_Permeate": 0.2,                   # bar
    "Area": 1e4,                                # m2
    "Permeance": [10000, 500, 150, 100],        # GPU
    "Pressure_Drop": True,
    }

Membrane_3 = {
    "Name": 'Membrane_3',
    "Solving_Method": 'CO',                     # 'CC' or 'CO' - CC is for counter-current, CO is for co-current
    "Temperature": 25+273.15,                   # Kelvin
    "Pressure_Feed": 2,                         # bar
    "Pressure_Permeate": 0.2,                   # bar
    "Area": 1e3,                                # m2
    "Permeance": [10000, 500, 150, 100],        # GPU
    "Pressure_Drop": True,
    }
 
Component_properties = {
    "Viscosity_param": ([0.0479,0.6112],[0.0466,3.8874],[0.0333,-0.23498],[0.0558,3.8970]), # Viscosity parameters for each component: slope and intercept for the viscosity correlation wiht temperature (in K) - from NIST
    "Molar_mass": [44.009, 28.0134, 18.01528, 31.999], # Molar mass of each component in kg/kmol"        
    }

Fibre_Dimensions = {
    "D_in" : 150 * 1e-6, # Inner diameter in m (from µm)
    "D_out" : 300 * 1e-6, # Outer diameter in m (from µm)
    }

def Run(Membrane):

    Membrane["Sweep_Flow"] = 0 # No sweep in this configuration
    Membrane["Sweep_Composition"] = [0] * len(Membrane_1["Permeance"])

    Export_to_mass_balance = Membrane, Component_properties, Fibre_Dimensions

    J = len(Membrane["Permeance"]) #number of components
     
    # No sweep in any of the membranes in this configuration
   
    results, profile = Hub_Connector(Export_to_mass_balance)
    Membrane["Retentate_Composition"],Membrane["Permeate_Composition"],Membrane["Retentate_Flow"],Membrane["Permeate_Flow"] = results

    print(f"Overall mass balance error of membrane {Membrane["Name"]}: Feed + Sweep  - Retentate - Permeate = {abs(Membrane["Feed_Flow"] + Membrane["Sweep_Flow"] - Membrane["Retentate_Flow"] - Membrane["Permeate_Flow"]):.3e}")
    
    if np.any(profile<0):
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
Placeholder_2={#Intermediade data storage for the recycling loop entering the second membrane used to check for convergence
    "Feed_Composition": [0] * J,
    "Feed_Flow": 0,              
    } 


for i in range(max_iter):

    print(f"Iteration {i+1}")

    if i ==0:
        Membrane_2["Feed_Composition"] = Feed["Feed_Composition"]
        Membrane_2["Feed_Flow"] = Feed["Feed_Flow"]
        
    results_2, profile_2 = Run(Membrane_2) # results "[0]: x", "[1]: y", "[2]: Q_ret", "[3]: Q_perm"

    Membrane_1["Feed_Composition"] = results_2[0] # Membrane 2 Retentate is feed to Membrane 1
    Membrane_1["Feed_Flow"] = results_2[-2]

    Membrane_3["Feed_Composition"] = results_2[1] # Membrane 2 Permeate is feed to Membrane 3
    Membrane_3["Feed_Flow"] = results_2[-1]

    results_1, profile_1 = Run(Membrane_1)
    results_3, profile_3 = Run(Membrane_3)

    Placeholder_2["Feed_Flow"] = results_1[3] + results_3[2] + Feed["Feed_Flow"] # Initial Feed + Membrane 1 Permeate + Membrane 3 Retentate flows
    for j in range(J):
        Placeholder_2["Feed_Composition"][j] = ( results_1[1][j] * results_1[3] + results_3[0][j] * results_3[2] + Feed["Feed_Composition"][j] * Feed["Feed_Flow"] ) / (Placeholder_2["Feed_Flow"]) #Component Mole Fractions


    #check for convergence
    if i > 0 and np.all(np.abs(np.array(Placeholder_2["Feed_Composition"]) - np.array(Membrane_2["Feed_Composition"]))) < tolerance and abs( ( (Placeholder_2["Feed_Flow"]) - (Membrane_2["Feed_Flow"] ) ) / (Membrane_2["Feed_Flow"] ) ) < tolerance:  
        print(f"Converged after {i} iterations")
        print()
        print("checking convergence for debugging:")
        print(f'Membrane 2 Feed Composition before iteration is {Placeholder_2["Feed_Composition"]} with flow {Placeholder_2["Feed_Flow"]:3e} mol/s')
        print(f'Membrane 2 Feed Composition after iteration is {Membrane_2["Feed_Composition"]} with flow {Membrane_2["Feed_Flow"]:3e} mol/s')
        print(f'Membrane 2 Profile')
        print(profile_2)
        print()
        print(f'Membrane 1 Profile')
        print(profile_1)
        print()
        print(f'Membrane 3 Profile')
        print(profile_3)
        break

    #Ready for next iteration
    Membrane_2["Feed_Composition"] = Placeholder_2["Feed_Composition"]
    Membrane_2["Feed_Flow"] = Placeholder_2["Feed_Flow"]

else: print("Max iterations reached")


print("Done - probably")
