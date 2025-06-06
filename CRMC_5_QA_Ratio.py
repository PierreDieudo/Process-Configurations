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
    "Solving_Method": 'CC',                     # 'CC' or 'CO' - CC is for counter-current, CO is for co-current
    "Temperature": 25+273.15,                   # Kelvin
    "Pressure_Feed": 2,                         # bar
    "Pressure_Permeate": 0.2,                   # bar
    "Q_A_ratio": 7.8,                             # ratio of the membrane feed flowrate to its area (in m3(stp)/m2.hr)
    "Permeance": [10000, 500, 150, 100],        # GPU
    "Pressure_Drop": True,
    }

Membrane_2 = {
    "Name": 'Membrane_2',
    "Solving_Method": 'CC',                   
    "Temperature": 25+273.15,                   
    "Pressure_Feed": 2,                       
    "Pressure_Permeate": 0.2,                  
    "Q_A_ratio": 8.5,                             # ratio of the membrane feed flowrate to its area (in m3(stp)/m2.hr)
    "Permeance": [10000, 500, 150, 100],        
    "Pressure_Drop": True,
    }

Membrane_3 = {
    "Name": 'Membrane_3',
    "Solving_Method": 'CC',                     
    "Temperature": 25+273.15,                    
    "Pressure_Feed": 2,                          
    "Pressure_Permeate": 0.2,                    
    "Q_A_ratio": 13,                             # ratio of the membrane feed flowrate to its area (in m3(stp)/m2.hr)
    "Permeance": [10000, 500, 150, 100],        
    "Pressure_Drop": True,
    }

Membrane_4 = {
    "Name": 'Membrane_4',
    "Solving_Method": 'CC',                     # 'CC' or 'CO' - CC is for counter-current, CO is for co-current
    "Temperature": 25+273.15,                   # Kelvin
    "Pressure_Feed": 2,                         # bar
    "Pressure_Permeate": 0.2,                   # bar
    "Q_A_ratio": 23,                             # ratio of the membrane feed flowrate to its area (in m3(stp)/m2.hr)
    "Permeance": [10000, 500, 150, 100],        # GPU
    "Pressure_Drop": True,
    }

Membrane_5 = {
    "Name": 'Membrane_5',
    "Solving_Method": 'CC',                     # 'CC' or 'CO' - CC is for counter-current, CO is for co-current
    "Temperature": 25+273.15,                   # Kelvin
    "Pressure_Feed": 2,                         # bar
    "Pressure_Permeate": 0.2,                   # bar
    "Q_A_ratio": 42,                             # ratio of the membrane feed flowrate to its area (in m3(stp)/m2.hr)
    "Permeance": [10000, 500, 150, 100],        # GPU
    "Pressure_Drop": True,
    }

Component_properties = {
    "Viscosity_param": ([0.0479,0.6112],[0.0466,3.8874],[0.0333,-0.23498],[0.0558,3.8970]), # Viscosity parameters for each component: slope and intercept for the viscosity correlation wiht temperature (in K) - from NIST
    "Molar_mass": [44.009, 28.0134, 18.01528, 31.999], # Molar mass of each component in kg/kmol"        
    }

Fibre_Dimensions = {
    "D_in" : 150 * 1e-6, # Inner diameter in m (from um)
    "D_out" : 300 * 1e-6, # Outer diameter in m (from um)
    }

def Run(Membrane):

    # Set membrane Area based on its feed flow and Q_A_ratio:
    Membrane["Area"] = (Membrane["Feed_Flow"] * 0.0224  * 3600) / Membrane["Q_A_ratio"] # (0.0224 is the molar volume of an ideal gas at STP in m3/mol)
    print(f'{Membrane["Name"]} Area: {Membrane["Area"]:.0f} m2')

    Membrane["Sweep_Flow"] = 0 # No sweep in this configuration
    Membrane["Sweep_Composition"] = [0] * len(Membrane_1["Permeance"])

    Export_to_mass_balance = Membrane, Component_properties, Fibre_Dimensions

    J = len(Membrane["Permeance"]) #number of components
     
    # No sweep in any of the membranes in this configuration
    results, profile = Hub_Connector(Export_to_mass_balance)
    Membrane["Retentate_Composition"],Membrane["Permeate_Composition"],Membrane["Retentate_Flow"],Membrane["Permeate_Flow"] = results

    #print(f"Overall mass balance error of membrane {Membrane["Name"]}: Feed + Sweep  - Retentate - Permeate = {abs(Membrane["Feed_Flow"] + Membrane["Sweep_Flow"] - Membrane["Retentate_Flow"] - Membrane["Permeate_Flow"]):.3e}")
    
    if np.any(profile<-1e-4):
       print(f"{Membrane["Name"]} profile")
       print(f'{Membrane["Name"]} feed flow: {Membrane["Feed_Flow"]:.3f}')
       print(profile)
       raise ValueError(f"Negative values in the {["Name"]} profile") #check for negative values in the profile

    #Reformat Permeance and Pressure values to the initial units
    Membrane["Permeance"] = [p / ( 3.348 * 1e-10 ) for p in Membrane["Permeance"]]  # convert from mol/m2.s.Pa to GPU
    Membrane["Pressure_Feed"] *= 1e-5  #convert to bar
    Membrane["Pressure_Permeate"] *= 1e-5  

    #print(profile)
    return results, profile

### Run iterations for process recycling loops - Specific to this configuration!

max_iter = 300
tolerance = 1-5
J = len(Membrane_1["Permeance"]) #number of components
Placeholder_2={#Intermediade data storage for the recycling loop entering the second membrane used to check for convergence
    "Feed_Composition": [0] * J,
    "Feed_Flow": 0,              
    } 
Placeholder_3={#Intermediade data storage for the recycling loop entering the second membrane used to check for convergence
    "Feed_Composition": [0] * J,
    "Feed_Flow": 0,              
    } 
Placeholder_4={#Intermediade data storage for the recycling loop entering the second membrane used to check for convergence
    "Feed_Composition": [0] * J,
    "Feed_Flow": 0,              
    } 

for i in range(max_iter):

    print(f"Iteration {i+1}")

    if i ==0:
        Membrane_2["Feed_Composition"] = Feed["Feed_Composition"]
        Membrane_2["Feed_Flow"] = Feed["Feed_Flow"]
        
    results_2, profile_2 = Run(Membrane_2) # results "[0]: x", "[1]: y", "[2]: Q_ret", "[3]: Q_perm"

    if i==0: #First iteration, no Placeholder_2 or Placeholder_3
        Membrane_3["Feed_Composition"] = results_2[1] # Membrane 2 Permeate is feed to Membrane 3
        Membrane_3["Feed_Flow"] = results_2[-1]
 
    Membrane_1["Feed_Composition"] = results_2[0] # Membrane 2 Retentate is feed to Membrane 1
    Membrane_1["Feed_Flow"] = results_2[-2]

    results_1, profile_1 = Run(Membrane_1)
    results_3, profile_3 = Run(Membrane_3)

    Placeholder_2["Feed_Flow"] = results_1[3] + results_3[2] + Feed["Feed_Flow"] # Initial Feed + Membrane 1 Permeate + Membrane 3 Retentate flows
    for j in range(J):
        Placeholder_2["Feed_Composition"][j] = ( results_1[1][j] * results_1[3] + results_3[0][j] * results_3[2] + Feed["Feed_Composition"][j] * Feed["Feed_Flow"] ) / (Placeholder_2["Feed_Flow"]) #Component Mole Fractions
    
    if i == 0: #First iteration, no Placeholder_4
        Membrane_4["Feed_Composition"]= results_3[1]
        Membrane_4["Feed_Flow"]= results_3[-1]

    results_4, profile_4 = Run(Membrane_4)

    Placeholder_3["Feed_Flow"] = results_4[2] + results_2[3] # Membrane 2 Permeate + Membrane 4 Retentate flows
    for j in range(J):
        Placeholder_3["Feed_Composition"][j] = ( results_4[0][j] * results_4[2] + results_2[1][j] * results_2[3]) / (Placeholder_3["Feed_Flow"]) #Component Mole Fractions


    Membrane_5["Feed_Composition"]= results_4[1]
    Membrane_5["Feed_Flow"]= results_4[-1]

    results_5, profile_5 = Run(Membrane_5)
    Placeholder_4["Feed_Flow"] = results_3[3] + results_5[2] # Membrane 4 Permeate + Membrane 5 Retentate flows
    for j in range(J):
        Placeholder_4["Feed_Composition"][j] = ( results_5[0][j] * results_5[2] + results_3[1][j] * results_3[3]) / (Placeholder_4["Feed_Flow"])
    

    Convergence_Composition3 = sum(abs(np.array(Placeholder_3["Feed_Composition"]) - np.array(Membrane_3["Feed_Composition"]))) #Convergence check for compositions
    Convergence_Composition4 = sum(abs(np.array(Placeholder_4["Feed_Composition"]) - np.array(Membrane_4["Feed_Composition"])))
    Convergence_Composition2 = sum(abs(np.array(Placeholder_2["Feed_Composition"]) - np.array(Membrane_2["Feed_Composition"])))
    Convergence_Composition = (Convergence_Composition2 + Convergence_Composition3 + Convergence_Composition4) 

    Convergence_Flowrate3 = abs( ( (Placeholder_3["Feed_Flow"]) - (Membrane_3["Feed_Flow"] ) ) / (Membrane_3["Feed_Flow"] ))
    Convergence_Flowrate4 = abs( ( (Placeholder_4["Feed_Flow"]) - (Membrane_4["Feed_Flow"] ) ) / (Membrane_4["Feed_Flow"] ))
    Convergence_Flowrate2 = abs( ( (Placeholder_2["Feed_Flow"]) - (Membrane_2["Feed_Flow"] ) ) / (Membrane_2["Feed_Flow"] ))
    Convergence_Flowrate = (Convergence_Flowrate2 + Convergence_Flowrate3 + Convergence_Flowrate4) / 100

    print(f'Convergence Composition: {Convergence_Composition:.3e}, Convergence Flowrate: {Convergence_Flowrate*100:.3e}')
    #check for convergence
    if i > 0 and Convergence_Composition < tolerance and Convergence_Flowrate < tolerance:  
        print(f"Converged after {i+1} iterations")
        print()
        print("checking convergence for debugging:")
        print(f'Membrane 3 Feed Composition before iteration is {Placeholder_3["Feed_Composition"]} with flow {Placeholder_3["Feed_Flow"]:3e} mol/s')
        print(f'Membrane 3 Feed Composition after iteration is {Membrane_3["Feed_Composition"]} with flow {Membrane_3["Feed_Flow"]:3e} mol/s')
        print(f'Membrane 4 Feed Composition before iteration is {Placeholder_4["Feed_Composition"]} with flow {Placeholder_4["Feed_Flow"]:3e} mol/s')
        print(f'Membrane 4 Feed Composition after iteration is {Membrane_4["Feed_Composition"]} with flow {Membrane_4["Feed_Flow"]:3e} mol/s')
        print(f'Membrane 2 Feed Composition before iteration is {Placeholder_2["Feed_Composition"]} with flow {Placeholder_2["Feed_Flow"]:3e} mol/s')
        print(f'Membrane 2 Feed Composition after iteration is {Membrane_2["Feed_Composition"]} with flow {Membrane_2["Feed_Flow"]:3e} mol/s')        
        print(f'Component 1 Recovery in P is {results_5[1][0]*results_5[3]/(Feed["Feed_Composition"][0]*Feed["Feed_Flow"])*100:.4f}%')
        print(f'Component 1 Purity in P is {results_5[1][0]*100:.4f}%')
        
        print(f'Membrane 2 Profile')
        print(profile_2)
        print(f'Membrane 2 retentate and permeate flow are: {results_2[2]:.3f} and {results_2[3]:.3f} with feed flow {Membrane_2["Feed_Flow"]:.3f} mol/s')
        print()
        print(f'Membrane 1 Profile')
        print(profile_1)
        print(f'Membrane 1 retentate and permeate flow are: {results_1[2]:.3f} and {results_1[3]:.3f} with feed flow {Membrane_1["Feed_Flow"]:.3f} mol/s')
        print()
        print(f'Membrane 3 Profile')
        print(profile_3)
        print(f'Membrane 3 retentate and permeate flow are: {results_3[2]:.3f} and {results_3[3]:.3f} with feed flow {Membrane_3["Feed_Flow"]:.3f} mol/s')
        print()
        print(f'Membrane 4 Profile')
        print(profile_4)
        print(f'Membrane 4 retentate and permeate flow are: {results_4[2]:.3f} and {results_4[3]:.3f} with feed flow {Membrane_4["Feed_Flow"]:.3f} mol/s')
        print()
        print(f'Membrane 5 Profile')
        print(profile_5)
        print(f'Membrane 5 retentate and permeate flow are: {results_5[2]:.3f} and {results_5[3]:.3f} with feed flow {Membrane_5["Feed_Flow"]:.3f} mol/s')


        break

    #Ready for next iteration
    Membrane_2["Feed_Composition"] = Placeholder_2["Feed_Composition"]
    Membrane_2["Feed_Flow"] = Placeholder_2["Feed_Flow"]
    Membrane_3["Feed_Composition"] = Placeholder_3["Feed_Composition"]
    Membrane_3["Feed_Flow"] = Placeholder_3["Feed_Flow"]
    Membrane_4["Feed_Composition"] = Placeholder_4["Feed_Composition"]
    Membrane_4["Feed_Flow"] = Placeholder_4["Feed_Flow"]

else: print("Max iterations reached")


print("Done - probably")

