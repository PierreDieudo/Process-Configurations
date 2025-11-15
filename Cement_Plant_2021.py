

import numpy as np
import os
from Hub import Hub_Connector
from UNISIMConnect import UNISIMConnector
import pandas as pd  
import matplotlib.pyplot as plt

''' General information here: 

    This configuration is from a 2021 publication from Ferrari, et al. https://doi.org/10.3390/en14164772

    Please refer to me (s1854031@ed.ac.uk) for any questions or issues
    
    Pierre
 '''


#-----------------------------------------#
#--------- User input parameters ---------#
#-----------------------------------------#

filename = 'Cement_4Comp_FerrariPaper_Flash.usc' #Unisim file name
directory = 'C:\\Users\\s1854031\\OneDrive - University of Edinburgh\\Python\\Cement_Plant_2021\\' #Directory of the unisim file

unisim_path = os.path.join(directory, filename)

Options = {
    "Plot_Profiles" : True,                     # Plots the profiles of membranes 1 and 2 once the process is solved
    "Export_Profiles": False,                   # Exports membrane profiles into a csv file
    "Permeance_From_Activation_Energy": True    # True will use the activation energies from the component_properties dictionary - False will use the permeances defined in the membranes dictionaries.
    }


Membrane_1 = {
    "Name": 'Membrane_1',
    "Solving_Method": 'CO',                 # 'CC' or 'CO' - CC is for counter-current, CO is for co-current
    "Temperature": 35+273.15,               # Kelvin
    "Pressure_Feed": 7.7,                  # bar
    "Pressure_Permeate": 1,                 # bar
    "Q_A_ratio": 20,                      # ratio of the membrane feed flowrate to its area (in m3(stp)/m2.hr)
    "Permeance": [360, 13, 60, 360],        # GPU
    "Pressure_Drop": False,
    }

Membrane_2 = {
    "Name": 'Membrane_2',
    "Solving_Method": 'CO',                   
    "Temperature": 35+273.15,                   
    "Pressure_Feed": 9,                       
    "Pressure_Permeate": 1,                  
    "Q_A_ratio": 15,                          
    "Permeance": [360, 13, 60, 360],        
    "Pressure_Drop": False,
    }

Process_param = {
"Recycling_Ratio" : 1,      # Ratio of the retentate flow from Membrane 2 that is recycled back to Membrane 1 feed    
"Target_Purity" : 0.95,     # Target purity of the dry permeate from Membrane 2
"Target_Recovery" : 0.9,    # Target recovery from Membrane 2 - for now not a hard limit, but a target to be achieved
"Replacement_rate": 4,      # Replacement rate of the membranes (in yr)
"Operating_hours": 8000,    # Operating hours per year
"Lifetime": 20,             # Lifetime of the plant (in yr)
"Base Plant Cost": 149.8 * 1e6,     # Total direct cost of plant (no CCS) in 2014 money
"Contingency": 0.3,         # or 0.4 (30% or 40% contingency for process design - based on TRL)
}

Component_properties = {
    "Viscosity_param": ([0.0479,0.6112],[0.0466,3.8874],[0.0558,3.8970], [0.03333, -0.23498]),  # Viscosity parameters for each component: slope and intercept for the viscosity correlation wiht temperature (in K) - from NIST
    "Molar_mass": [44.009, 28.0134, 31.999,18.01528],                                           # Molar mass of each component in g/mol"
    "Activation_Energy_Aged": ([12750,321019],[25310,2186946],[15770,196980],[12750,321019]),   # ([Activation energy - J/mol],[pre-exponential factor - GPU])
    "Activation_Energy_Fresh": ([2880,16806],[16520,226481],[3770,3599],[2880,16806]),          #Valid only if temperature is -20 C or under - not considered for now
    }


Fibre_Dimensions = {
"D_in" : 150 * 1e-6,    # Inner diameter in m (from um)
"D_out" : 300 * 1e-6,   # Outer diameter in m (from um)
}
  
J = len(Membrane_1["Permeance"]) #number of components



with UNISIMConnector(unisim_path, close_on_completion=False) as unisim:

    class ConvergenceError(Exception): #allowing to skip iteration if convergence error appears in one of the mass balance
        pass
    Feed = {}
    Flue_Inlet = unisim.get_spreadsheet('Flue_Inlet')
    Feed['Feed_Flow'] = Flue_Inlet.get_cell_value('C3') / 3.6 # feed flow rate from UNISIM in mol/s (from kmol/h)
    Feed["Feed_Composition"] = [Flue_Inlet.get_cell_value(f'C{i+4}') for i in range(J)] #feed mole fractions from UNISIM

    #----------------------------------------#
    #----------- Get spreadsheets -----------#
    #----------------------------------------#

    Membrane1 = unisim.get_spreadsheet('Membrane 1')
    Membrane2 = unisim.get_spreadsheet('Membrane 2')
    Recycle_Membrane_2 = unisim.get_spreadsheet('Recycle Membrane 2')
    Duties = unisim.get_spreadsheet('Duties')


    #-----------------------------------------#
    #------------- Initial setup -------------#
    #-----------------------------------------#

    # Reset streams
    for i in range(J):
        Membrane1.set_cell_value(f'D{i+14}', 0) # Reset Membrane 1 permeate component flows
        Membrane1.set_cell_value(f'D{i+21}', 0) # Reset Membrane 1 retentate component flows
        Membrane2.set_cell_value(f'D{i+14}', 0) # Reset Membrane 2 permeate component flows
        Membrane2.set_cell_value(f'D{i+21}', 0) # Reset Membrane 2 retentate component flows
    
    unisim.wait_solution(timeout=10, check_pop_ups=2, check_consistency_error=3)

    # Setup Recycling Ratio
    Recycle_Membrane_2.set_cell_value('C3', Process_param["Recycling_Ratio"] ) # Set recycling ratio in the spreadsheet

    # Setup temperatures and pressures
    Membrane1.set_cell_value('D3', Membrane_1["Temperature"])  # Set temperature in Kelvin
    Membrane1.set_cell_value('D4', Membrane_1["Pressure_Feed"])  # Set feed pressure in bar
    Membrane1.set_cell_value('D6', Membrane_1["Pressure_Permeate"])  # Set permeate pressure in bar
    Membrane2.set_cell_value('D3', Membrane_2["Temperature"]) 
    Membrane2.set_cell_value('D4', Membrane_2["Pressure_Feed"]) 
    Membrane2.set_cell_value('D6', Membrane_2["Pressure_Permeate"])  # Set permeate pressure in bar

    unisim.wait_solution(timeout=10, check_pop_ups=2, check_consistency_error=3)
        
        
    #------------------------------------------------#
    #--------- Get correct feed from Unisim ---------#
    #------------------------------------------------#
    '''Function to select the correct membrane compression train - Determined to be the one with the lowest non-null compressor duty'''

    def Mem_Train_Choice(Membrane):
        
        unisim.wait_solution()

        Train_data = []
        for i in range(3):  
            if Membrane == Membrane_1:
                Train_data.append([  
                    Duties.get_cell_value(f'H{i+9}'),  # Compressor Duty (kW)
                    Duties.get_cell_value(f'I{i+9}'),  # Hex Area (m2)
                    Duties.get_cell_value(f'J{i+9}'),  # Water Flowrate (kg/hr)
                    Duties.get_cell_value(f'K{i+9}')/1e6   # Cryogenic Cooler Duty (GJ/hr)
                ])
            elif Membrane == Membrane_2:
                Train_data.append([  
                    Duties.get_cell_value(f'H{i+15}'),  # Compressor Duty (kW)
                    Duties.get_cell_value(f'I{i+15}'),  # Hex Area (m2)
                    Duties.get_cell_value(f'J{i+15}'),   # Water Flowrate (kg/hr)
                    Duties.get_cell_value(f'K{i+15}')/1e6   # Cryogenic Cooler Duty (GJ/hr)
                ])

            else: raise ValueError ("Incorrect membrane denomination")

        # Filter out trains with None or non-positive compressor duty
        valid_train_indices = [i for i, train in enumerate(Train_data) if train[0] is not None and train[0] > 0 and train[1] is not None and train[1] >0]
        if not valid_train_indices:
            raise ValueError("No valid trains found with positive compressor duty.")
        
        # Find the index of the train with the lowest compressor duty in Train_data
        lowest_duty_train_index = min(valid_train_indices, key=lambda i: Train_data[i][0])
        lowest_duty_train = Train_data[lowest_duty_train_index]

        # Read the spreadsheet of the corresponding train  
        if Membrane == Membrane_1:
            Mem_train = unisim.get_spreadsheet(f'Train 10{lowest_duty_train_index + 1}')  
            Membrane_1['Feed_Flow'] = Mem_train.get_cell_value('C3') / 3.6 # feed flow rate from UNISIM in mol/s (from kmol/h)
            Membrane_1["Feed_Composition"] = [Mem_train.get_cell_value(f'C{i+4}') for i in range(J)] #feed mole fractions from UNISIM
            Membrane_1["Train_Data"] = lowest_duty_train # Store the train data in the membrane dictionary
        elif Membrane == Membrane_2:
            Mem_train = unisim.get_spreadsheet(f'Train 20{lowest_duty_train_index + 1}')
            Membrane_2['Feed_Flow'] = Mem_train.get_cell_value('C3') / 3.6 # feed flow rate from UNISIM in mol/s (from kmol/h)
            Membrane_2["Feed_Composition"] = [Mem_train.get_cell_value(f'C{i+4}') for i in range(J)] #feed mole fractions from UNISIM
            Membrane_2["Train_Data"] = lowest_duty_train # Store the train data in the membrane dictionary

    #------------------------------------------#
    #--------- Function to run module ---------#
    #------------------------------------------#
    

    def plot_composition_profiles(profile,name):  

        df = profile.copy()
        global J
        z = df["norm_z"]

        fig1, axes1 = plt.subplots(1, 2, figsize=(16, 5))  

        # Retentate composition plot  
        for j in range(J):  
            col = f'x{j+1}'  
            if col in profile.columns:  
                axes1[0].plot(z, profile[col] * 100, label=f'Component {j+1}')  
        axes1[0].set_xlabel('Normalised Length')  
        axes1[0].set_ylabel('Retentate Composition (%)')  
        axes1[0].set_title(f'{name} Retentate Composition Profile')  
        axes1[0].legend()  
        axes1[0].grid(True)  

        # Permeate composition plot  
        for j in range(J):  
            col = f'y{j+1}'  
            if col in profile.columns:  
                axes1[1].plot(z, profile[col] * 100, label=f'Component {j+1}')  
        axes1[1].set_xlabel('Normalised Length')  
        axes1[1].set_ylabel('Permeate Composition (%)')  
        axes1[1].set_title(f'{name} Permeate Composition Profile')  
        axes1[1].legend()  
        axes1[1].grid(True)  

        plt.tight_layout()  
        plt.show(block=False)

        fig2, axes2 = plt.subplots(1, 2, figsize=(16, 5))  

        # Retentate component flows plot
        for j in range(J):  
            col = f'x{j+1}'  
            if col in profile.columns:  
                # Multiply component fraction by normalised retentate flow
                axes2[0].plot(z, profile[col] * profile['cut_r/Qr'], label=f'Component {j+1}')  

        axes2[0].set_xlabel('Normalised Length')  
        axes2[0].set_ylabel('Retentate Normalised Component Flow (-)')  
        axes2[0].set_title(f'{name} Retentate Flow Profile')  
        axes2[0].legend()  
        axes2[0].grid(True)  

        # Permeate component flows plot
        for j in range(J):  
            col = f'y{j+1}'  
            if col in profile.columns:  
                # Multiply component fraction by normalised permeate flow
                axes2[1].plot(z, profile[col] * profile['cut_p/Qp'], label=f'Component {j+1}')  

        axes2[1].set_xlabel('Normalised Length')  
        axes2[1].set_ylabel('Permeate Normalised Component Flow (-)')  
        axes2[1].set_title(f'{name} Permeate Flow Profile')  
        axes2[1].legend()  
        axes2[1].grid(True)  

        plt.tight_layout()  
        plt.show(block=True)


    def Profile_Export():

        desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop', 'Membrane_Profiles_Cement_Plant_2021.xlsx')

        with pd.ExcelWriter(desktop_path) as writer:
            profile_1.to_excel(writer, sheet_name='Profile_1')
            profile_2.to_excel(writer, sheet_name='Profile_2')

        return





    def Run(Membrane):

        # Set membrane Area based on its feed flow and Q_A_ratio:
        Membrane["Area"] = (Membrane["Feed_Flow"] * 0.0224  * 3600) / Membrane["Q_A_ratio"] # (0.0224 is the molar volume of an ideal gas at STP in m3/mol)

        Membrane["Sweep_Flow"] = 0 # No sweep in this configuration
        Membrane["Sweep_Composition"] = [0] * len(Membrane_1["Permeance"])

        Export_to_mass_balance = Membrane, Component_properties, Fibre_Dimensions

        J = len(Membrane["Permeance"]) #number of components
            
        if Options["Permeance_From_Activation_Energy"]:
            # Obtain Permeance with temperature:
            for i in range(J):
                Membrane["Permeance"][i] = Component_properties["Activation_Energy_Aged"][i][1] * np.exp(-Component_properties["Activation_Energy_Aged"][i][0] / (8.314 * Membrane["Temperature"]))

        results, profile = Hub_Connector(Export_to_mass_balance)
        Membrane["Retentate_Composition"],Membrane["Permeate_Composition"],Membrane["Retentate_Flow"],Membrane["Permeate_Flow"] = results

        #debug
        overall_error = abs(Membrane["Feed_Flow"] + Membrane["Sweep_Flow"] - Membrane["Retentate_Flow"] - Membrane["Permeate_Flow"])/Membrane["Total_Flow"]
        if overall_error > 1e-5:
            print(profile)
            plot_composition_profiles(profile)
            print()
            print(results)
            print()
            print(Membrane)
            raise ValueError(f"Overall mass balance error of membrane {Membrane["Name"]}: Feed + Sweep  - Retentate - Permeate = {overall_error:.2e}")

        #Reformat Permeance and Pressure values to the initial units - will find a smarter way to do this later
        Membrane["Permeance"] = [p / ( 3.348 * 1e-10 ) for p in Membrane["Permeance"]]  # convert from mol/m2.s.Pa to GPU
        Membrane["Pressure_Feed"] *= 1e-5  #convert to bar
        Membrane["Pressure_Permeate"] *= 1e-5  

        if np.any(profile<-1e-5):
            #print(profile)
            #print("Negative values in the membrane profile") #check for negative values in the profile
            raise ConvergenceError  # Return a scalar and an empty dictionary as a placeholder


        #print(profile)
        return results, profile 

    #-----------------------------------------------------------------------------#
    # Run iterations for process recycling loop - Specific to this configuration! #
    #-----------------------------------------------------------------------------#

    max_iter = 150      # maximum number of iterations for the recycling loop
    tolerance = 5e-5    # convergence tolerance for the recycling loop

    Placeholder_1={ #Intermediade data storage for the recycling loop entering the first membrane used to check for convergence
        "Feed_Composition": [0] * J,
        "Feed_Flow": 0,              
        } 

    for j in range(max_iter):

        Mem_Train_Choice(Membrane_1) # Get the correct membrane compression train for Membrane 1
        try:
            results_1 , profile_1 = Run(Membrane_1) # Run the first membrane
        except ConvergenceError:
            raise ValueError ("Convergence error in Membrane 1 ")

        for i in range(J): #results "[0]: x", "[1]: y", "[2]: Q_ret", "[3]: Q_perm"
            Membrane1.set_cell_value(f'D{i+14}', results_1[1][i] * results_1[3] * 3.6) # convert from mol/s to kmol/h and send to unisim
            Membrane1.set_cell_value(f'D{i+21}', results_1[0][i] * results_1[2] * 3.6)

        Mem_Train_Choice(Membrane_2) # Get the correct membrane compression train for Membrane 2
            
        try: 
            results_2, profile_2 = Run(Membrane_2) # Run the second membrane
        except ConvergenceError:
            raise ValueError ("Convergence error in membrane 2")

        for i in range(J): #results "[0]: x", "[1]: y", "[2]: Q_ret", "[3]: Q_perm"
            Membrane2.set_cell_value(f'D{i+14}', results_2[1][i] * results_2[3] * 3.6)
            Membrane2.set_cell_value(f'D{i+21}', results_2[0][i] * results_2[2] * 3.6) 

        Convergence_Composition = sum(abs(np.array(Placeholder_1["Feed_Composition"]) - np.array(Membrane_1["Feed_Composition"])))
        Convergence_Flowrate = abs( ( (Placeholder_1["Feed_Flow"]) - (Membrane_1["Feed_Flow"] ) ) / (Membrane_1["Feed_Flow"] ) / 100 )
        print(f'Convergence Composition: {Convergence_Composition:.3e}, Convergence Flowrate: {Convergence_Flowrate*100:.3e}')
            
        #check for convergence
        if j > 0 and Convergence_Composition < tolerance and Convergence_Flowrate < tolerance:  
            break

        #Ready for next iteration
        Placeholder_1["Feed_Composition"] = Membrane_1["Feed_Composition"]
        Placeholder_1["Feed_Flow"] = Membrane_1["Feed_Flow"]

        
        unisim.wait_solution()


    else: 
        print("Max iterations for recycling reached")
                           

    #-------------------------------------------#
    #----------- Export/Import UniSim ----------#
    #-------------------------------------------#

    def Duty_Gather():  # Gather Duties of the trains from the solved process

        def get_lowest_duty_train(train_data):

            # Filter out trains with None or non-positive compressor duty
            valid_trains = [i for i, train in enumerate(train_data) if train[0] is not None and train[0] > 0 and train[1] is not None and train[1] >0]
            if not valid_trains:
                raise ValueError("No valid trains found with positive compressor duty.")
            
            # Find the train with the lowest compressor duty
            lowest_duty_train_index = min(valid_trains, key=lambda i: train_data[i][0])
            lowest_duty_train = train_data[lowest_duty_train_index]
            lowest_duty_train.append(lowest_duty_train_index) #append the index of the train with the lowest duty to know the equipment count
            return lowest_duty_train

        def gather_train_data(start_row):
            return [
                [
                    Duties.get_cell_value(f'H{i+start_row}'),  # Compressor Duty (kW)
                    Duties.get_cell_value(f'I{i+start_row}'),  # Hex Area (m2)
                    Duties.get_cell_value(f'J{i+start_row}'),  # Water Flowrate (kg/hr)
                    Duties.get_cell_value(f'K{i+start_row}') / 1e6 
                    if Duties.get_cell_value(f'K{i+start_row}') is not None and Duties.get_cell_value(f'K{i+start_row}') > 0 
                    else 0  # Cryogenic Cooler Duty (MJ/hr)
                ]
                for i in range(3)
            ]


        Train1 = gather_train_data(9)
        Train2 = gather_train_data(15)
        Liquefaction = gather_train_data(3)

        # Get the train with the lowest compressor duty for each category
        Train1_lowest = get_lowest_duty_train(Train1)
        Train2_lowest = get_lowest_duty_train(Train2)
        Liquefaction_lowest = get_lowest_duty_train(Liquefaction)

        return Train1_lowest, Train2_lowest, Liquefaction_lowest

    Train1, Train2, Liquefaction = Duty_Gather() # Gather the duties from the solved process

    # Gather the energy recovery form the retentate. Assume flue gas at 1 bar and a maximum temperature of 120 C to match original flue gas.
    Expanders = (Duties.get_cell_value('H21'), Duties.get_cell_value('H24'), Duties.get_cell_value('H27')) # Get the retentate expanders duties (kW)
    Heaters = (Duties.get_cell_value('I21'), Duties.get_cell_value('I24')) # Get the retentate heaters duties (kJ/hr)

    # Gather the cryogenic cooler duties - if any
    Cryogenics = ( (Train1[3], Membrane_1["Temperature"]) , (Train2[3], Membrane_2["Temperature"]) ) # Get the cryogenic cooler duties (MJ/hr) for each membrane train

    # Add information to the compression trains about their number of compressors and coolers
    Train1.append(Train1[4]+1) # Append the number of compressors in the train
    Train1.append(Train1[4]+2)  # Extra heat exchanger for retentate heat recovery

    Train2.append(Train2[4]+1)
    if Process_param["Recycling_Ratio"] == 1: # If the recycling ratio is 1, no extra heat exchanger is needed for retentate heat recovery
        Train2.append(Train2[4]+1)
    else: 
        Train2.append(Train2[4]+2)  # Extra heat exchange for retentate heat recovery

    Liquefaction.append(Liquefaction[4]+3)  # Append the number of compressors and heat exchangers in the liquefaction train
    Liquefaction.append(Liquefaction[4]+3)

    '''
    Process_specs = {
    ...
    "Compressor_trains" : ([duty1, number_of_compressors1], ... , [dutyi, number_of_compressorsi]), # Compressor trains data]
    "Cooler_trains" : ([area1, waterflow1, number_of_coolers1], ... , [areai, waterflowi, number_of_coolersi]), # Cooler trains data
    "Membranes" : (Membrane_1, ..., Membrane_i), # Membrane data)
    "Expanders" : ([expander1_duty], ...[expanderi_duty]), # Expander data
    "Heaters" : ([heater1_duty], ...[heateri_duty]), # Heater data
    "Cryogenics" = ([cooling_power1, temperature1], ... [cooling_poweri, temperaturei]), # Cryogenic cooling data
    }
    '''

    Process_specs = {
        "Feed": Feed,
        "Purity": (results_2[1][0]/(1-results_2[1][-1])),
        "Recovery": results_2[1][0]*results_2[3]/(Feed["Feed_Composition"][0]*Feed["Feed_Flow"]),
        "Compressor_trains": ( (Train1[0], Train1[-2]), (Train2[0],Train2[-2]), (Liquefaction[0], Liquefaction[-2]) ),  # Compressor trains data
        "Cooler_trains": ( (Train1[1], Train2[2], Train1[-1]), (Train2[1],Train2[2], Train2[-1]), (Liquefaction[1], Liquefaction[2], Liquefaction[-1]) ),  # Cooler trains data"
        "Membranes": (Membrane_1, Membrane_2),
        "Expanders": Expanders,  # Expander data
        "Heaters": Heaters,  # Heater data
        "Cryogenics": Cryogenics,
    }

    from Costing import Costing
    Economics = Costing(Process_specs, Process_param, Component_properties)
    print()
    print ("----- Final Results -----")
    print (Economics)

    if Options["Plot_Profiles"]: 
        plot_composition_profiles(profile_1,"Membrane 1")
        plot_composition_profiles(profile_2,"Membrane 2")
    if Options["Export_Profiles"]: Profile_Export()










