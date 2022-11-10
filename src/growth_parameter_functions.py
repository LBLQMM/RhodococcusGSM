import numpy as np
from scipy import stats
import math
import matplotlib.pyplot as plt

# These functions are used in notebook D

# This function takes in time course OD (optical density) and substrate concentration, and returns the average
# growth parameters (growth rate, yield coefficient, and substrate uptake rate)
def get_average_growth_parameters(od_df, sub_df, trial_1, trial_2, trial_3, molar_mass, substrate='', max_time=0):
    
    # remove data from after the specified maximum time
    if max_time != 0:
        od_df = od_df[od_df['Hours'] < max_time]
        sub_df = sub_df[sub_df['Hours'] < max_time]
        
    # isolate the data from each trial using the dataframes and trial names
    od_1 = od_df[od_df['Line Name'] == trial_1]
    sub_1 = sub_df[sub_df['Line Name'] == trial_1]

    od_2 = od_df[od_df['Line Name'] == trial_2]
    sub_2 = sub_df[sub_df['Line Name'] == trial_2]

    od_3 = od_df[od_df['Line Name'] == trial_3]
    sub_3 = sub_df[sub_df['Line Name'] == trial_3]

    # get the growth parameters for each trial
    gr_1, yc_1, scr_1 = get_trial_growth_parameters(od_1, sub_1, molar_mass, substrate=substrate)
    gr_2, yc_2, scr_2 = get_trial_growth_parameters(od_2, sub_2, molar_mass, substrate=substrate)
    gr_3, yc_3, scr_3 = get_trial_growth_parameters(od_3, sub_3, molar_mass, substrate=substrate)
    
    # calculate the average parameter values
    growth_rate = np.average([gr_1, gr_2, gr_3])
    yield_coefficient = np.average([yc_1, yc_2, yc_3])
    substrate_uptake_rate = np.average([scr_1, scr_2, scr_3])
    
    # calculate the standard deviation of parameter values
    growth_rate_std = np.std([gr_1, gr_2, gr_3])
    yield_coefficient_std = np.std([yc_1, yc_2, yc_3])
    substrate_uptake_rate_std = np.std([scr_1, scr_2, scr_3])
    
    # return the average and standard deviations of the 
    return growth_rate, yield_coefficient, substrate_uptake_rate, growth_rate_std, yield_coefficient_std, substrate_uptake_rate_std 

    
# This function takes in time course growth and consumption data 
# and calculates the growth parameters from a single trial
def get_trial_growth_parameters(growth_data, substrate_data, molar_mass, max_time=0, substrate=''):
    
    # get the biomass concentrations and time points from the biomass dataframe
    biomass_values = growth_data['Biomass Conc']
    biomass_times = growth_data['Hours']
    
    # get the substrate concentrations and time points from the substrate dataframe
    substrate_values = substrate_data['Value']*1000/molar_mass
    substrate_times = substrate_data['Hours']
    
    # get the initial substrate and biomass concentrations
    starting_substrate = list(substrate_values)[0]
    starting_biomass = list(biomass_values)[0]
    
    # caluculate the growth rate
    # growth rate is equal to the slope of log(biomass) vs. time
    growth_rate, _, _, _, _ = stats.linregress(biomass_times, [math.log(val) for val in biomass_values])
    
    # determine the fitted biomass concentration at each time point in the substrate concentration data set
    # This is needed to ensure that there is a biomass value for every substrate measurement
    # biomass X = X0*e^(μ*t)
    fitted_biomass_conc = [starting_biomass * math.exp(growth_rate * time) for time in substrate_times]
    
    # calculate the amount of substrate consumed at each time point
    substrate_consumed = [starting_substrate - substrate_conc for substrate_conc in substrate_values]
    
    # calculate the amount of biomass produced at each time point
    fitted_biomass_produced = [biomass_conc - starting_biomass for biomass_conc in fitted_biomass_conc]
    
    # calculate the yield coefficient
    # the yield coefficient is equal to the slope of biomass produced vs. substrate consumed
    yield_coefficient, _, _, _, _ = stats.linregress(substrate_consumed, fitted_biomass_produced)
    
    # calculate the substrate uptake rate
    # the substrate uptake rate is equal to the inverse of the yield coefficient times the growth rate
    substrate_uptake_rate = (1/yield_coefficient) * growth_rate

    # calculate the fitted amount of substrate consumed at each time point using yield coeefficient    
    fitted_substrate_conc = [starting_substrate - (1 / yield_coefficient) * biomass_produced 
                                 for biomass_produced in fitted_biomass_produced]

    # define a plotting area with two subplots
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 5))
    
    # plot biomass data on the left
    axes[0].set_title('Biomass Growth', fontsize=16)
    axes[0].set_xlabel('Time (hr)', fontsize=14)
    axes[0].set_ylabel('Biomass (g/L)', fontsize=14)
    
    # plot experimental biomass concentration data points
    axes[0].plot(biomass_times, biomass_values, 'o', color='black')
    # plot fitted biomass concentration curve
    axes[0].plot(substrate_times, fitted_biomass_conc, '-', color='black')
    
    # plot substrate consumption data on the right
    axes[1].set_title(f'{substrate.capitalize()} Consumption', fontsize=16)
    axes[1].set_ylabel(f'{substrate.capitalize()} (mmol/L)', fontsize=14)
    axes[1].set_xlabel('Time (hr)', fontsize=14)

    # plot experimental substrate concentration data points
    axes[1].plot(substrate_times, substrate_values, 'o', color='blue')
    # plot fitted substrate concentration curve
    axes[1].plot(substrate_times, fitted_substrate_conc, '-', color='blue')

    return growth_rate, yield_coefficient, substrate_uptake_rate