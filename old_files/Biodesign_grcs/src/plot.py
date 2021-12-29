from utils import mae_func
import sys
import pandas as pd
import numpy as np
import scipy.stats

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)
from matplotlib.cbook import get_sample_data
import matplotlib.image as mpimg
import matplotlib.cm as cm
from sklearn.metrics import r2_score
from scipy.stats import linregress
import math

#define function to compare growth rates in scatter plot
def scatter_plot_compare(observed, predicted, labels, method, output_dir='./'):
    fig, ax = plt.subplots(figsize=(8, 8))
    
#     lims = [
#                 np.min([observed, predicted]),  # min of both axes
#                 np.max([observed, predicted]),  # max of both axes
#             ]
    lims = [0, 0.4]
    ax.set_xlim(lims)
    
    # Plot Diagonal Dashed Line
    ax.plot(lims, lims, ls="--", color=".8", zorder=0)
    for i in range(0, len(observed)):
        ax.scatter(observed[i], predicted[i])
        ax.annotate(str(labels[i]),(observed[i],predicted[i]), fontsize=14)
        
    #calculate statistical quantities:
    r2 = r2_score(observed,predicted)
    #mse = np.round(mse_func(predicted, observed),2)
    #rmse = np.round(rmse_func(predicted, observed),2)
    mae_score = np.round(mae_func(observed, predicted),2)

    
    plt.xlabel(r'Observations', fontsize=18)#ed growth rate [$mmol/gDW/hr$]',fontsize=14)
    plt.ylabel(r'Predictions', fontsize=18)#ed growth rate [$mmol/gDW/hr$]', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}, MAE={mae_score}", fontsize=18)#Growth rates: Observed vs. Predicted ('+strtitle+'), \n'
    #plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}, MAE={mae_score}, MSE = {mse}, RMSE={rmse}", fontsize=18)#Growth rates: Observed vs. Predicted ('+strtitle+'), \n'
    plt.savefig(str(output_dir)+'Plot_growthrates_'+ str(method) + '.png')
    plt.savefig(str(output_dir)+'Plot_growthrates_'+ str(method) + '.svg')
    plt.show()

# def mae_func(observed, predicted):
#     """Mean Absolute Error.
#     Multioutput case included."""

#     if observed.ndim == 1:
#         return np.mean(np.abs([y_o - y_p for y_o, y_p in zip(observed, predicted)]))
#     else:
#         return [
#             np.mean(
#                 np.abs([y_o - y_p for y_o, y_p in zip(observed[:, i], predicted[:, i])])
#             )
#             for i in range(observed.shape[1])
#         ]

def scatterplotcomp_obs_vs_pred(obspred_fluxes, substrate, method, strain, output_dir='./'):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
    pathway_list = sorted(list(set(obspred_fluxes['Pathway']))) #sorted list such that same colors with every run
    for pathway in pathway_list:
        pathway_df = obspred_fluxes[obspred_fluxes['Pathway'] == pathway]

        measured_flux_list = list(pathway_df['Flux'])
        simulated_flux_list = list(pathway_df[str(method) + ' ' + str(strain) + ' Value'])
        
        if pathway == 'Pentose Phosphate Pathway':
            pathway = 'PP Pathway'

        ax.scatter(measured_flux_list, simulated_flux_list, label=pathway)


    # Dashed line
    x = np.linspace(*ax.get_xlim())
    ax.plot(x, x, ls="--", c=".3")
    
    if substrate=='phenol':
        sub = 'Phenol'
    elif substrate=='glucose':
        sub = 'Glucose'
    else:
        print("Unknown substrate")
    predicted1 = obspred_fluxes.loc[~obspred_fluxes['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']), str(method)+' '+ str(strain) + ' ' + 'Value']
    observed1 =  obspred_fluxes.loc[~obspred_fluxes['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']),'Flux']
    
    predicted2 = obspred_fluxes.loc[:, str(method)+' '+ str(strain) + ' ' + 'Value']
    observed2 =  obspred_fluxes.loc[:,'Flux']
    
    r2_scikit_1 = r2_score(observed1,predicted1)
    r2_scikit_2 = r2_score(observed2,predicted2)
    
    mae_score_1 = np.round(mae_func(observed1, predicted1),2)
    mae_score_2 = np.round(mae_func(observed2, predicted2),2)
    


    if substrate=='phenol':
        plt.title(r''+ r"$\bf{" + str(method) + "}$"  + ': ' + f"$R^2$={r2_scikit_2:.2F} ({r2_scikit_1:.2F}$^\star$), MAE={mae_score_2} ({mae_score_1}$^\star$)", fontsize=18) #star: without 'ATP -> ATP.ext', 'NADH <-> NADPH'
    else:
        plt.title(r''+ r"$\bf{" + str(method) + "}$" + ': '+ f"$R^2$={r2_scikit_2:.2F}, MAE={mae_score_2}", fontsize=18)#r''+str(sub)+  ' 13C MFA vs. '+ str(method) + ' Fluxes for ' +linename+ '\n' + f"$R^2$={r2_scikit_2:.2F} (all reactions)", fontsize=18)#, MAE={mae_score}, MSE = {mse}, RMSE={rmse}
    plt.xlabel(r'Observations', fontsize=18)#13C MFA fluxes (per 100 mmol '+str(sub)+  ' uptake)')
    plt.ylabel(r'Predictions', fontsize=18)#+ str(method) + ' flux (per 100 mmol of '+str(sub)+  ' uptake)')
    plt.legend(fontsize=14)
    plt.savefig(str(output_dir)+'Plot_13c_'+str(method)+'_'+str(substrate)+'_'+str(strain)+'.svg')
    plt.savefig(str(output_dir)+'Plot_13c_'+str(method)+'_'+str(substrate)+'_'+str(strain)+'.png')
    plt.show()
    
def scatterplotcomp_obs_vs_pred_withstd(obspred_fluxes, substrate, method, strain, output_dir='./'):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
    pathway_list = sorted(list(set(obspred_fluxes['Pathway'])))
    for pathway in pathway_list:                          
        pathway_df = obspred_fluxes[obspred_fluxes['Pathway'] == pathway]
        
        measured_flux_list = list(pathway_df['Flux'])
        simulated_flux_list = list(pathway_df[str(method) + ' ' +  str(strain) + ' Value'])
        if str(method) + ' ' +  str(strain) + ' std Value' in pathway_df.columns:
            simulated_std_list = list(pathway_df[str(method) + ' ' +  str(strain) + ' std Value'])
        measured_std_list = list(pathway_df['90% Confidence Upper Bound']-pathway_df['Flux'])
        
        if pathway == 'Pentose Phosphate Pathway':
            pathway = 'PP Pathway'
        sc = ax.scatter(measured_flux_list, simulated_flux_list, label=pathway)
        if str(method) + ' ' +  str(strain) + ' std Value' in pathway_df.columns:
            ax.errorbar(
                measured_flux_list, simulated_flux_list, xerr=[std1 for std1 in measured_std_list], yerr=[1.9*std for std in simulated_std_list],
                    ecolor="gray", ls='none',
                    alpha=0.8)
        else: 
            ax.errorbar(
                measured_flux_list, simulated_flux_list, xerr=[std1 for std1 in measured_std_list],
                    ecolor="gray", ls='none',
                    alpha=0.8)

    # Dashed line
    x = np.linspace(*ax.get_xlim())
    ax.plot(x, x, ls="--", c=".3")
    if substrate=='phenol':
        sub = 'Phenol'
    elif substrate=='glucose':
        sub = 'Glucose'
    else:
        print("Unknown substrate")
    predicted1 = obspred_fluxes.loc[~obspred_fluxes['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']), str(method)+' '+ str(strain) + ' ' + 'Value']
    observed1 =  obspred_fluxes.loc[~obspred_fluxes['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']),'Flux']
    
    predicted2 = obspred_fluxes.loc[:, str(method)+' '+ str(strain) + ' ' + 'Value']
    observed2 =  obspred_fluxes.loc[:,'Flux']
    
    r2_scikit_1 = r2_score(observed1,predicted1)
    r2_scikit_2 = r2_score(observed2,predicted2)
       
    mae_score_1 = np.round(mae_func(observed1, predicted1),2)
    mae_score_2 = np.round(mae_func(observed2, predicted2),2)
    


    if substrate=='phenol':
        plt.title(r''+ r"$\bf{" + str(method) + "}$"  + ': ' + f"$R^2$={r2_scikit_2:.2F} ({r2_scikit_1:.2F}$^\star$), MAE={mae_score_2} ({mae_score_1}$^\star$)", fontsize=18) #star: without 'ATP -> ATP.ext', 'NADH <-> NADPH'
    else:
        plt.title(r''+ r"$\bf{" + str(method) + "}$" + ': '+ f"$R^2$={r2_scikit_2:.2F}, MAE={mae_score_2}", fontsize=18)#r''+str(sub)+  ' 13C MFA vs. '+ str(method) + ' Fluxes for ' +linename+ '\n' + f"$R^2$={r2_scikit_2:.2F} (all reactions)", fontsize=18)#, MAE={mae_score}, MSE = {mse}, RMSE={rmse}
    plt.ylabel(r'Predictions', fontsize=18)#+ str(method) + ' flux (per 100 mmol of '+str(sub)+  ' uptake)')
    plt.legend(fontsize=14)
    plt.savefig(str(output_dir)+'Plot_13cwithstd_'+str(method)+'_'+str(substrate)+'_'+str(strain)+'.svg')
    plt.savefig(str(output_dir)+'Plot_13cwithstd_'+str(method)+'_'+str(substrate)+'_'+str(strain)+'.png')
    plt.show()
    
    
def map_flux_results(data_df, flux_column, output_dir='./'):
    fig, ax = plt.subplots(figsize=(15, 20), dpi=50)
    xy = (0.5, 0.5)
    arr_img = plt.imread('./unlabeled_flux_map.png')
    imagebox = OffsetImage(arr_img)
    imagebox.image.axes = ax
    ab = AnnotationBbox(imagebox, xy, frameon=False)
    ax.add_artist(ab)

    for _, row in data_df.iterrows():
        if not pd.isnull(row['Location on map']):
            location =  row['Location on map'].replace('(', '').replace(')', '')
            location_list = location.split(',')
            location_tuple = tuple((int(location_list[0]), int(location_list[1])))

            offsetbox = TextArea(f'{row[flux_column]:.1f}',textprops=dict(fontsize=22))
            ab = AnnotationBbox(offsetbox, xy,
                                xybox=location_tuple,
                                xycoords='data',
                                boxcoords="offset points",
                                frameon=False)
            ax.add_artist(ab)

    # Fix the display limits to see everything
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    plt.savefig(str(output_dir)+'Plot_fluxmap'+str(flux_column)+'s'+'CPM.svg')
    plt.savefig(str(output_dir)+'Plot_fluxmap'+str(flux_column)+'s'+'CPM.png')
    plt.show()
    
def stats_for_trial(growth_data, substrate_data, molar_mass, display=False, max_time=0, substrate=''):
    
    biomass_values = growth_data['Biomass Conc']
    biomass_times = growth_data['Hours']
    biomass_init = list(biomass_values)[0]

    substrate_values = substrate_data['Value']*1000/molar_mass
    substrate_times = substrate_data['Hours']
    substrate_init = list(substrate_values)[0]
    
    # growth is the slope of log(biomass) vs. time
    growth_rate, _, _, _, _ = linregress(biomass_times, [math.log(val) for val in biomass_values])
    
    # biomass X = X0*e^(μ*t)
    # This is different from above to ensure that there is a biomass value for every substrate measurement
    biomass_sim = [biomass_init*math.exp(growth_rate*time) for time in substrate_times]
    
    # actual consumption = S0 - S
    sub_consumed = [substrate_init - sub_value for sub_value in substrate_values]
    
    # new biomass X = X0 - X
    biomass_sim_growth = [sim_value - biomass_init for sim_value in biomass_sim ]
    
    # yield is the amount of biomass that can be made from a mmol of substrate
    yield_coeff, _, _, _, _ = linregress(sub_consumed, biomass_sim_growth)

    # S = S0 - (1/yield)*X
    substrate_sim = [substrate_init - 1/yield_coeff*val for val in biomass_sim_growth]
    
    substrate_consumption_rate = (1/yield_coeff) * growth_rate

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 5))
    axes[0].plot(biomass_times, biomass_values, 'o', color='black')
    axes[0].plot(substrate_times, biomass_sim, '-', color='black')
    axes[1].plot(substrate_times, substrate_values, 'o', color='blue')
    axes[1].plot(substrate_times, substrate_sim, '-', color='blue')
    axes[0].set_title('Biomass growth')
    axes[1].set_title(f'{substrate} consumption')
    axes[0].set_xlabel('Time (hr)')
    axes[1].set_xlabel('Time (hr)')
    axes[0].set_ylabel('Biomass (g/L)')
    axes[1].set_ylabel(f'{substrate} (mmol/L)')
    fig.tight_layout()
    
    if display:
        print(f'growth_rate = {growth_rate:.3f} hr-1')
        print(f'yield coefficient = {yield_coeff:.3f} g biomass / mmol substrate')
        print(f'substrate consumption rate = {substrate_consumption_rate:.3f} mmol substrate/gram biomass * hr')
        return growth_rate, yield_coeff, substrate_consumption_rate
    else:
        return growth_rate, yield_coeff, substrate_consumption_rate
    
    
    
def stats_for_condtion(od_df, sub_df, trial_1, trial_2, trial_3, molar_mass, substrate='', max_time=0):
    
    if max_time != 0:
        od_df = od_df[od_df['Hours'] < max_time]
        sub_df = sub_df[sub_df['Hours'] < max_time]
        
    od_1 = od_df[od_df['Line Name'] == trial_1]
    sub_1 = sub_df[sub_df['Line Name'] == trial_1]

    od_2 = od_df[od_df['Line Name'] == trial_2]
    sub_2 = sub_df[sub_df['Line Name'] == trial_2]

    od_3 = od_df[od_df['Line Name'] == trial_3]
    sub_3 = sub_df[sub_df['Line Name'] == trial_3]

    gr_1, yc_1, scr_1 = stats_for_trial(od_1, sub_1, molar_mass, substrate=substrate)
    gr_2, yc_2, scr_2 = stats_for_trial(od_2, sub_2, molar_mass, substrate=substrate)
    gr_3, yc_3, scr_3 = stats_for_trial(od_3, sub_3, molar_mass, substrate=substrate)
    
    growth_rate = np.average([gr_1, gr_2, gr_3])
    yield_coeff = np.average([yc_1, yc_2, yc_3])
    substrate_consumption_rate = np.average([scr_1, scr_2, scr_3])
    growth_rate_std = np.std([gr_1, gr_2, gr_3])
    yield_coeff_std = np.std([yc_1, yc_2, yc_3])
    substrate_consumption_rate_std = np.std([scr_1, scr_2, scr_3])
    
    print(f'growth_rate = {growth_rate:.3f} ± {growth_rate_std:.3f} hr-1')
    print(f'yield coefficient = {yield_coeff:.3f} ± {yield_coeff_std:.3f} g biomass / mmol substrate')
    print(f'substrate consumption rate = {substrate_consumption_rate:.3f} ± {substrate_consumption_rate_std:.3f} mmol substrate/gram biomass * hr')
    return growth_rate, yield_coeff, substrate_consumption_rate, growth_rate_std, yield_coeff_std, substrate_consumption_rate_std 