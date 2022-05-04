#from utils import mae_func
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

from utils import *

#define function to compare growth rates in scatter plot
def comparison_scatter_plot(observed, predicted, labels, method, output_dir='./'):
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

    
    plt.xlabel(r'Observed (hr-1)', fontsize=18)#ed growth rate [$mmol/gDW/hr$]',fontsize=14)
    plt.ylabel(r'Predicted (hr-1)', fontsize=18)#ed growth rate [$mmol/gDW/hr$]', fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}, MAE={mae_score}", fontsize=18)#Growth rates: Observed vs. Predicted ('+strtitle+'), \n'
    #plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}, MAE={mae_score}, MSE = {mse}, RMSE={rmse}", fontsize=18)#Growth rates: Observed vs. Predicted ('+strtitle+'), \n'
    plt.savefig(str(output_dir)+'growth_rate_scatter_plot_'+str(method)+'.png')
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
    plt.savefig(str(output_dir)+'Plot_13c_'+str(method)+'_'+str(substrate)+'_'+str(strain)+'.png')
    plt.show()
    
def obs_vs_pred_scatter_plot_with_std(obspred_fluxes, substrate, method, strain, output_dir='./'):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
    pathway_list = sorted(list(set(obspred_fluxes['Pathway'])))
    pathway_list.remove('Biomass Equation') # don't plot this point
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
    


#     if substrate=='phenol':
    plt.title(r''+ r"$\bf{" + str(method) + "}$"  + ': ' + f"$R^2$={r2_scikit_2:.2F} ({r2_scikit_1:.2F}$^\star$), MAE={mae_score_2} ({mae_score_1}$^\star$)", fontsize=18) #star: without 'ATP -> ATP.ext', 'NADH <-> NADPH'
#     else:
#         plt.title(r''+ r"$\bf{" + str(method) + "}$" + ': '+ f"$R^2$={r2_scikit_2:.2F}, MAE={mae_score_2}", fontsize=18)#r''+str(sub)+  ' 13C MFA vs. '+ str(method) + ' Fluxes for ' +linename+ '\n' + f"$R^2$={r2_scikit_2:.2F} (all reactions)", fontsize=18)#, MAE={mae_score}, MSE = {mse}, RMSE={rmse}
#     plt.ylabel(r'Predicted', fontsize=18)#+ str(method) + ' flux (per 100 mmol of '+str(sub)+  ' uptake)')
    plt.ylabel(r'Predicted Flux (per 100 mmol of '+str(sub)+  ' uptake)', fontsize=14)
    plt.xlabel(r'13C MFA Flux (per 100 mmol of '+str(sub)+  ' uptake)', fontsize=14)
    plt.legend(fontsize=14)
    plt.savefig(str(output_dir)+str(substrate)+'_'+str(method)+'_'+str(strain)+'_flux_scatter_plot.png')
#     plt.savefig(str(output_dir)+'Plot_13cwithstd_'+str(method)+'_'+str(substrate)+'_'+str(strain)+'.png')
    plt.show()
    
    
def generate_flux_map(data_df, flux_column, substrate='', method='', strain='', output_dir='./'):
    fig, ax = plt.subplots(figsize=(15, 20), dpi=50)
    xy = (0.5, 0.5)
    arr_img = plt.imread('../src/unlabeled_flux_map.png')
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
    plt.savefig(str(output_dir)+str(substrate)+'_'+str(method)+'_'+str(strain)+'_flux_map.png')
#     plt.savefig(str(output_dir)+'Plot_fluxmap'+str(flux_column)+'s'+'CPM.png')
    plt.show()