import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import linregress
from sklearn.metrics import r2_score
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox, TextArea)

def flux_prediction_scatterplot(fluxes_df, substrate, method, strain, output_dir='./'):
    # define column names
    prediction_column_name = f'{method} {strain} Flux'
    prediction_std_column_name = f'{method} {strain} Flux Std'
    
    # define plot area
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # create a list of pathways from the dataframe, this is sorted so the colors match on all plots
    pathway_list = sorted(list(set(fluxes_df['Pathway'])))
    pathway_list.remove('Biomass Equation') # don't plot this point (I could remove this later)
    
    # loop over each pathway
    for pathway in pathway_list:
        # make a pathway specific filtered dataframe
        pathway_df = fluxes_df[fluxes_df['Pathway'] == pathway]
        
        # make lists of measured and predicted fluxes
        measured_flux_list = list(pathway_df['13C Flux'])
        predicted_flux_list = list(pathway_df[prediction_column_name])
                    
        # shorten the label for pentose phosphate pathway (could update the original csv's)
        if pathway == 'Pentose Phosphate Pathway':
            pathway = 'PP Pathway'
            
        # add the points for each pathway to the plot
        sc = ax.scatter(measured_flux_list, predicted_flux_list, label=pathway)
        
        # this is incorrect (update later with custom upper and lower bound errors)
        # calculate horizontal error bars
        measured_std_list = list(pathway_df['13C Upper Bound']-pathway_df['13C Flux'])
        
        # if prediction standard deviations are included, then add vertical and horizontal error bars
        if prediction_std_column_name in pathway_df.columns:
            # make a list of predicted flux standard deviations
            predicted_std_list = list(pathway_df[prediction_std_column_name])
            
            ax.errorbar(
                measured_flux_list, 
                predicted_flux_list, 
                xerr=[std1 for std1 in measured_std_list], 
                yerr=[1.9*std for std in predicted_std_list], 
                ecolor="gray", 
                ls='none', 
                alpha=0.8
            )
            
        # if prediction standard deviations are NOT included, then only add horizontal error bars
        else: 
            ax.errorbar(
                measured_flux_list, 
                predicted_flux_list, 
                xerr=[std1 for std1 in measured_std_list], 
                ecolor="gray", 
                ls='none', 
                alpha=0.8
            )

    # add 45 degree dashed line
    x = np.linspace(*ax.get_xlim())
    ax.plot(x, x, ls="--", c=".3")
    
    # make lists of the measured fluxes and the predicted fluxes for statistic calculations
    unfiltered_measurements =  fluxes_df.loc[:, '13C Flux']
    unfiltered_predictions = fluxes_df.loc[:, prediction_column_name]
    
    # make list with outliers filtered out
    filtered_measurements =  fluxes_df.loc[~fluxes_df['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']), '13C Flux']
    filtered_predictions = fluxes_df.loc[~fluxes_df['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']), prediction_column_name]
    
    # calculate r-squared values using scikitLearn
    # unfiltered_r2 = r2_score(unfiltered_measurements, unfiltered_predictions)
    # filtered_r2 = r2_score(filtered_measurements, filtered_predictions)
    
    # calculate r-squared values using 
    _, _, unfiltered_r, _, _ = linregress(unfiltered_measurements, unfiltered_predictions)
    unfiltered_r2 = unfiltered_r * unfiltered_r
    
    _, _, filtered_r, _, _ = linregress(filtered_measurements, filtered_predictions)
    filtered_r2 = filtered_r * filtered_r
    
    # add labels to the plot
    plt.title(r''+ r"$\bf{" + str(method) + "}$"  + ': ' + f"$R^2$={filtered_r2:.2F} ({unfiltered_r2:.2F}$^\star$)", fontsize=18)
    plt.ylabel(f'Predicted Flux (per 100 mmol of {substrate.capitalize()} Uptake)', fontsize=14)
    plt.xlabel(f'13C MFA Flux (per 100 mmol of {substrate.capitalize()} Uptake)', fontsize=14)
    plt.legend(fontsize=14)
    
    # save and show the plot
    plt.savefig(str(output_dir)+str(substrate)+'_'+str(method)+'_'+str(strain)+'_flux_scatter_plot.png')
    plt.show()
    
# This function takes in a flux data frame and column name, and  
# generates a flux map using the flux values in the specified column
def generate_flux_map(fluxes_df, column_name, substrate='', method='', strain='', output_dir='./'):
    # define plot area
    fig, ax = plt.subplots(figsize=(15, 20), dpi=50)
    
    
    # load unlabeled image and set as plot background
    unlabeled_image = plt.imread('../src/unlabeled_flux_map.png')
    imagebox = OffsetImage(unlabeled_image)
    imagebox.image.axes = ax
    xy = (0.5, 0.5)
    ab = AnnotationBbox(imagebox, xy, frameon=False)
    ax.add_artist(ab)

    # loop over each reaction in the dataframe
    for _, row in fluxes_df.iterrows():
        # check that there is a location for the reaction's flux
        if not pd.isnull(row['Location on map']):
            # convert the location string to an integer tuple
            location =  row['Location on map'].replace('(', '').replace(')', '')
            location_list = location.split(',')
            location_tuple = tuple((int(location_list[0]), int(location_list[1])))

            # put flux value in a text area, and put the textarea on the plot area
            offsetbox = TextArea(f'{row[column_name]:.1f}',textprops=dict(fontsize=22))
            ab = AnnotationBbox(offsetbox, xy,
                                xybox=location_tuple,
                                xycoords='data',
                                boxcoords="offset points",
                                frameon=False)
            ax.add_artist(ab)

    # ensure that the axes have minimal styles 
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis("off")
    
    # save and show the plot
    plt.savefig(f'{output_dir}{substrate}_{method}_{strain}_flux_map.png')
    plt.show()
    
# This function takes in the predicted_growth_parameters data frame, the method name, and x and y boundaries,
# and it displays a scatterplot of predicted and measured growth rates
def growth_rate_scatterplot(predicted_growth_parameters_df, method, xlim, ylim, output_dir='./'):
    # define plot area
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # plot diagonal dashed line
    ax.plot(xlim, ylim, ls="--", color=".8", zorder=0)
    
    # get the names of the conditions
    conditions = list(predicted_growth_parameters_df.index)
    
    # get measured and predicted rates from the dataframe
    measured_rates = list(predicted_growth_parameters_df.loc[conditions, 'growth rate'])
    predicted_rates = list(predicted_growth_parameters_df.loc[conditions, 'FBA growth rate'])
    
    # loop over the data points, and add each point to the plot
    for measured_rate, predicted_rate, condition in zip(measured_rates, predicted_rates, conditions):
        ax.scatter(measured_rate, predicted_rate)
        ax.annotate(condition, (measured_rate, predicted_rate), fontsize=14)
    
    # calculate r-squared value
    r2 = r2_score(measured_rates, predicted_rates)
    # _, _, r2, _, _ = linregress(measured_rates, predicted_rates)
    
    # calculate the mean absolute error between measured and predicted rates
    errors = [abs(measured_rate - predicted_rates) for measured_rates, predicted_rates in zip(measured_rates, predicted_rates)]
    mae = round(sum(errors) / len(measured_rates), 3)
    
    # define plot axes limits
    plt.xlim(xlim)
    plt.ylim(ylim)

    # add labels to the plot
    plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}, " + f"MAE={mae:.2F} " , fontsize=18)
    plt.ylabel('Predicted (hr-1)', fontsize=18)
    plt.xlabel('Measured (hr-1)', fontsize=18)
    
    # add styles to the plot
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    # save and show the plot
    plt.savefig(f'{output_dir} growth_rate_scatter_plot_{method}.png')
    plt.show()
    
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------    
# ----------------------------------------------------------------------------------------------
    
# def growth_rate_scatterplot(observed, predicted, labels, method, xlim, ylim, output_dir='./'):
#     fig, ax = plt.subplots(figsize=(8, 8))
    
#     # Plot Diagonal Dashed Line
#     ax.plot(xlim, ylim, ls="--", color=".8", zorder=0)
#     for i in range(0, len(observed)):
#         ax.scatter(observed[i], predicted[i])
#         ax.annotate(str(labels[i]),(observed[i],predicted[i]), fontsize=14)
        
#     #calculate statistical quantities:
# #     slope, intercept, r, p, se = linregress(observed, predicted)
# #     r2 = r*r
#     mae = np.round(mae_func(observed, predicted),2)
#     r2 = r2_score(observed,predicted)

    
#     plt.xlabel(r'Observed (hr-1)', fontsize=18)#ed growth rate [$mmol/gDW/hr$]',fontsize=14)
#     plt.ylabel(r'Predicted (hr-1)', fontsize=18)#ed growth rate [$mmol/gDW/hr$]', fontsize=14)
#     plt.xticks(fontsize=14)
#     plt.yticks(fontsize=14)
#     plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}", fontsize=18)
#     # plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}, " + f"MAE={mae:.2F} " , fontsize=18)
#     plt.xlim(xlim)
#     plt.ylim(ylim)
#     plt.savefig(str(output_dir)+'growth_rate_scatter_plot_'+str(method)+'.png')
#     plt.show()

# def generate_flux_map(data_df, flux_column, substrate='', method='', strain='', output_dir='./'):
#     fig, ax = plt.subplots(figsize=(15, 20), dpi=50)
#     xy = (0.5, 0.5)
#     arr_img = plt.imread('../src/unlabeled_flux_map.png')
#     imagebox = OffsetImage(arr_img)
#     imagebox.image.axes = ax
#     ab = AnnotationBbox(imagebox, xy, frameon=False)
#     ax.add_artist(ab)

#     for _, row in data_df.iterrows():
#         if not pd.isnull(row['Location on map']):
#             location =  row['Location on map'].replace('(', '').replace(')', '')
#             location_list = location.split(',')
#             location_tuple = tuple((int(location_list[0]), int(location_list[1])))

#             offsetbox = TextArea(f'{row[flux_column]:.1f}',textprops=dict(fontsize=22))
#             ab = AnnotationBbox(offsetbox, xy,
#                                 xybox=location_tuple,
#                                 xycoords='data',
#                                 boxcoords="offset points",
#                                 frameon=False)
#             ax.add_artist(ab)

#     # Fix the display limits to see everything
#     ax.set_xlim(0, 1)
#     ax.set_ylim(0, 1)
#     ax.set_xticks([])
#     ax.set_yticks([])
#     ax.axis("off")
#     plt.savefig(str(output_dir)+str(substrate)+'_'+str(method)+'_'+str(strain)+'_flux_map.png')
#     plt.show()
    
    
#define function to compare growth rates in scatter plot
# def comparison_scatter_plot(observed, predicted, labels, method, xlim, ylim, output_dir='./'):
#     fig, ax = plt.subplots(figsize=(8, 8))
    
#     # Plot Diagonal Dashed Line
#     ax.plot(xlim, ylim, ls="--", color=".8", zorder=0)
#     for i in range(0, len(observed)):
#         ax.scatter(observed[i], predicted[i])
#         ax.annotate(str(labels[i]),(observed[i],predicted[i]), fontsize=14)
        
#     #calculate statistical quantities:
# #     slope, intercept, r, p, se = linregress(observed, predicted)
# #     r2 = r*r
#     mae = np.round(mae_func(observed, predicted),2)
#     r2 = r2_score(observed,predicted)

    
#     plt.xlabel(r'Observed (hr-1)', fontsize=18)#ed growth rate [$mmol/gDW/hr$]',fontsize=14)
#     plt.ylabel(r'Predicted (hr-1)', fontsize=18)#ed growth rate [$mmol/gDW/hr$]', fontsize=14)
#     plt.xticks(fontsize=14)
#     plt.yticks(fontsize=14)
#     plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}", fontsize=18)
#     # plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}, " + f"MAE={mae:.2F} " , fontsize=18)
#     plt.xlim(xlim)
#     plt.ylim(ylim)
#     plt.savefig(str(output_dir)+'growth_rate_scatter_plot_'+str(method)+'.png')
#     plt.show()

# def scatterplotcomp_obs_vs_pred(obspred_fluxes, substrate, method, strain, output_dir='./'):
#     fig, ax = plt.subplots(figsize=(8, 8))
#     ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
#     pathway_list = sorted(list(set(obspred_fluxes['Pathway']))) #sorted list such that same colors with every run
#     for pathway in pathway_list:
#         pathway_df = obspred_fluxes[obspred_fluxes['Pathway'] == pathway]

#         measured_flux_list = list(pathway_df['13C Flux'])
#         simulated_flux_list = list(pathway_df[str(method) + ' ' + str(strain) + ' Flux'])
        
#         if pathway == 'Pentose Phosphate Pathway':
#             pathway = 'PP Pathway'

#         ax.scatter(measured_flux_list, simulated_flux_list, label=pathway)


#     # Dashed line
#     x = np.linspace(*ax.get_xlim())
#     ax.plot(x, x, ls="--", c=".3")
    
#     if substrate=='phenol':
#         sub = 'Phenol'
#     elif substrate=='glucose':
#         sub = 'Glucose'
#     else:
#         print("Unknown substrate")
#     predicted1 = obspred_fluxes.loc[~obspred_fluxes['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']), str(method)+' '+ str(strain) + ' Flux']
#     observed1 =  obspred_fluxes.loc[~obspred_fluxes['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']),'13C Flux']
    
#     predicted2 = obspred_fluxes.loc[:, str(method)+' '+ str(strain) + ' Flux']
#     observed2 =  obspred_fluxes.loc[:,'13C Flux']
    
#     slope, intercept, r_filtered, p, se = linregress(observed1, predicted1)
#     r2_scikit_1 = r_filtered*r_filtered
    
#     slope, intercept, r_unfiltered, p, se = linregress(observed2, predicted2)
#     r2_scikit_2 = r_unfiltered*r_unfiltered
    
    
# #     r2_scikit_1 = r2_score(observed1,predicted1)
# #     r2_scikit_2 = r2_score(observed2,predicted2)
    
#     mae_score_1 = np.round(mae_func(observed1, predicted1),2)
#     mae_score_2 = np.round(mae_func(observed2, predicted2),2)
    


#     if substrate=='phenol':
#         plt.title(r''+ r"$\bf{" + str(method) + "}$"  + ': ' + f"$R^2$={r2_scikit_2:.2F} ({r2_scikit_1:.2F}$^\star$), MAE={mae_score_2} ({mae_score_1}$^\star$)", fontsize=18) #star: without 'ATP -> ATP.ext', 'NADH <-> NADPH'
#     else:
#         plt.title(r''+ r"$\bf{" + str(method) + "}$" + ': '+ f"$R^2$={r2_scikit_2:.2F}, MAE={mae_score_2}", fontsize=18)#r''+str(sub)+  ' 13C MFA vs. '+ str(method) + ' Fluxes for ' +linename+ '\n' + f"$R^2$={r2_scikit_2:.2F} (all reactions)", fontsize=18)#, MAE={mae_score}, MSE = {mse}, RMSE={rmse}
#     plt.xlabel(r'Observations', fontsize=18)#13C MFA fluxes (per 100 mmol '+str(sub)+  ' uptake)')
#     plt.ylabel(r'Predictions', fontsize=18)#+ str(method) + ' flux (per 100 mmol of '+str(sub)+  ' uptake)')
#     plt.legend(fontsize=14)
#     plt.savefig(str(output_dir)+'Plot_13c_'+str(method)+'_'+str(substrate)+'_'+str(strain)+'.png')
#     plt.show()
    

    
def mae_func(observed, predicted):
    """Mean Absolute Error.
    Multioutput case included."""

    if observed.ndim == 1:
        return np.mean(np.abs([y_o - y_p for y_o, y_p in zip(observed, predicted)]))
    else:
        return [
            np.mean(
                np.abs([y_o - y_p for y_o, y_p in zip(observed[:, i], predicted[:, i])])
            )
            for i in range(observed.shape[1])
        ]

# save old version for now
# def obs_vs_pred_scatter_plot_with_std(obspred_fluxes, substrate, method, strain, output_dir='./'):
#     fig, ax = plt.subplots(figsize=(8, 8))
#     ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")
#     pathway_list = sorted(list(set(obspred_fluxes['Pathway'])))
#     pathway_list.remove('Biomass Equation') # don't plot this point
#     for pathway in pathway_list:                          
#         pathway_df = obspred_fluxes[obspred_fluxes['Pathway'] == pathway]
        
#         measured_flux_list = list(pathway_df['13C Flux'])
#         simulated_flux_list = list(pathway_df[str(method) + ' ' +  str(strain) + ' Flux'])
#         if str(method) + ' ' +  str(strain) + ' Flux Std' in pathway_df.columns:
#             simulated_std_list = list(pathway_df[str(method) + ' ' +  str(strain) + ' Flux Std'])
#         measured_std_list = list(pathway_df['13C Upper Bound']-pathway_df['13C Flux'])
        
#         if pathway == 'Pentose Phosphate Pathway':
#             pathway = 'PP Pathway'
#         sc = ax.scatter(measured_flux_list, simulated_flux_list, label=pathway)
#         if str(method) + ' ' +  str(strain) + ' Flux Std' in pathway_df.columns:
#             ax.errorbar(
#                 measured_flux_list, simulated_flux_list, xerr=[std1 for std1 in measured_std_list], yerr=[1.9*std for std in simulated_std_list],
#                     ecolor="gray", ls='none',
#                     alpha=0.8)
#         else: 
#             ax.errorbar(
#                 measured_flux_list, simulated_flux_list, xerr=[std1 for std1 in measured_std_list],
#                     ecolor="gray", ls='none',
#                     alpha=0.8)

#     # Dashed line
#     x = np.linspace(*ax.get_xlim())
#     ax.plot(x, x, ls="--", c=".3")
#     if substrate=='phenol':
#         sub = 'Phenol'
#     elif substrate=='glucose':
#         sub = 'Glucose'
#     else:
#         print("Unknown substrate")
#     predicted1 = obspred_fluxes.loc[~obspred_fluxes['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']), str(method)+' '+ str(strain) + ' Flux']
#     observed1 =  obspred_fluxes.loc[~obspred_fluxes['Reaction'].isin(['ATP -> ATP.ext', 'NADH <-> NADPH']),'13C Flux']
    
#     predicted2 = obspred_fluxes.loc[:, str(method)+' '+ str(strain) + ' Flux']
#     observed2 =  obspred_fluxes.loc[:,'13C Flux']
    
#     slope, intercept, r_filtered, p, se = linregress(observed1, predicted1)
#     r2_scikit_1 = r_filtered*r_filtered
    
#     slope, intercept, r_unfiltered, p, se = linregress(observed2, predicted2)
#     r2_scikit_2 = r_unfiltered*r_unfiltered
    
# #     r2_scikit_1 = r2_score(observed1,predicted1)
# #     r2_scikit_2 = r2_score(observed2,predicted2)
       
#     mae_score_1 = np.round(mae_func(observed1, predicted1),2)
#     mae_score_2 = np.round(mae_func(observed2, predicted2),2)
    


# #     if substrate=='phenol':
# #     plt.title(r''+ r"$\bf{" + str(method) + "}$"  + ': ' + f"$R^2$={r2_scikit_2:.2F} ({r2_scikit_1:.2F}$^\star$), MAE={mae_score_2} ({mae_score_1}$^\star$)", fontsize=18) 
    
#     plt.title(r''+ r"$\bf{" + str(method) + "}$"  + ': ' + f"$R^2$={r2_scikit_2:.2F} ({r2_scikit_1:.2F}$^\star$)", fontsize=18)
    
#     #star: without 'ATP -> ATP.ext', 'NADH <-> NADPH'
# #     else:
# #         plt.title(r''+ r"$\bf{" + str(method) + "}$" + ': '+ f"$R^2$={r2_scikit_2:.2F}, MAE={mae_score_2}", fontsize=18)#r''+str(sub)+  ' 13C MFA vs. '+ str(method) + ' Fluxes for ' +linename+ '\n' + f"$R^2$={r2_scikit_2:.2F} (all reactions)", fontsize=18)#, MAE={mae_score}, MSE = {mse}, RMSE={rmse}
# #     plt.ylabel(r'Predicted', fontsize=18)#+ str(method) + ' flux (per 100 mmol of '+str(sub)+  ' uptake)')
#     plt.ylabel(r'Predicted Flux (per 100 mmol of '+str(sub)+  ' uptake)', fontsize=14)
#     plt.xlabel(r'13C MFA Flux (per 100 mmol of '+str(sub)+  ' uptake)', fontsize=14)
#     plt.legend(fontsize=14)
#     plt.savefig(str(output_dir)+str(substrate)+'_'+str(method)+'_'+str(strain)+'_flux_scatter_plot.png')
# #     plt.savefig(str(output_dir)+'Plot_13cwithstd_'+str(method)+'_'+str(substrate)+'_'+str(strain)+'.png')
#     plt.show()