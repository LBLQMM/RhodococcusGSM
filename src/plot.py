import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import linregress
from sklearn.metrics import r2_score
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox, TextArea)

# This function takes in a dataframe of flux values, and makes a scatter plot of
#  predicted values vs. measured values
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
    plt.title(r''+ r"$\bf{" + str(method) + "}$"  + ': ' + f"$R^2$={filtered_r2:.2F} ({unfiltered_r2:.2F}$^\star$)", fontsize=24)
    plt.ylabel(f'Predicted Flux', fontsize=22)
    plt.xlabel(f'13C-MFA Flux', fontsize=22)
    plt.legend(fontsize=14)
    
    # add styles to the plot
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    # save and show the plot
    plt.savefig(str(output_dir)+str(substrate)+'_'+str(method)+'_'+str(strain)+'_flux_scatter_plot.png')
    plt.show()
    
# This function takes in a flux data frame and column name, and  
# generates a flux map using the flux values in the specified column
def generate_flux_map(fluxes_df, column_name, substrate='', method='', strain='', output_dir='./'):
    # define plot area
    fig, ax = plt.subplots(figsize=(15, 20), dpi=50)
    
    
    # load unlabeled image and set as plot background
    unlabeled_image = plt.imread('../src/images/unlabeled_flux_map.png')
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
    predicted_rates = list(predicted_growth_parameters_df.loc[conditions, f'{method} growth rate'])
    
    # loop over the data points, and add each point to the plot
    for measured_rate, predicted_rate, condition in zip(measured_rates, predicted_rates, conditions):
        ax.scatter(measured_rate, predicted_rate)
        ax.annotate(condition, (measured_rate, predicted_rate), fontsize=18)
    
    # calculate r-squared value
    r2 = r2_score(measured_rates, predicted_rates)
    # _, _, r2, _, _ = linregress(measured_rates, predicted_rates)
    
    # make a list of absolute errors between measured and predicted rates
    absolute_errors = [
        abs(predicted_rate - measured_rate) 
        for predicted_rate, measured_rate 
        in zip(predicted_rates, measured_rates)
    ]
    
    # make a list of squared errors between measured and predicted rates
    squared_errors = [
        pow(predicted_rate - measured_rate, 2) 
        for predicted_rate, measured_rate 
        in zip(predicted_rates, measured_rates)
    ]
    
    # calculate mean absolute error
    mae = round(sum(absolute_errors) / len(measured_rates), 3)
    
    # calculate sum of squared errors
    ssr = sum(squared_errors) / len(measured_rates)
    print('ssr:', ssr)
    
    # define plot axes limits
    plt.xlim(xlim)
    plt.ylim(ylim)

    # add labels to the plot
    plt.title(r''+ r"$\bf{" + str(method) + "}$" +': '+ f"$R^2$={r2:.2F}, " + f"SSR={ssr:.4F} " , fontsize=24)
    plt.ylabel('Predicted (hr-1)', fontsize=22)
    plt.xlabel('Measured (hr-1)', fontsize=22)
    
    # add styles to the plot
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
    # save and show the plot
    plt.savefig(f'{output_dir}growth_rate_scatter_plot_{method}.png')
    plt.show()