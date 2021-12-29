#Import python packages:
import pandas as pd
import numpy as np
from ensemblemethods import SPOT
import cobra

#################################################################
########Utilities needed for prediction methods and data formatting##############
#################################################################
"""
    Transforms data to data frame with just index as gene identifiers and one column for values, this data frame can then be used for EFLUX2 and SPOT prediction functions.
    
        Parameters
        ----------
        transdata: data frame for data imported from EDD.
        linename: Line Name in EDD
        
        Returns
        -------
        transdataWTPR1new2: reformatted data frame
        
"""

#!!!!TODO: Generalize for multiple time points
def construct_trans_df(transdata, linename):
    transdataWTPR1 = transdata[transdata['Line Name']==linename]
    transdataWTPR1new = transdataWTPR1.filter(['Value', 'Measurement Type'])
    transdataWTPR1new2 = transdataWTPR1new.set_index('Measurement Type')
    return transdataWTPR1new2

"""
    Make FBA predictons.
    
        Parameters
        ----------
        model: cobrapy model.
        substrate: carbon source used.
        sub_uptake_rate: uptake rate for carbon source used.
        
        Returns
        -------
        fbasol: FBA prediction solution.
        
"""

def FBA_pred(model, substrate, sub_uptake_rate=100):
    with model:
        medium = model.medium 
        if substrate=='phenol':
            model.objective = 'Growth_Phenol'
            medium = {key:1000 for (key,value) in model.medium.items()}
            model.reactions.get_by_id('Growth_Glucose').upper_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            
            #remove all non-phenol carbon sources:
            medium["EX_glc__D_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            medium['EX_tag'] = 0
            medium["EX_phenol_e"] = sub_uptake_rate
            
            model.reactions.get_by_id('EX_glc__D_e').upper_bound = 0
            model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
            
        elif substrate=='glucose':
            model.objective = 'Growth_Glucose'
            medium = {key:1000 for (key,value) in model.medium.items()}
            model.reactions.get_by_id('Growth_Phenol').upper_bound = 0
            model.reactions.get_by_id('Growth_Phenol').lower_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            
            #remove all non-glucose carbon sources:
            medium["EX_phenol_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            medium['EX_tag'] = 0
            medium["EX_glc__D_e"] = sub_uptake_rate
            
            model.reactions.get_by_id('EX_phenol_e').upper_bound = 0
            model.reactions.get_by_id('EX_phenol_e').lower_bound = 0
        else:
            print('Unknown substrate: Please choose among phenol and glucose')
        model.medium = medium
        fbasol = model.optimize()
        display(model.medium)
    return fbasol

"""
    Make FBA predictons.
    
        Parameters
        ----------
        model: cobrapy model.
        substrate: carbon source used.
        sub_uptake_rate: uptake rate for carbon source used.
        
        Returns
        -------
        pFBA_solution_all: pFBA prediction solution.
        
"""

def pFBA_pred(model, substrate, sub_uptake_rate=100):
    with model:
        medium = model.medium 
        if substrate=='phenol':
            model.objective = 'Growth_Phenol'
            growth = 'Growth_Phenol'
            medium = {key:1000 for (key,value) in model.medium.items()}
            model.reactions.get_by_id('Growth_Glucose').upper_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            #remove all non-phenol carbon sources:
            medium["EX_glc__D_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            medium['EX_tag'] = 0
            medium["EX_phenol_e"] = sub_uptake_rate
            
            model.reactions.get_by_id('EX_glc__D_e').upper_bound = 0
            model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0

        elif substrate=='glucose':
            model.objective = 'Growth_Glucose'
            growth = 'Growth_Glucose'
            medium = {key:1000 for (key,value) in model.medium.items()}
            model.reactions.get_by_id('Growth_Phenol').upper_bound = 0
            model.reactions.get_by_id('Growth_Phenol').lower_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            
            #remove all non-glucose carbon sources:
            medium["EX_phenol_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            medium['EX_tag'] = 0
            medium["EX_glc__D_e"] = sub_uptake_rate

            model.reactions.get_by_id('EX_phenol_e').upper_bound = 0
            model.reactions.get_by_id('EX_phenol_e').lower_bound = 0
        else:
            print('Unknown substrate: Please choose among phenol and glucose')
        model.medium = medium
        try:
            pFBA_solution_all = cobra.flux_analysis.pfba(model)
        except:
            pFBA_solution_all = 0
            print("warning because of substrate cons rate for "+ index) 
        display(model.medium)
        
    return pFBA_solution_all

"""
    Make SPOT predictons.
    
        Parameters
        ----------
        model: cobrapy model.
        transcriptdf: data frame with transcript data.
        linename: linename in EDD used, WT or mutation. 
        substrate: carbon source used.
        sub_uptake_rate: uptake rate for carbon source used.
        
        Returns
        -------
        spotsol: SPOT prediction solution.
        
"""


#Function for EFLUX2 and SPOT Predictions:
def spot_pred(model, transcriptdf, linename, substrate, sub_uptake_rate=100):    
    with model:
        medium = model.medium
        if substrate=='phenol':
            model.objective = 'Growth_Phenol'
            medium = {key:1000 for (key,value) in model.medium.items()}
            model.reactions.get_by_id('Growth_Glucose').upper_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            medium["EX_glc__D_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            medium['EX_tag'] = 0
            #medium["EX_phenol_e"] = sub_uptake_rate
            #model.reactions.get_by_id('EX_phenol_e').upper_bound = -sub_uptake_rate
            #model.reactions.get_by_id('EX_phenol_e').lower_bound = -sub_uptake_rate
            model.reactions.get_by_id('EX_glc__D_e').upper_bound = 0
            model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').upper_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').lower_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').upper_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').lower_bound = 0
            model.reactions.get_by_id('EX_tag').upper_bound = 0
            model.reactions.get_by_id('EX_tag').lower_bound = 0
            #medium["EX_phenol_e"] = sub_uptake_rate
        elif substrate=='glucose':
            model.objective = 'Growth_Glucose'
            medium = {key:1000 for (key,value) in model.medium.items()}
            model.reactions.get_by_id('Growth_Phenol').upper_bound = 0
            model.reactions.get_by_id('Growth_Phenol').lower_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            #medium["EX_glc__D_e"] = sub_uptake_rate
            medium["EX_phenol_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            medium['EX_tag'] = 0
            model.reactions.get_by_id('EX_glc__D_e').upper_bound = 0#-sub_uptake_rate
            #model.reactions.get_by_id('EX_glc__D_e').lower_bound = -sub_uptake_rate
            model.reactions.get_by_id('EX_phenol_e').upper_bound = 0
            model.reactions.get_by_id('EX_phenol_e').lower_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').upper_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').lower_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').upper_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').lower_bound = 0
            model.reactions.get_by_id('EX_tag').upper_bound = 0
            model.reactions.get_by_id('EX_tag').lower_bound = 0
        else:
            print('Unknown substrate: Please choose among phenol and glucose')
        model.medium = medium
        spotsol = SPOT(model, transcriptdf)
    return spotsol

"""
    Make SPOT predictons for three replicates, average the solutions and calculate standard deviations.
    
        Parameters
        ----------
        model: cobrapy model.
        transcriptdf: data frame with transcript data.
        linename1: linename for replicate 1.
        linename2: linename for replicate 2.
        linename3: linename for replicate 3.
        substrate: carbon source used.
        
        Returns
        -------
        spotsol: SPOT prediction solution.
        spotsol_std: standard deviations of SPOT predictions. 
        
"""
# Function for predictions for three replicates and averaging the solutions and calculating the standard deviation:
def spot_pred_for_three_reps(model, transcriptdf, linename1, linename2, linename3, substrate):
    #call prediction functions for individual spot predictions for all 3 replicates:
    transdata_R1 = construct_trans_df(transcriptdf, linename1)
    transdata_R2 = construct_trans_df(transcriptdf, linename2)
    transdata_R3 = construct_trans_df(transcriptdf, linename3)
 
    spotsol_R1 = spot_pred(model, transdata_R1, linename1, substrate)
    spotsol_R2 = spot_pred(model, transdata_R2, linename2, substrate)
    spotsol_R3 = spot_pred(model, transdata_R3, linename3, substrate)
    
    #spot calculations:
    spotsol_R1_df = pd.DataFrame(spotsol_R1, columns=['fluxes'])
    spotsol_R2_df = pd.DataFrame(spotsol_R2, columns=['fluxes'])
    spotsol_R3_df = pd.DataFrame(spotsol_R3, columns=['fluxes'])
    spotsol_all = pd.concat([spotsol_R1_df, spotsol_R2_df, spotsol_R3_df], axis=1)
    
    spotsol = pd.DataFrame(spotsol_all.mean(axis=1), columns=['fluxes'])
    spotsol_std = spotsol_all.std(axis=1)
    
    return spotsol, spotsol_std

"""
    Make E-Flux2 predictons for three replicates, average the solutions and calculate standard deviations.
    
        Parameters
        ----------
        model: cobrapy model.
        transcriptdf: data frame with transcript data.
        linename1: linename for replicate 1.
        linename2: linename for replicate 2.
        linename3: linename for replicate 3.
        substrate: carbon source used.
        
        Returns
        -------
        eflux2sol: E-Flux2 prediction solution.
        eflux2sol_std: standard deviations of E-Flux2 predictions. 
        
"""

#Function for predictions for three replicates and averaging the solutions and calculating the standard deviation:
def eflux2_pred_for_three_reps(model, transcriptdf, linename1, linename2, linename3, substrate):
    #call prediction functions for individual EFLUX2 predictions for all 3 replicates:
    transdata_R1 = construct_trans_df(transcriptdf, linename1)
    transdata_R2 = construct_trans_df(transcriptdf, linename2)
    transdata_R3 = construct_trans_df(transcriptdf, linename3)
 
    eflux2sol_R1 = eflux2_pred(model, transdata_R1, linename1, substrate)
    eflux2sol_R2 = eflux2_pred(model, transdata_R2, linename2, substrate)
    eflux2sol_R3 = eflux2_pred(model, transdata_R3, linename3, substrate)
    
    #EFLUX2 calculations:
    eflux2sol_R1_df = pd.DataFrame(eflux2sol_R1.fluxes, columns=['fluxes'])
    eflux2sol_R2_df = pd.DataFrame(eflux2sol_R2.fluxes, columns=['fluxes'])
    eflux2sol_R3_df = pd.DataFrame(eflux2sol_R3.fluxes, columns=['fluxes'])
    eflux2sol_all = pd.concat([eflux2sol_R1_df, eflux2sol_R2_df, eflux2sol_R3_df], axis=1)

    eflux2sol = pd.DataFrame(eflux2sol_all.mean(axis=1), columns=['fluxes'])
    eflux2sol_std = eflux2sol_all.std(axis=1)
    
    return eflux2sol, eflux2sol_std

"""
    Reformat fluxes, such that reverse_ replaced with -1*flux.
    
        Parameters
        ----------
        reaction_id: reaction id that is being reformatted.
        solution: solution of fluxes.
        
        Returns
        -------
        reformatted solution.
        
"""
def get_flux_value(reaction_id, solution):
    if reaction_id.startswith('reverse_'):
        reaction_id = reaction_id.split('reverse_')[1]
        return -1*solution.fluxes[reaction_id]
    else:
        return solution.fluxes[reaction_id]    
"""
    Reformat standard deviations, such that for reverse_ fluxes returns solution of standard deviation.
    
        Parameters
        ----------
        reaction_id: reaction id that is being reformatted.
        solution: solution of fluxes.
        
        Returns
        -------
        reformatted solution.
        
"""
def get_std_value(reaction_id, solution):
    if reaction_id.startswith('reverse_'):
        reaction_id = reaction_id.split('reverse_')[1]
    return solution.stds[reaction_id]

"""
    Scales predicted fluxes and standard deviation in relation to carbon input.
    
        Parameters
        ----------
        observed_fluxes: data frame of observed fluxes
        predictions: data frame of predicted fluxes
        stdpredictions: data frame of standard deviations of predictions
        substrate: string of carbon input source
        method: string of prediction method
        
        Returns
        -------
        scalepred_fluxes: data frame of scaled fluxes.
        scalepred_stds: data frame of scaled standard deviations.
        
"""
# def scale_predictions(observed_fluxes, predictions, stdpredictions, substrate, method):
#     scalepred_stds = pd.DataFrame(index=stdpredictions.index, columns= ['stds'], dtype=np.float64)
#     scalepred_fluxes = pd.DataFrame(index=predictions.index, columns= ['fluxes'], dtype=np.float64)
#     if substrate == 'phenol':
#         phenoluptakerow = observed_fluxes[observed_fluxes['Pathway']=='Phenol Uptake']
#         sourceuptake = float(phenoluptakerow['Flux'])
#         scalepred_fluxes = predictions*(sourceuptake/(-1*predictions.loc['EX_phenol_e']))
#     elif substrate == 'glucose':
#         glucoseuptakerow = observed_fluxes[observed_fluxes['Pathway']=='Glucose Uptake']
#         sourceuptake = float(glucoseuptakerow['Flux'])
#         scalepred_fluxes = predictions*(sourceuptake/(-1*predictions.loc['EX_glc__D_e']))
#     else:   
#         print('Unknown Substrate')
#     for ind in stdpredictions.index:
#         if abs(stdpredictions.loc[ind,'stds'])<1e-5:
#             scalepred_stds.loc[ind,'stds'] = stdpredictions.loc[ind,'stds']
#         else:
#             scalepred_stds.loc[ind,'stds'] = (stdpredictions.loc[ind,'stds']/predictions.loc[ind,'fluxes'])*scalepred_fluxes.loc[ind, 'fluxes']
#     return scalepred_fluxes, scalepred_stds

"""
    Scales predicted fluxes and standard deviation in relation to carbon input.
    
        Parameters
        ----------
        observed_fluxes: data frame of observed fluxes
        predictions: data frame of predicted fluxes
        stdpredictions: data frame of standard deviations of predictions
        substrate: string of carbon input source
        method: string of prediction method
        
        Returns
        -------
        scalepred_fluxes: data frame of scaled fluxes.
        scalepred_stds: data frame of scaled standard deviations.
        
"""

def scale_predictions(observed_fluxes, predictions, stdpredictions, substrate, method):
    scalepred_stds = pd.DataFrame(index=stdpredictions.index, columns= ['stds'], dtype=np.float64)
    scalepred_fluxes = pd.DataFrame(index=predictions.index, columns= ['fluxes'], dtype=np.float64)
    if substrate == 'phenol':
        phenoluptakerow = observed_fluxes[observed_fluxes['Pathway']=='Phenol Uptake']
        sourceuptake = float(phenoluptakerow['Flux'])
        scale_factor = (sourceuptake/(-1*predictions.loc['EX_phenol_e']))
    elif substrate == 'glucose':
        glucoseuptakerow = observed_fluxes[observed_fluxes['Pathway']=='Glucose Uptake']
        sourceuptake = float(glucoseuptakerow['Flux'])
        scale_factor = (sourceuptake/(-1*predictions.loc['EX_glc__D_e']))
    else:   
        print('Unknown Substrate')
    scalepred_fluxes = predictions*scale_factor.values
    scalepred_stds = stdpredictions*scale_factor.values
    return scalepred_fluxes, scalepred_stds

"""
   Add predictions without standard deviations to 13C data frame.
    
        Parameters
        ----------
        df_13c: 13C data frame
        sol: solution of predicted fluxes
        method: string of prediction method
        strain: string of strain
        
        Returns
        -------
        df_13c: new 13C data frame that includes the predictions. 
        
"""
def add_pred_fluxes_to_13c_df_without_std(df_13C, sol, method, strain):
    FBA_fluxes = []
    for _, row in df_13C.iterrows():
        reactions = row['Forward Reactions']
        flux_value = 0
        for x in [x.strip('() ') for x in reactions.split(' or ')]:
            and_split = [y.strip('() ') for y in x.split(' and ')]
            flux_value += min([get_flux_value(v, sol) for v in and_split])
        FBA_fluxes.append(flux_value)

    df_13C[str(method) + ' ' + str(strain) + ' ' + 'Value'] = FBA_fluxes
    return df_13C

"""
   Add predictions and standard deviations to 13C data frame.
    
        Parameters
        ----------
        observed_fluxes: 13C data frame.
        predictions: solution of predicted fluxes.
        stdpredictions: standard deviations of predictions.
        substrate: string of carbon input source.
        method: string of prediction method.
        strain: string of strain.
        
        Returns
        -------
        observed_fluxes: new 13C data frame that includes the predictions. 
        
"""

def add_pred_fluxes_to_13c_df(observed_fluxes, predictions, stdpredictions, substrate, method, strain):
    predicted_fluxes = []
    predicted_stds = []
    scalepred_fluxes, scalepred_stds = scale_predictions(observed_fluxes, predictions, stdpredictions, substrate, method)
    for _, row in observed_fluxes.iterrows():
        reactions = row['Forward Reactions']
        flux_value_pred = 0
        std_value_pred = 0
        for x in [x.strip('() ') for x in reactions.split(' or ')]:
            and_split = [y.strip('() ') for y in x.split(' and ')]
            flux_value_pred += min([get_flux_value(v, scalepred_fluxes) for v in and_split])
            std_value_pred += min([get_std_value(v,scalepred_stds) for v in and_split])
        predicted_fluxes.append(flux_value_pred)
        predicted_stds.append(std_value_pred)

    observed_fluxes[str(method) + ' ' + str(strain) + ' Value'] = predicted_fluxes
    observed_fluxes[str(method) + ' ' + str(strain) + ' std Value'] = predicted_stds
    
    return observed_fluxes

"""
   Calculate mean absolute error.
    
        Parameters
        ----------
        observed: measurements.
        predicted: solution of predicted fluxes.
        
        Returns
        -------
        mean absolute error.
        
"""

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

"""
   Calculate root mean squared error.
    
        Parameters
        ----------
        observed: measurements.
        predicted: solution of predicted fluxes.
        
        Returns
        -------
        root mean squared error.
        
"""
    
def rmse_func(predicted, observed):
    return np.sqrt(((predicted - observed) ** 2).mean())



# def add_pred_fluxes_to_13c_df(observed_fluxes, predictions, stdpredictions, substrate, method, strain):
#     predicted_fluxes = []
#     predicted_stds = []
#     scalepred_fluxes, scalepred_stds = scale_predictions(observed_fluxes, predictions, stdpredictions, substrate, method)
#     for _, row in observed_fluxes.iterrows():
#         reactions = row['Forward Reactions']
#         flux_value_pred = 0
#         std_value_pred = 0
#         for x in [x.strip('() ') for x in reactions.split(' or ')]:
#             and_split = [y.strip('() ') for y in x.split(' and ')]
#             flux_value_pred += min([get_flux_value(v, scalepred_fluxes) for v in and_split])
#             std_value_pred += min([get_std_value(v,scalepred_stds) for v in and_split])
#         predicted_fluxes.append(flux_value_pred)
#         predicted_stds.append(std_value_pred)

#     observed_fluxes[str(method) + ' ' + str(strain) + ' Value'] = predicted_fluxes
#     observed_fluxes[str(method) + ' ' + str(strain) + ' std Value'] = predicted_stds
    
#     return observed_fluxes

"""
   Scale growth rates for comparison of predicted and observed growth rates. Scaled by multiplying with observed substrate uptake/predicted substrate uptake.
    
        Parameters
        ----------
        solgrowth: growth rate solution.
        soluptake: predicted substrate uptake
        sub_uptake_2comp: observed substrate uptake.
        
        Returns
        -------
        solgrowthnew: scaled growth rate.
        
"""

def scale_growth_to_sub(solgrowth, soluptake, sub_uptake_2comp):
    if soluptake<=1e-6:
        solgrowthnew = solgrowth
    else:
        factor = abs(sub_uptake_2comp/(-soluptake))
        solgrowthnew = solgrowth*factor
    return solgrowthnew