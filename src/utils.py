import pandas as pd
import numpy as np
from ensemblemethods import SPOT, EFlux2
import cobra
from scipy import mean


#Transform data to dataframe with just index as gene identifiers and one column for values
#!!!!TODO: Generalize for multiple time points
#Function to construct df from E-Flux2 functions: Needs to be modified for multiple time points!!!!!
def construct_trans_df(transdata, linename):
    transdataWTPR1 = transdata[transdata['Line Name']==linename]
    transdataWTPR1new = transdataWTPR1.filter(['Count', 'Measurement Type'])
    transdataWTPR1new2 = transdataWTPR1new.set_index('Measurement Type')
    return transdataWTPR1new2



def FBA_pred(model, substrate, sub_uptake_rate=100, verbose=True):
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
        if verbose:
            display(model.medium)
    return fbasol



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
#             medium['EX_tag'] = 0
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
#             medium['EX_tag'] = 0
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

#Function for E-Flux2 and SPOT Predictions:
def eflux2_pred(model, transcriptdf, linename, substrate, sub_uptake_rate=100):  
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
            medium['EX_phenol_e'] = 1000
#             medium['EX_tag'] = 0
            
            model.reactions.get_by_id('EX_phenol_e').upper_bound = 0 #-sub_uptake_rate
#             model.reactions.get_by_id('EX_phenol_e').lower_bound = -sub_uptake_rate
            model.reactions.get_by_id('EX_glc__D_e').upper_bound = 0
            model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').upper_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').lower_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').upper_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').lower_bound = 0
#             model.reactions.get_by_id('EX_tag').upper_bound = 0
#             model.reactions.get_by_id('EX_tag').lower_bound = 0
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
#             medium['EX_tag'] = 0
            model.reactions.get_by_id('EX_glc__D_e').upper_bound = 0#-sub_uptake_rate
            #model.reactions.get_by_id('EX_glc__D_e').lower_bound = -sub_uptake_rate
            model.reactions.get_by_id('EX_phenol_e').upper_bound = 0
            model.reactions.get_by_id('EX_phenol_e').lower_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').upper_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').lower_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').upper_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').lower_bound = 0
#             model.reactions.get_by_id('EX_tag').upper_bound = 0
#             model.reactions.get_by_id('EX_tag').lower_bound = 0
        else:
            print('Unknown substrate: Please choose among phenol and glucose')
        model.medium = medium
        eflux2sol = EFlux2(model, transcriptdf)
    return eflux2sol

#Function for predictions for three replicates and averaging the solutions and calculating the standard deviation:
def eflux2_pred_for_three_reps(model, transcriptdf, linename1, linename2, linename3, substrate):
    
    #call prediction functions for individual E-Flux2 predictions for all 3 replicates:
    transdata_R1 = construct_trans_df(transcriptdf, linename1)
    transdata_R2 = construct_trans_df(transcriptdf, linename2)
    transdata_R3 = construct_trans_df(transcriptdf, linename3)
 
    print('running first replicate')
    eflux2sol_R1 = eflux2_pred(model, transdata_R1, linename1, substrate)
    print('running second replicate')
    eflux2sol_R2 = eflux2_pred(model, transdata_R2, linename2, substrate)
    print('running third replicate')
    eflux2sol_R3 = eflux2_pred(model, transdata_R3, linename3, substrate)
    
    #EFLUX2 calculations:
    eflux2sol_R1_df = pd.DataFrame(eflux2sol_R1.fluxes, columns=['fluxes'])
    eflux2sol_R2_df = pd.DataFrame(eflux2sol_R2.fluxes, columns=['fluxes'])
    eflux2sol_R3_df = pd.DataFrame(eflux2sol_R3.fluxes, columns=['fluxes'])
    eflux2sol_all = pd.concat([eflux2sol_R1_df, eflux2sol_R2_df, eflux2sol_R3_df], axis=1)

    eflux2sol = pd.DataFrame(eflux2sol_all.mean(axis=1), columns=['fluxes'])
    eflux2sol_std = eflux2sol_all.std(axis=1)
    
    return eflux2sol, eflux2sol_std

#Function for E-Flux2 and SPOT Predictions:
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
#             medium['EX_tag'] = 0
            #medium["EX_phenol_e"] = sub_uptake_rate
            #model.reactions.get_by_id('EX_phenol_e').upper_bound = -sub_uptake_rate
#             model.reactions.get_by_id('EX_phenol_e').lower_bound = -sub_uptake_rate
            model.reactions.get_by_id('EX_glc__D_e').upper_bound = 0
            model.reactions.get_by_id('EX_glc__D_e').lower_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').upper_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').lower_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').upper_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').lower_bound = 0
#             model.reactions.get_by_id('EX_tag').upper_bound = 0
#             model.reactions.get_by_id('EX_tag').lower_bound = 0
#             medium["EX_phenol_e"] = sub_uptake_rate
            medium["EX_phenol_e"] = np.inf
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
#             medium['EX_tag'] = 0
            model.reactions.get_by_id('EX_glc__D_e').upper_bound = 0#-sub_uptake_rate
            #model.reactions.get_by_id('EX_glc__D_e').lower_bound = -sub_uptake_rate
            model.reactions.get_by_id('EX_phenol_e').upper_bound = 0
            model.reactions.get_by_id('EX_phenol_e').lower_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').upper_bound = 0
            model.reactions.get_by_id('EX_guaiacol_e').lower_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').upper_bound = 0
            model.reactions.get_by_id('EX_vanlt_e').lower_bound = 0
#             model.reactions.get_by_id('EX_tag').upper_bound = 0
#             model.reactions.get_by_id('EX_tag').lower_bound = 0
        else:
            print('Unknown substrate: Please choose among phenol and glucose')
        model.medium = medium
        spotsol = SPOT(model, transcriptdf)
    return spotsol


# Function for predictions for three replicates and averaging the solutions and calculating the standard deviation:
def spot_pred_for_three_reps(model, transcriptdf, linename1, linename2, linename3, substrate):
    #call prediction functions for individual spot predictions for all 3 replicates:
    
    transdata_R1 = construct_trans_df(transcriptdf, linename1)
    transdata_R2 = construct_trans_df(transcriptdf, linename2)
    transdata_R3 = construct_trans_df(transcriptdf, linename3)
 
    print('running first replicate')
    spotsol_R1 = spot_pred(model, transdata_R1, linename1, substrate)
    print('running second replicate')
    spotsol_R2 = spot_pred(model, transdata_R2, linename2, substrate)
    print('running first third')
    spotsol_R3 = spot_pred(model, transdata_R3, linename3, substrate)
    
    #spot calculations:
    spotsol_R1_df = pd.DataFrame(spotsol_R1, columns=['fluxes'])
    spotsol_R2_df = pd.DataFrame(spotsol_R2, columns=['fluxes'])
    spotsol_R3_df = pd.DataFrame(spotsol_R3, columns=['fluxes'])
    spotsol_all = pd.concat([spotsol_R1_df, spotsol_R2_df, spotsol_R3_df], axis=1)
    
    spotsol = pd.DataFrame(spotsol_all.mean(axis=1), columns=['fluxes'])
    spotsol_std = spotsol_all.std(axis=1)
    
    return spotsol, spotsol_std

def get_flux_value(reaction_id, solution):
    if reaction_id.startswith('reverse_'):
        reaction_id = reaction_id.split('reverse_')[1]
        return -1*solution.fluxes[reaction_id]
    else:
        return solution.fluxes[reaction_id]    
    
def get_std_value(reaction_id, solution):
    if reaction_id.startswith('reverse_'):
        reaction_id = reaction_id.split('reverse_')[1]
    return solution.stds[reaction_id]
    
def scale_predictions(observed_fluxes, predictions, stdpredictions, substrate, method):
    scalepred_stds = pd.DataFrame(index=stdpredictions.index, columns= ['stds'], dtype=np.float64)
    scalepred_fluxes = pd.DataFrame(index=predictions.index, columns= ['fluxes'], dtype=np.float64)
    if substrate == 'phenol':
        phenoluptakerow = observed_fluxes[observed_fluxes['Pathway']=='Phenol Uptake']
        sourceuptake = float(phenoluptakerow['Flux'])
        scalepred_fluxes = predictions*(sourceuptake/(-1*predictions.loc['EX_phenol_e']))
    elif substrate == 'glucose':
        glucoseuptakerow = observed_fluxes[observed_fluxes['Pathway']=='Glucose Uptake']
        sourceuptake = float(glucoseuptakerow['Flux'])
        scalepred_fluxes = predictions*(sourceuptake/(-1*predictions.loc['EX_glc__D_e']))
    else:   
        print('Unknown Substrate')
    for ind in stdpredictions.index:
        if abs(stdpredictions.loc[ind,'stds'])<1e-5:
            scalepred_stds.loc[ind,'stds'] = stdpredictions.loc[ind,'stds']
        else:
            scalepred_stds.loc[ind,'stds'] = (stdpredictions.loc[ind,'stds']/predictions.loc[ind,'fluxes'])*scalepred_fluxes.loc[ind, 'fluxes']
    return scalepred_fluxes, scalepred_stds

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
    
#define function to calculate root mean squared error
def rmse_func(predicted, observed):
    return np.sqrt(((predicted - observed) ** 2).mean())

def scale_predictions(observed_fluxes, predictions, stdpredictions, substrate, method):
    scalepred_stds = pd.DataFrame(index=stdpredictions.index, columns= ['stds'], dtype=np.float64)
    scalepred_fluxes = pd.DataFrame(index=predictions.index, columns= ['fluxes'], dtype=np.float64)
    
    if substrate == 'phenol':
        phenoluptakerow = observed_fluxes[observed_fluxes['Pathway']=='Substrate Uptake']
        sourceuptake = float(phenoluptakerow['Flux'])
        scale_factor = (sourceuptake/(-1*predictions.loc['EX_phenol_e']))
        print('scale_factor', scale_factor)
    elif substrate == 'glucose':
        glucoseuptakerow = observed_fluxes[observed_fluxes['Pathway']=='Substrate Uptake']
        sourceuptake = float(glucoseuptakerow['Flux'])
        scale_factor = (sourceuptake/(-1*predictions.loc['EX_glc__D_e']))
        print('scale_factor', scale_factor)
    else:   
        print('Unknown Substrate')
    scalepred_fluxes = predictions*scale_factor.values
    scalepred_stds = stdpredictions*scale_factor.values
    return scalepred_fluxes, scalepred_stds

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


#For comparison of predicted and observed growth rates: scale predicted growth rate by multiplying with (observed substrate uptake / predicted substrate uptake)
def scale_growth_to_sub(solgrowth, soluptake, sub_uptake_2comp):
    if soluptake<=1e-6:
        solgrowthnew = solgrowth
    else:
        factor = abs(sub_uptake_2comp/(-soluptake))
        solgrowthnew = solgrowth*factor
    return solgrowthnew
    
def fba_solution_to_df(model, solution):
    fluxes = []
    for rxn_id, flux in solution.fluxes.items():
        fluxes.append({
            'reaction_id': rxn_id,
            'reaction_name': model.reactions.get_by_id(rxn_id).name,
            'reaction_reaction': model.reactions.get_by_id(rxn_id).reaction,
            'flux': flux
        })
    return pd.DataFrame(fluxes)
    
def eflux_solution_to_df(model, solution):
    fluxes = []
    for rxn_id, flux in solution.fluxes.items():
        fluxes.append({
            'reaction_id': rxn_id,
            'reaction_name': model.reactions.get_by_id(rxn_id).name,
            'reaction_reaction': model.reactions.get_by_id(rxn_id).reaction,
            'flux': flux
        })
    return pd.DataFrame(fluxes)

def spot_solution_to_df(model, solution):
    fluxes = []
    for rxn_id, flux in solution.items():
        fluxes.append({
            'reaction_id': rxn_id,
            'reaction_name': model.reactions.get_by_id(rxn_id).name,
            'reaction_reaction': model.reactions.get_by_id(rxn_id).reaction,
            'flux': flux
        })
    return pd.DataFrame(fluxes)

# def get_flux_value(reaction_id, solution):
#     if reaction_id.startswith('reverse_'):
#         reaction_id = reaction_id.split('reverse_')[1]
#         return -1*solution.fluxes[reaction_id]
#     else:
#         return solution.fluxes[reaction_id]
    
    
# def get_std_value(reaction_id, solution):
#     if reaction_id.startswith('reverse_'):
#         reaction_id = reaction_id.split('reverse_')[1]
#     return solution.stds[reaction_id]