import numpy as np
import cobra
from t_to_f_functions import SPOT, EFlux2

# this function takes in a model, a substrate, and a substrate uptake rate
# and returns a cobrapy solution object
def get_FBA_solution(model, substrate, substrate_uptake_rate=100, verbose=False):
    with model:
        medium = model.medium
        
        # set the biomass reaction and media composition based on carbon source
        if substrate=='phenol':
            # block glucose biomass reaction and the default CarveMe biomass reactions
            model.reactions.get_by_id('Growth_Glucose').upper_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            
            # make maximizing the phenol biomass reaction the objective function
            model.objective = 'Growth_Phenol'
            
            # first, set all media components to surplus levels
            medium = {key:1000 for (key,value) in model.medium.items()}
            
            # set the phenol uptake rate the specified value
            medium["EX_phenol_e"] = substrate_uptake_rate
            
            # remove all non-phenol carbon sources
            medium["EX_glc__D_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            
        elif substrate=='glucose':
            # block phenol biomass reaction and the default CarveMe biomass reactions
            model.reactions.get_by_id('Growth_Phenol').upper_bound = 0
            model.reactions.get_by_id('Growth_Phenol').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
        
            # make maximizing the glucose biomass reaction the objective function
            model.objective = 'Growth_Glucose'
            
            # first, set all media components to surplus levels
            medium = {key:1000 for (key,value) in model.medium.items()}
            
            # set the glucose uptake rate the specified value
            medium["EX_glc__D_e"] = substrate_uptake_rate
            
            #remove all non-glucose carbon sources:
            medium["EX_phenol_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
        else:
            print('Unknown substrate: Please choose among phenol and glucose')
            
        # conditionally, display the medium
        if verbose:
            display(model.medium)
            
        # update model medium and solve the FBA problem
        model.medium = medium
        fba_solution = model.optimize()
            
    return fba_solution

# this function takes in a model, a substrate, and a substrate uptake rate
# and returns a cobrapy solution object
def get_pFBA_solution(model, substrate, substrate_uptake_rate=100, verbose=False):
    with model:
        medium = model.medium
        
        # set the biomass reaction and media composition based on carbon source
        if substrate=='phenol':
            # block glucose biomass reaction and the default CarveMe biomass reactions
            model.reactions.get_by_id('Growth_Glucose').upper_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            
            # make maximizing the phenol biomass reaction the objective function
            model.objective = 'Growth_Phenol'
            
            # first, set all media components to surplus levels
            medium = {key:1000 for (key,value) in model.medium.items()}
            
            # set the phenol uptake rate the specified value
            medium["EX_phenol_e"] = substrate_uptake_rate
            
            # remove all non-phenol carbon sources
            medium["EX_glc__D_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            
        elif substrate=='glucose':
            # block phenol biomass reaction and the default CarveMe biomass reactions
            model.reactions.get_by_id('Growth_Phenol').upper_bound = 0
            model.reactions.get_by_id('Growth_Phenol').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
        
            # make maximizing the glucose biomass reaction the objective function
            model.objective = 'Growth_Glucose'
            
            # first, set all media components to surplus levels
            medium = {key:1000 for (key,value) in model.medium.items()}
            
            # set the glucose uptake rate the specified value
            medium["EX_glc__D_e"] = substrate_uptake_rate
            
            #remove all non-glucose carbon sources:
            medium["EX_phenol_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
        else:
            print('Unknown substrate: Please choose among phenol and glucose')
        
        # conditionally, display the medium
        if verbose:
            display(model.medium)
        
        # update model medium and solve the FBA problem
        model.medium = medium
        pFBA_solution = cobra.flux_analysis.pfba(model)

    return pFBA_solution

# this function takes in a model, transcript dataframe, a substrate string, 
# and a substrate uptake rate. It returns a cobrapy solution object
def get_EFlux2_solution(model, transcript_df, substrate, substrate_uptake_rate=100):  
    with model:
        # set the biomass reaction and media composition based on carbon source
        medium = model.medium
        if substrate == 'phenol':
            # block glucose biomass reaction and the default CarveMe biomass reactions
            model.reactions.get_by_id('Growth_Glucose').upper_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            
            # make maximizing the phenol biomass reaction the objective function
            model.objective = 'Growth_Phenol'
            
            # first, set all media components to surplus levels
            medium = {key:1000 for (key,value) in model.medium.items()}
            
            # set the phenol uptake rate the specified value
            medium["EX_phenol_e"] = substrate_uptake_rate
            
            # remove all non-phenol carbon sources
            medium["EX_glc__D_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            
            # allow for excess phenol
            medium['EX_phenol_e'] = 1000

        elif substrate == 'glucose':
            # block phenol biomass reaction and the default CarveMe biomass reactions
            model.reactions.get_by_id('Growth_Phenol').upper_bound = 0
            model.reactions.get_by_id('Growth_Phenol').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            # model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
        
            # make maximizing the glucose biomass reaction the objective function
            model.objective = 'Growth_Glucose'
            
            # first, set all media components to surplus levels
            medium = {key:1000 for (key,value) in model.medium.items()}
            
            #remove all non-glucose carbon sources:
            medium["EX_phenol_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            
            # allow for excess glucose
            medium['EX_glc__D_e'] = 1000
        else:
            print('Unknown substrate: Please choose among phenol and glucose')
            
        # update model medium and solve the E-Flux2 problem
        model.medium = medium
        EFlux2_solution = EFlux2(model, transcript_df)
        
    return EFlux2_solution

# this function takes in a model, transcript dataframe, a substrate string, 
# and a substrate uptake rate. It returns a cobrapy solution object  
def get_SPOT_solution(model, transcript_df, substrate, substrate_uptake_rate=100): 
    with model:
        medium = model.medium
        if substrate == 'phenol':
            # block glucose biomass reaction and the default CarveMe biomass reactions
            model.reactions.get_by_id('Growth_Glucose').upper_bound = 0
            model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            
            # make maximizing the phenol biomass reaction the objective function
            model.objective = 'Growth_Phenol'
            
            # first, set all media components to surplus levels
            medium = {key:1000 for (key,value) in model.medium.items()}
            
            # remove all non-phenol carbon sources
            medium["EX_glc__D_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            
            # allow for excess phenol
            medium["EX_phenol_e"] = np.inf
            
        elif substrate == 'glucose':
        
            # block phenol biomass reaction and the default CarveMe biomass reactions
            model.reactions.get_by_id('Growth_Phenol').upper_bound = 0
            model.reactions.get_by_id('Growth_Phenol').lower_bound = 0
            model.reactions.get_by_id('Growth').upper_bound = 0
            model.reactions.get_by_id('Growth').lower_bound = 0
            # model.reactions.get_by_id('Growth_Glucose').lower_bound = 0
        
            # make maximizing the glucose biomass reaction the objective function
            model.objective = 'Growth_Glucose'
            
            # first, set all media components to surplus levels
            medium = {key:1000 for (key,value) in model.medium.items()}
            
            #remove all non-glucose carbon sources:
            medium["EX_phenol_e"] = 0
            medium['EX_guaiacol_e'] = 0
            medium['EX_vanlt_e'] = 0
            
            # allow for excess glucose
            medium['EX_glc__D_e'] = np.inf
        else:
            print('Unknown substrate: Please choose among phenol and glucose')
        model.medium = medium
        spot_solution = SPOT(model, transcript_df)
    return spot_solution