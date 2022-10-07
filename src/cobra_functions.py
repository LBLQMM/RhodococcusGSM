import cobra

# this function takes in a model, a substrate, and a substrate uptake rate
# and returns a cobrapy solution object
def get_FBA_solution(model, substrate, sub_uptake_rate=100, verbose=False):
    with model:
        # set the biomass reaction and media composition based on carbon source
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

# this function takes in a model, a substrate, and a substrate uptake rate
# and returns a cobrapy solution object
def get_pFBA_solution(model, substrate, sub_uptake_rate=100, verbose=False):
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
        
        if verbose:
            display(model.medium)
        
    return pFBA_solution_all