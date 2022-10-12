import pandas as pd

# This function takes in a cobrapy model and a cobrapy solution, and 
# returns a dataframe with all the fluxes labeled by id, name, and reaction.
def cobra_solution_to_df(model, solution):
    fluxes = []
    for rxn_id, flux in solution.fluxes.items():
        fluxes.append({
            'reaction_id': rxn_id,
            'reaction_name': model.reactions.get_by_id(rxn_id).name,
            'reaction_reaction': model.reactions.get_by_id(rxn_id).reaction,
            'flux': flux
        })
    return pd.DataFrame(fluxes)

# This function takes in a reaction string and a cobra 
# solution, and it returns the flux value for that reaction
def reaction_id_to_flux(reaction_id, solution):
    if reaction_id.startswith('reverse_'):
        reaction_id = reaction_id.split('reverse_')[1]
        return -1*solution.fluxes[reaction_id]
    else:
        return solution.fluxes[reaction_id]    

# This function takes in a row from a central flux dataframe and a string 
# of reaction ids, and it returns a flux value that corrosponds to the reaction ids.
def reaction_ids_to_flux_value(solution, reaction_ids):
    total_flux = 0
    for x in [x.strip('() ') for x in reaction_ids.split(' or ')]:
        and_split = [y.strip('() ') for y in x.split(' and ')]
        total_flux += min([reaction_id_to_flux(v, solution) for v in and_split])
            
    return total_flux
    
# This function takes in a 13C flux dataframe and a cobra solution. It determines the 
# flux value using a reaction ids string that maps the 13C-MFA reaction to GSM reactions/
# It uses the reaction_ids_to_flux_value function to calculate the GSM flux.
def add_column_to_13C_flux_df(central_flux_df, solution, column_name):
    updated_df = central_flux_df.copy()
    
    # create a blank list to hold values to add to the column
    column_values = []
    
    # loop over rows in the central flux dataframe
    for _, row in central_flux_df.iterrows():
        # Add the flux value for each row to the column values list
        reaction_ids = row['Reaction Ids']
        flux_value = reaction_ids_to_flux_value(solution, reaction_ids)
        column_values.append(flux_value)

    # add the column to the dataframe
    updated_df[column_name] = column_values
    
    return updated_df