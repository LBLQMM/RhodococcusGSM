import numpy as np
from optlang.symbolics import add
"""
    Provides EFLUX2 predictions as explained in Machado et. al (2014) 
    
        Parameters
        ----------
        model : cobrapy model.
        Transcriptomics : pandas dataframe with transcriptomics data.
        
        Returns
        -------
        eflux2_sol as output from eflux2_model.optimize().
        
"""
#Code only works for GPRs written in disjunctive normal form (DNF). Majority of models have them in DNF but there are some exceptions. 

def EFlux2(model, Transcriptomics):
    eflux2_model = model.copy()
    # Parse GPR into a dict containing isozymes (separated by 'or')
    # Each isozyme has a set of subunits (separated by 'and')
    #'and' and 'or' can occur at the same time, or can occur by itself.
    gpr_dict = dict()
    for r in eflux2_model.reactions:
        if r.gene_reaction_rule:
            temp = set()
            for x in [x.strip('() ') for x in r.gene_reaction_rule.split(' or ')]:
                temp.add(frozenset(y.strip('() ') for y in x.split(' and ')))
            gpr_dict[r.id] = temp
    # Set the bounds using the transcriptomics data
    for r in eflux2_model.reactions:
        if r.gene_reaction_rule:
            #If a reaction R1 has the GPR of 'A and B', it would be parsed to { {A, B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B) ] ) = min(A, B).
            #If a reaction R1 has the GPR of 'A or B', it would be parsed to { {A}, {B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A), min(B) ] ) = sum( [A, B] ).
            #If a reaction R1 has the GPR of '(A and B) or (C and D)', it would be parsed to { {A, B}, {C, D} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B), min(C, D) ] ).
            t = np.sum([np.min([Transcriptomics.loc[g] if g in Transcriptomics.index 
                                else np.array([np.Inf]) for g in p])
                        for p in gpr_dict[r.id]])
            if r.lower_bound < 0.0:
                r.lower_bound = -t
            else:
                pass
            if r.upper_bound > 0.0:
                r.upper_bound = t
            else:
                pass
        else: 
            #When there is no GPR, the arbitrary bounds are removed. 
            #Common arbitrary bound value of 1000 for E.coli, might be different depending on the model, e.g., 99999.0 for iMM904 yeast model in BiGG
            if r.lower_bound <= -1000.0:
                r.lower_bound = -np.Inf
            if r.upper_bound >= 1000.0:
                r.upper_bound = np.Inf
    # solve FBA to calculate the maximum biomass
    eflux2_model.tolerance = 1e-9
    fba_sol = eflux2_model.optimize()
    print('FBA status', fba_sol.status)
    print('FBA solution', fba_sol.objective_value)
    # Constrain the biomass to the optimal value
    for r in eflux2_model.reactions:
        if r.objective_coefficient:
            r.lower_bound = fba_sol.objective_value
    # minimize the sum of squared flux values
    eflux2_model.objective = eflux2_model.problem.Objective(add([r.flux_expression**2 for r in eflux2_model.reactions]), direction='min')
    eflux2_sol = eflux2_model.optimize()
    print('EFlux2 status', eflux2_sol.status)
    print('EFlux2 solution', eflux2_sol.objective_value)
    # return eflux2 solution
    return(eflux2_sol)