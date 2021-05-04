import numpy as np
import pandas as pd
import cplex

"""
    Provides SPOT predictions as explained in Machado et. al (2014) using cplex python library. If optlang enhanced such that quadratic constraints are supported, we should be able to switch to formulation without using cplex (spot.py code). 
    
        Parameters
        ----------
        model : cobrapy model.
        Transcriptomics : pandas dataframe with transcriptomics data.
        
        Returns
        -------
        sol as output from optimization of SPOT model via CPLEX formulation.
        
"""
def SPOT(model, Transcriptomics):
    
    mets = [m.id for m in model.metabolites]
    rxns = [r.id for r in model.reactions]
    nrow = len(mets)
    ncol = len(rxns)

    rev_rxns = ['rev_'+r.id for r in model.reactions if r.reversibility]
    rev_ncol = len(rev_rxns)

    gpr_dict = dict()
    for r in model.reactions:
        # Parse GPR into a dict containing isozymes (separated by 'or')
        # Each isozyme has a set of subunits (separated by 'and')
        #'and' and 'or' can occur at the same time, or can occur by itself.
        if r.gene_reaction_rule:
            temp = set()
            for x in [x.strip('() ') for x in r.gene_reaction_rule.split(' or ')]:
                temp.add(frozenset(y.strip('() ') for y in x.split(' and ')))
            gpr_dict[r.id] = temp

    lb = [0.0 if r.reversibility else r.lower_bound for r in model.reactions] + [0.0 for r in model.reactions if r.reversibility]
    ub = [r.upper_bound for r in model.reactions] + [-r.lower_bound for r in model.reactions if r.reversibility]
        
    c = []
    for r in model.reactions:
        if r.gene_reaction_rule:
        #If a reaction R1 has the GPR of 'A and B', it would be parsed to { {A, B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B) ] ) = min(A, B).
        #If a reaction R1 has the GPR of 'A or B', it would be parsed to { {A}, {B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A), min(B) ] ) = sum( [A, B] ).
        #If a reaction R1 has the GPR of '(A and B) or (C and D)', it would be parsed to { {A, B}, {C, D} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B), min(C, D) ] ).
            t = np.sum([np.min([Transcriptomics.loc[g] if g in Transcriptomics.index 
                                else np.array([np.Inf]) for g in p])
                        for p in gpr_dict[r.id]])
            if t == np.Inf:
                t = 0
            c.append(t)
        else:
            c.append(0.0)
    for r in model.reactions:
        if r.reversibility:
            if r.gene_reaction_rule:
                t = np.sum([np.min([Transcriptomics.loc[g] if g in Transcriptomics.index 
                                    else np.array([np.Inf]) for g in p])
                            for p in gpr_dict[r.id]])
                if t == np.Inf:
                    t = 0
                c.append(t)
            else:
                c.append(0.0)

    SPOT = cplex.Cplex()
    SPOT.set_results_stream(None)
    SPOT.parameters.simplex.tolerances.optimality.set(1e-9)
    SPOT.parameters.simplex.tolerances.feasibility.set(1e-9)

    SPOT.linear_constraints.add(rhs=[0]*nrow, senses='E'*nrow, names=mets)
    SPOT.variables.add(obj=c, lb=lb, ub=ub, names=rxns+rev_rxns)
    for r in model.reactions:
        for m, v in r.metabolites.items():
            SPOT.linear_constraints.set_coefficients(m.id, r.id, v)
    for r in model.reactions:
        if r.reversibility:
            for m, v in r.metabolites.items():
                SPOT.linear_constraints.set_coefficients(m.id, 'rev_'+r.id, -v)
    SPOT.quadratic_constraints.add(quad_expr=[rxns+rev_rxns, rxns+rev_rxns, [1]*len(c)],
                                   sense='L', rhs=1.0, name='L2norm')#L indicating <=
    SPOT.objective.set_sense(SPOT.objective.sense.maximize)
    SPOT.solve()
    SPOT_sol = SPOT.solution.get_objective_value()

    sol = type('',(),{})()
    temp = pd.Series(data=SPOT.solution.get_values(), index=rxns+rev_rxns)
    flux = temp.loc[rxns]
    flux_rev = temp.loc[rev_rxns]
    for r in model.reactions:
        if r.reversibility:
            flux.loc[r.id] = flux.loc[r.id] - flux_rev.loc['rev_'+r.id]
    sol = flux
    sol.objective_value = SPOT.solution.get_objective_value()
    sol.status = SPOT.solution.get_status_string()
    
    return(sol)