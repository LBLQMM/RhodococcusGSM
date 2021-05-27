#Import python packages:
import numpy as np
from optlang.symbolics import add
import pandas as pd
import cplex

#################################################################
########Functions needed for all prediction methods##############
#################################################################
"""
    Creates dictionary of isozymes by parsing GPR:
    Parse GPR into a dict containing isozymes (separated by 'or'). Each isozyme has a set of subunits (separated by 'and') 'and' and 'or' can occur at the same time, or can occur by itself.
    
        Parameters
        ----------
        model : cobrapy model.
        
        
        Returns
        -------
        gpr_dict: dictionary with isozymes.
        
"""
#Code only works for GPRs written in disjunctive normal form (DNF). Majority of models have them in DNF but there are some exceptions. 

def create_gprdict(model):   
    gpr_dict = dict()
    for rxn in model.reactions:
        if rxn.gene_reaction_rule:
            temp = set()
            for x in [x.strip('() ') for x in rxn.gene_reaction_rule.split(' or ')]:
                temp.add(frozenset(y.strip('() ') for y in x.split(' and ')))
            gpr_dict[rxn.id] = temp
    return gpr_dict

"""
    Calculates bound value based on transcriptomics data for reactions in gene reaction rule
    
    NOTE: 
    If a reaction R1 has the GPR of 'A and B', it would be parsed to { {A, B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B) ] ) = min(A, B).
    If a reaction R1 has the GPR of 'A or B', it would be parsed to { {A}, {B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A), min(B) ] ) = sum( [A, B] ).
    If a reaction R1 has the GPR of '(A and B) or (C and D)', it would be parsed to { {A, B}, {C, D} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B), min(C, D) ] ).
    
        Parameters
        ----------
        model : cobrapy model.
        Transcriptomics : pandas dataframe with transcriptomics data.Data frame has gene identifiers as index and just one column with transcript values.  
        rxn : cobrapy model reaction
        
        
        Returns
        -------
        transscript bound value: float.
"""

def findtransboundval_forgprrxns(model, Transcriptomics,rxn):
    finaltransval = 0
    for parallel_gene in create_gprdict(model)[rxn.id]:
        transvals = []
        for gene in parallel_gene:
            if gene in Transcriptomics.index:
                transvals.append(Transcriptomics.loc[gene].values)
            else:
                transvals.append(np.Inf)
            mintransval=np.min(transvals)
        finaltransval = finaltransval + mintransval
    return finaltransval

#############################################
##################EFLUX2#######################
############################################
"""
    Provides EFLUX2 predictions as explained in Machado et. al (2014) 
    
        Parameters
        ----------
        model : cobrapy model.
        Transcriptomics : pandas dataframe with transcriptomics data.Data frame has gene identifiers as index and just one column with transcript values.  
        
        Returns
        -------
        eflux2_sol as output from eflux2_model.optimize().
        
"""
def EFlux2(model, Transcriptomics):
    eflux2_model = model.copy()
    
    # Set the bounds using the transcriptomics data    
    for rxn in eflux2_model.reactions:
        if rxn.gene_reaction_rule:
                
            if rxn.lower_bound < 0.0:
                rxn.lower_bound = -findtransboundval_forgprrxns(model, Transcriptomics,rxn)
            else:
                pass
            if rxn.upper_bound > 0.0:
                rxn.upper_bound = findtransboundval_forgprrxns(model, Transcriptomics,rxn)
            else:
                pass
        else:
            """When there is no GPR, the arbitrary bounds are removed. 
            Common arbitrary bound value of 1000 for E.coli, might be different depending on the model, e.g., 99999.0 for iMM904 yeast model in BiGG"""
            if rxn.lower_bound <= -1000:
                rxn.lower_bound = -np.Inf
            if rxn.upper_bound >= 1000:
                rxn.upper_bound = np.Inf 
    
    # Solve FBA to calculate the maximum biomass
    eflux2_model.tolerance = 1e-9
    fba_sol = eflux2_model.optimize()
    print('FBA status', fba_sol.status)
    print('FBA solution', fba_sol.objective_value)
    display(eflux2_model.objective)
    
    # Constrain the biomass to the optimal value
    for r in eflux2_model.reactions:
        if r.objective_coefficient:
            r.lower_bound = fba_sol.objective_value
            
    # Minimize the sum of squared flux values
    """Note: Because of quadratic objective still have to use cplex objective formulation.
    Optlang does not support quadratic type of constraints and objectives yet."""
    eflux2_model.objective = eflux2_model.problem.Objective(add([rxn.flux_expression**2 for rxn in eflux2_model.reactions]), direction='min')
    eflux2_sol = eflux2_model.optimize()
    print('EFlux2 status', eflux2_sol.status)
    print('EFlux2 solution', eflux2_sol.objective_value)
    
    return eflux2_sol


#############################################
##################SPOT#######################
############################################

"""
    Provides SPOT predictions as explained in Machado et. al (2014) using cplex python library. If optlang enhanced such that quadratic constraints are supported, we should be able to switch to formulation without using cplex (spot.py code). 
    Calls create_gprdict function.
    
        Parameters
        ----------
        model : cobrapy model.
        Transcriptomics : pandas dataframe with transcriptomics data. Data frame has gene identifiers as index and just one column with transcript values.  
        
        Returns
        -------
        sol as output from optimization of SPOT model via CPLEX formulation.
        
"""
def SPOT(model, Transcriptomics):
    
    mets = [met.id for met in model.metabolites]
    rxns = [rxn.id for rxn in model.reactions]
    nrow = len(mets)
    ncol = len(rxns)

    rev_rxns = ['rev_'+rxn.id for r in model.reactions if rxn.reversibility]
    rev_ncol = len(rev_rxns)

    gpr_dict = dict()
    for r in model.reactions:
        """Parse GPR into a dict containing isozymes (separated by 'or')
        # Each isozyme has a set of subunits (separated by 'and')
        #'and' and 'or' can occur at the same time, or can occur by itself."""

    lb = [0.0 if r.reversibility else r.lower_bound for r in model.reactions] + [0.0 for r in model.reactions if r.reversibility]
    ub = [r.upper_bound for r in model.reactions] + [-r.lower_bound for r in model.reactions if r.reversibility]
        
    c = []
    for r in model.reactions:
        if r.gene_reaction_rule:
        #If a reaction R1 has the GPR of 'A and B', it would be parsed to { {A, B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B) ] ) = min(A, B).
        #If a reaction R1 has the GPR of 'A or B', it would be parsed to { {A}, {B} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A), min(B) ] ) = sum( [A, B] ).
        #If a reaction R1 has the GPR of '(A and B) or (C and D)', it would be parsed to { {A, B}, {C, D} } in gpr_dict['R1']. Then t for R1 would be sum( [ min(A, B), min(C, D) ] ).
        
#             t = np.sum([np.min([Transcriptomics.loc[g] if g in Transcriptomics.index 
#                                 else np.array([np.Inf]) for g in p])
#                         for p in create_gprdict(model)[r.id]])
            transboundval = findtransboundval_forgprrxns(model, Transcriptomics,rxn)
            if transboundval == np.Inf:
                transboundval = 0
            c.append(transboundval)
        else:
            c.append(0.0)
    for r in model.reactions:
        if r.reversibility:
            if r.gene_reaction_rule:
                transboundval = findtransboundval_forgprrxns(model, Transcriptomics,rxn)
                if transboundval == np.Inf:
                    transboundval = 0
                c.append(transboundval)
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