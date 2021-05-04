import numpy as np
from optlang.symbolics import add

def SPOT(model, Transcriptomics):
    spot_model = model.copy()
    # Parse GPR into a dict containing isozymes (separated by 'or')
    # Each isozyme has a set of subunits (separated by 'and')
    gpr_dict = dict()
    for r in spot_model.reactions:
        if r.gene_reaction_rule:
            temp = set()
            for x in [x.strip('() ') for x in r.gene_reaction_rule.split(' or ')]:
                temp.add(frozenset(y.strip('() ') for y in x.split(' and ')))
            gpr_dict[r.id] = temp
    # calculate the lumped expression values for each reaction
    g_coef = dict()
    for r in spot_model.reactions:
        if r.gene_reaction_rule:
            t = np.sum([np.min([Transcriptomics.loc[g] if g in Transcriptomics.index 
                                else np.array([np.Inf]) for g in p])
                        for p in gpr_dict[r.id]])
            if t == np.Inf:
                g_coef[r.id] = 0
            else:
                g_coef[r.id] = t
        else:
            g_coef[r.id] = 0
    # constrain the norm of flux to be less than or equal to 1
    flux_norm = spot_model.problem.Constraint(add([r.forward_variable**2 for r in spot_model.reactions] + 
                                                  [r.reverse_variable**2 for r in spot_model.reactions]),
                                              ub = 1.0)
    spot_model.add_cons_vars(flux_norm)
    # maximize the pearson product-moment correlation
    spot_model.objective = spot_model.problem.Objective(add([g_coef[r.id]*r.forward_variable for r in spot_model.reactions] + 
                                                            [g_coef[r.id]*r.reverse_variable for r in spot_model.reactions]),
                                                        direction='max')
    # solve spot
    spot_model.tolerance = 1e-9
    spot_sol = spot_model.optimize()
    print('SPOT status', spot_sol.status)
    print('SPOT solution', spot_sol.objective_value)
    # return spot solution
    return(spot_sol)