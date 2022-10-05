# -*- coding: utf-8 -*-

"""Provide Strain Design Algorithms"""

import sys
import numpy as np
from optlang.symbolics import add
from cobra.flux_analysis import find_essential_reactions, find_blocked_reactions, flux_variability_analysis
from cobra.util.solver import interface_to_str

def optforce(model, biomass, target, maximum_knockout, minimum_growth_fraction, minimum_target_fraction, peripheral_reactions,
             high_order_must=False, number_of_solutions=10, time_limit=60, threads=None, verbose=False):
    """
    Compute multiple solutions of OptForce
   
    Parameters
    ----------
    model : cobra.Model
        The model to compute an OptForce solution for
    biomass : str
        Biomass reaction id
    target : str
        Target reaction id
    maximum_knockout : int
        Maximum knockout limit
    minimum_growth_fraction : float
        Minimum growth requirement as a fraction of maximum
    minimum_target_fraction : float
        Minimum target production as a fraction of maximum
    peripheral_reactions : array
        Reactions to NOT consider in optforce solutions
    high_order_must : bool, optional
        Whether to identify high order must sets (default False)
    number_of_solutions : int, optional
        Maximum number of OptForce solutions (default 10)
    time_limit : any nonnegative value, optional
        Maximum time in seconds for each optimization (default 60)
    threads : int, optional
        Maximum number of parallel threads (default None)
    verbose : bool, optional
        Whether to print the solver log (default False)

    Returns
    -------
    solution
        A dictionary containing the must set solution and force set solution
    See Also
    --------
    find_must : find the must set solution
    find_force : find the force set solution
    """
    # Check the solver
    solver_str = interface_to_str(model.problem)
    if not solver_str in ['gurobi','cplex']:
        raise ValueError('%s is not supported. Use gurobi or cplex.' % solver_str)
    # Check biomass and target reactions
    if not biomass in model.reactions:
         raise ValueError('%s is not in the model.' % biomass)
    if not target in model.reactions:
        raise ValueError('%s is not in the model.' % target)
    # Check minimum biomass and target values
    with model:
        model.objective = biomass
        sol = model.optimize()
        maximum_growth = sol.objective_value
        if verbose:
            print('Maximum growth rate:', round(maximum_growth,3))
        minimum_growth = minimum_growth_fraction*maximum_growth
        model.objective = target
        sol = model.optimize()
        maximum_target = sol.objective_value
        if verbose:
            print('Maximum target production:', round(maximum_target,3))
        if maximum_target < 1e-9:
            raise ValueError('Maximum target production is less than 1e-9.')
        minimum_target = minimum_target_fraction*maximum_target
        model.objective = biomass
        model.reactions.get_by_id(target).lower_bound = minimum_target
        sol = model.optimize()
        if verbose:
            print('Maximum growth with minimum target production:', round(sol.objective_value,3))
            print('Minimum target production:', round(minimum_target,3))
    # Get the solver name and set parameters
    solver_str = interface_to_str(model.problem)
    if solver_str == 'gurobi':
        model.solver.problem.params.FeasibilityTol = 1e-9
        model.solver.problem.params.OptimalityTol = 1e-9
    elif solver_str == 'cplex':
        model.solver.problem.parameters.simplex.tolerances.feasibility.set(1e-9)
        model.solver.problem.parameters.simplex.tolerances.optimality.set(1e-9)
    # Find essential and blocked reactions
    model.tolerance = 1e-9
    essential_reactions = [r.id for r in find_essential_reactions(model)]
    blocked_reactions = find_blocked_reactions(model)
    # Remove blocked reactions
    model.remove_reactions(blocked_reactions, remove_orphans=True)
    # Calculate flux ranges for the parent strain
    with model:
        #model.reactions.get_by_id(target).upper_bound = 0.05*maximum_target
        fva_solution_WT = flux_variability_analysis(model, fraction_of_optimum=0.9)
    # Calculate flux ranges with the minimum growth and target production
    with model:
        model.reactions.get_by_id(biomass).lower_bound = minimum_growth
        model.reactions.get_by_id(target).lower_bound = minimum_target
        fva_solution_MT = flux_variability_analysis(model, fraction_of_optimum=0.0)
    # Find Must set
    excluded_reactions = peripheral_reactions + [biomass,'ATPM'] + [r.id for r in model.reactions if r.boundary or not r.gene_reaction_rule]
    print('excluded reactions')
    print(excluded_reactions)
    
    included_reactions = [r.id for r in model.reactions if r.id not in excluded_reactions]
    print('included reactions')
    print(included_reactions)
    
    
    Must_solution, MustL_set, MustU_set = find_must(model, biomass, target, minimum_growth, minimum_target, fva_solution_WT,
                                                    excluded_reactions, high_order_must, time_limit, threads, verbose)
    # Find Force set
    maximum_L = 2
    maximum_U = 2
    excluded_reactions = excluded_reactions + essential_reactions + list(MustL_set) + list(MustU_set)
    knockout_penalty = 0.01*maximum_target
    Force_solution  = find_force(model, biomass, target, maximum_knockout, minimum_growth,
                                 fva_solution_MT, MustL_set, MustU_set, maximum_L, maximum_U,
                                 excluded_reactions, knockout_penalty,
                                 number_of_solutions, time_limit, threads, verbose)
    solution = {'Must_solution': Must_solution, 'Force_solution': Force_solution}
    return solution

def find_must(model, biomass, target, minimum_growth, minimum_target, fva_solution_WT,
              excluded_reactions=[], high_order_must=False, time_limit=60, threads=None, verbose=False):
    """
    Identify must set
   
    Parameters
    ----------
    model : cobra.Model
        The model to compute an OptForce solution for
    biomass : str
        Biomass reaction id
    target : str
        Target reaction id
    minimum_growth : float
        Minimum growth requirement
    minimum_target : float
        Minimum target production
    fva_solution_WT : ouput from cobra.flux_analysis.flux_variability_analysis
        Minimum and maximum possible flux value for each reaction in WT
    exlcluded_reactions : list, optional
        List of reactions excluded from must set (default [])
    high_order_must : bool, optional
        Whether to identify high order must sets (default False)
    time_limit : any nonnegative value, optional
        Maximum time in seconds for each optimization (default 60)
    threads : int, optional
        Maximum number of parallel threads (default None)
    verbose : bool, optional
        Whether to print the solver log (default False)

    Returns
    -------
    solution
        A dictionary containing the must set solution
    MustL
        A set of mustL reactions
    MustU
        A set of mustU reactions
    See Also
    --------
    add_must : add variables, constraints, and objective for must set identification
    """
    # Make a copy of the model and add variables, constraints, and objective 
    optmodel = model.copy()
    solver_str = interface_to_str(optmodel.problem)
    maximum_L = 1
    maximum_U = 1
    add_must(optmodel, biomass, target, minimum_growth, minimum_target, fva_solution_WT, maximum_L, maximum_U,
             excluded_reactions)
    # Set solver specific parameters
    if solver_str == 'gurobi':
        if time_limit:
            optmodel.solver.problem.params.TimeLimit = time_limit
        if threads:
            optmodel.solver.problem.params.Threads = threads
        if verbose:
            optmodel.solver.problem.params.OutputFlag = 1
    elif solver_str == 'cplex':
        if time_limit:
            optmodel.solver.problem.parameters.timelimit.set(time_limit)
        if threads:
            optmodel.solver.problem.parameters.threads.set(threads)
        if verbose:
            optmodel.solver.problem.set_results_stream(sys.stdout)
    # Solve Must problems and read solution
    # optlang's optimize function involves extra steps
    # Gurobi interface resets if the existing status is not optimal
    variable_names = [v.name for v in optmodel.variables]
    # MustL
    if verbose:
        print('MustL')
    optmodel.constraints.get('maximum_U_constraint').ub = 0
    optmodel.constraints.get('maximum_L_constraint').ub = 1
    MustL = []
    temp = 0
    while temp < 100:
        optmodel.optimize()
        if optmodel.solver.status == 'optimal':
            if optmodel.solver.objective.value < 0.01:
                break
            variable_values = [v.primal for v in optmodel.variables]
            for j in range(len(variable_names)):
                if variable_names[j].startswith("yL_") and variable_values[j]:
                    if verbose:
                        print('L', variable_names[j][3:])
                    MustL.append(variable_names[j][3:])
                    optmodel.variables.get(variable_names[j]).ub = 0
            temp = temp + 1
        else:
            break
    # MustU
    if verbose:
        print('MustU')
    optmodel.constraints.get('maximum_L_constraint').ub = 0
    optmodel.constraints.get('maximum_U_constraint').ub = 1
    MustU = []
    temp = 0
    while temp < 100:
        optmodel.optimize()
        if optmodel.solver.status == 'optimal':
            if optmodel.solver.objective.value < 0.01:
                break
            variable_values = [v.primal for v in optmodel.variables]
            for j in range(len(variable_names)):
                if variable_names[j].startswith("yU_") and variable_values[j]:
                    if verbose:
                        print('U', variable_names[j][3:])
                    MustU.append(variable_names[j][3:])
                    optmodel.variables.get(variable_names[j]).ub = 0
            temp = temp + 1
        else:
            break
    # MustLL
    if verbose:
        print('MustLL')
    optmodel.constraints.get('maximum_U_constraint').ub = 0
    optmodel.constraints.get('maximum_L_constraint').ub = 2
    MustLL = []
    temp = 0
    while temp < 100:
        optmodel.optimize()
        if optmodel.solver.status == 'optimal':
            if optmodel.solver.objective.value < 0.01:
                break
            variable_values = [v.primal for v in optmodel.variables]
            MustLL.append([])
            for j in range(len(variable_names)):
                if variable_names[j].startswith("yL_") and variable_values[j]:
                    if verbose:
                        print('L', variable_names[j][3:])
                    MustLL[temp].append(variable_names[j][3:])
            intcut = optmodel.problem.Constraint(add([optmodel.variables.get('yL_'+r) for r in MustLL[temp]]),
                                                 ub=1, name='intcut_LL'+str(temp))
            optmodel.add_cons_vars(intcut)
            temp = temp + 1
        else:
            break
    # MustUU
    if verbose:
        print('MustUU')
    optmodel.constraints.get('maximum_U_constraint').ub = 2
    optmodel.constraints.get('maximum_L_constraint').ub = 0
    MustUU = []
    temp = 0
    while temp < 100:
        optmodel.optimize()
        if optmodel.solver.status == 'optimal':
            if optmodel.solver.objective.value < 0.01:
                break
            variable_values = [v.primal for v in optmodel.variables]
            MustUU.append([])
            for j in range(len(variable_names)):
                if variable_names[j].startswith("yU_") and variable_values[j]:
                    if verbose:
                        print('U', variable_names[j][3:])
                    MustUU[temp].append(variable_names[j][3:])
            intcut = optmodel.problem.Constraint(add([optmodel.variables.get('yU_'+r) for r in MustUU[temp]]),
                                                 ub=1, name='intcut_UU'+str(temp))
            optmodel.add_cons_vars(intcut)
            temp = temp + 1
        else:
            break
    # MustLU
    if verbose:
        print('MustLU')
    optmodel.constraints.get('maximum_U_constraint').ub = 1
    optmodel.constraints.get('maximum_L_constraint').ub = 1
    MustLU = []
    temp = 0
    while temp < 100:
        optmodel.optimize()
        if optmodel.solver.status == 'optimal':
            if optmodel.solver.objective.value < 0.01:
                break
            variable_values = [v.primal for v in optmodel.variables]
            sol_yL = []
            sol_yU = []
            for j in range(len(variable_names)):
                if variable_names[j].startswith("yL_") and variable_values[j]:
                    if verbose:
                        print('L', variable_names[j][3:])
                    sol_yL.append(variable_names[j][3:])
                if variable_names[j].startswith("yU_") and variable_values[j]:
                    if verbose:
                        print('U', variable_names[j][3:])
                    sol_yU.append(variable_names[j][3:])
            MustLU.append((sol_yL,sol_yU))
            intcut = optmodel.problem.Constraint(add([optmodel.variables.get('yL_'+r) for r in sol_yL]) + 
                                                 add([optmodel.variables.get('yU_'+r) for r in sol_yU]),
                                                 ub=1, name='intcut_LU'+str(temp))
            optmodel.add_cons_vars(intcut)
            temp = temp + 1
        else:
            break
    if high_order_must:
        # MustLU3
        if verbose:
            print('MustLU3')
        optmodel.constraints.get('maximum_U_constraint').ub = 2
        optmodel.constraints.get('maximum_L_constraint').ub = 2
        optmodel.constraints.get('maximum_LU_constraint').ub = 3
        MustLU3 = []
        temp = 0
        while temp < 100:
            optmodel.optimize()
            if optmodel.solver.status == 'optimal':
                if optmodel.solver.objective.value < 0.01:
                    break
                variable_values = [v.primal for v in optmodel.variables]
                sol_yL = []
                sol_yU = []
                for j in range(len(variable_names)):
                    if variable_names[j].startswith("yL_") and variable_values[j]:
                        if verbose:
                            print('L', variable_names[j][3:])
                        sol_yL.append(variable_names[j][3:])
                    if variable_names[j].startswith("yU_") and variable_values[j]:
                        if verbose:
                            print('U', variable_names[j][3:])
                        sol_yU.append(variable_names[j][3:])
                MustLU3.append((sol_yL,sol_yU))
                intcut = optmodel.problem.Constraint(add([optmodel.variables.get('yL_'+r) for r in sol_yL]) + 
                                                    add([optmodel.variables.get('yU_'+r) for r in sol_yU]),
                                                    ub=2, name='intcut_LU3'+str(temp))
                optmodel.add_cons_vars(intcut)
                temp = temp + 1
            else:
                break
        # MustLU4
        if verbose:
            print('MustLU4')
        optmodel.constraints.get('maximum_U_constraint').ub = 2
        optmodel.constraints.get('maximum_L_constraint').ub = 2
        optmodel.constraints.get('maximum_LU_constraint').ub = 4
        MustLU4 = []
        temp = 0
        while temp < 100:
            optmodel.optimize()
            if optmodel.solver.status == 'optimal':
                if optmodel.solver.objective.value < 0.01:
                    break
                variable_values = [v.primal for v in optmodel.variables]
                sol_yL = []
                sol_yU = []
                for j in range(len(variable_names)):
                    if variable_names[j].startswith("yL_") and variable_values[j]:
                        if verbose:
                            print('L', variable_names[j][3:])
                        sol_yL.append(variable_names[j][3:])
                    if variable_names[j].startswith("yU_") and variable_values[j]:
                        if verbose:
                            print('U', variable_names[j][3:])
                        sol_yU.append(variable_names[j][3:])
                MustLU4.append((sol_yL,sol_yU))
                intcut = optmodel.problem.Constraint(add([optmodel.variables.get('yL_'+r) for r in sol_yL]) + 
                                                    add([optmodel.variables.get('yU_'+r) for r in sol_yU]),
                                                    ub=3, name='intcut_LU4'+str(temp))
                optmodel.add_cons_vars(intcut)
                temp = temp + 1
            else:
                break
        solution = {'MustL': MustL, 'MustU': MustU, 'MustLL': MustLL, 'MustUU': MustUU, 'MustLU': MustLU,
                    'MustLU3': MustLU3, 'MustLU4': MustLU4}
        MustL_set = set(MustL + sum(MustLL,[]) + sum([x[0] for x in MustLU],[]) + 
                        sum([x[0] for x in MustLU3],[]) + sum([x[0] for x in MustLU4],[]))
        MustU_set = set(MustU + sum(MustUU,[]) + sum([x[1] for x in MustLU],[]) + 
                        sum([x[1] for x in MustLU3],[]) + sum([x[1] for x in MustLU4],[]))
    else:
        solution = {'MustL': MustL, 'MustU': MustU, 'MustLL': MustLL, 'MustUU': MustUU, 'MustLU': MustLU}
        MustL_set = set(MustL + sum(MustLL,[]) + sum([x[0] for x in MustLU],[]))
        MustU_set = set(MustU + sum(MustUU,[]) + sum([x[1] for x in MustLU],[]))
    return solution, MustL_set, MustU_set

def add_must(model, biomass, target, minimum_growth, minimum_target, fva_solution_WT, maximum_L, maximum_U,
             excluded_reactions=[], primal_bound=True):
    """
    Add variables, constraints, and objective for OptForce

    Parameters
    ----------
    model : cobra.Model
        The model to add OptForce constraints and objective to
    biomass : str
        Biomass reaction id
    target : str
        Target reaction id
    minimum_growth : float
        Minimum growth requirement
    minimum_target : float
        Minimum target production
    fva_solution_WT : ouput from cobra.flux_analysis.flux_variability_analysis
        Minimum and maximum possible flux value for each reaction in WT
    maximum_L : int
        Maximum number of reactions with decreased flux
    maximum_U : int
        Maximum number of reactions with increased flux
    exlcluded_reactions : list, optional
        List of reactions excluded from must set (default [])
    primal_bound : bool, optional
        Whether to place [-1000.0, 1000.0] bounds on the primal variables (default True)
    """
    # Get the solver name
    prob = model.problem
    solver_str = interface_to_str(prob)
    model.tolerance = 1e-9
    # Remove arbitrary bounds (usually -1000.0 and 1000.0)
    if not primal_bound:
        for r in model.reactions:
            if r.lower_bound <= -1000.0:
                r.lower_bound = -np.inf
            if r.upper_bound >= 1000.0:
                r.upper_bound = np.inf
    model.reactions.get_by_id(biomass).lower_bound = minimum_growth
    model.reactions.get_by_id(target).lower_bound = minimum_target
    v_min_WT = fva_solution_WT['minimum'].to_dict()
    v_max_WT = fva_solution_WT['maximum'].to_dict()
    subset_reactions = [r.id for r in model.reactions if r.id not in excluded_reactions]
    # Define variables
    var_mu = dict()
    var_l = dict()
    var_u = dict()
    var_wL = dict()
    var_yL = dict()
    var_wU = dict()
    var_yU = dict()
    for m in model.metabolites:
        var_mu[m.id] = prob.Variable('mu_'+m.id)
    for r in model.reactions:
        if r.lower_bound > -np.inf:
            var_l[r.id] = prob.Variable('l_'+r.id, lb=0)
        if r.upper_bound < np.inf:
            var_u[r.id] = prob.Variable('u_'+r.id, lb=0)
        if r.id in subset_reactions:
            var_wL[r.id] = prob.Variable('wL_'+r.id)
            var_yL[r.id] = prob.Variable('yL_'+r.id, type='binary')
            var_wU[r.id] = prob.Variable('wU_'+r.id)
            var_yU[r.id] = prob.Variable('yU_'+r.id, type='binary')
    # Add constraints
    maximum_L_constraint = prob.Constraint(add([var_yL[r.id] for r in model.reactions if r.id in subset_reactions]),
                                           ub=maximum_L, name='maximum_L_constraint')
    maximum_U_constraint = prob.Constraint(add([var_yU[r.id] for r in model.reactions if r.id in subset_reactions]),
                                           ub=maximum_U, name='maximum_U_constraint')
    maximum_LU_constraint = prob.Constraint(add([var_yL[r.id] for r in model.reactions if r.id in subset_reactions]) +
                                            add([var_yU[r.id] for r in model.reactions if r.id in subset_reactions]),
                                            ub=maximum_L+maximum_U, name='maximum_LU_constraint')
    strong_duality = prob.Constraint(add([var_wL[r.id] for r in model.reactions if r.id in subset_reactions]) -
                                     add([var_wU[r.id] for r in model.reactions if r.id in subset_reactions]) +
                                     add([r.lower_bound*var_l[r.id] for r in model.reactions if r.lower_bound > -np.inf]) -
                                     add([r.upper_bound*var_u[r.id] for r in model.reactions if r.upper_bound < np.inf]),
                                     lb=0, ub=0)
    dual_constraint = dict()
    for r in model.reactions:
        dual_constraint[r.id] = prob.Constraint(add([s*var_mu[m.id] for m, s in r.metabolites.items()]) -
                                                (var_l[r.id] if r.lower_bound > -np.inf else 0) +
                                                (var_u[r.id] if r.upper_bound < np.inf else 0) -
                                                (var_yL[r.id] if r.id in subset_reactions else 0) +
                                                (var_yU[r.id] if r.id in subset_reactions else 0),
                                                lb=0, ub=0)
    model.add_cons_vars([maximum_L_constraint, maximum_U_constraint, maximum_LU_constraint, 
                         strong_duality, *dual_constraint.values()])
    # optlang gurobi interface does not support indicator
    if solver_str == 'gurobi':
        model.solver.update()
        for r in model.reactions:
            if r.id in subset_reactions:
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yL[r.id].name), True,
                                                            model.solver.problem.getVarByName(var_wL[r.id].name) -
                                                            model.solver.problem.getVarByName(r.forward_variable.name) + 
                                                            model.solver.problem.getVarByName(r.reverse_variable.name) == 0)
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yL[r.id].name), False,
                                                            model.solver.problem.getVarByName(var_wL[r.id].name) == 0)
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yU[r.id].name), True,
                                                            model.solver.problem.getVarByName(var_wU[r.id].name) -
                                                            model.solver.problem.getVarByName(r.forward_variable.name) + 
                                                            model.solver.problem.getVarByName(r.reverse_variable.name) == 0)
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yU[r.id].name), False,
                                                            model.solver.problem.getVarByName(var_wU[r.id].name) == 0)
    else:
        indicator_constraintL1 = dict()
        indicator_constraintL2 = dict()
        indicator_constraintU1 = dict()
        indicator_constraintU2 = dict()
        for r in model.reactions:
            if r.id in subset_reactions:
                indicator_constraintL1[r.id] = prob.Constraint(var_wL[r.id] - r.flux_expression, lb=0, ub=0,
                                                                indicator_variable=var_yL[r.id], active_when=1)
                indicator_constraintL2[r.id] = prob.Constraint(var_wL[r.id], lb=0, ub=0,
                                                                indicator_variable=var_yL[r.id], active_when=0)
                indicator_constraintU1[r.id] = prob.Constraint(var_wU[r.id] - r.flux_expression, lb=0, ub=0,
                                                                indicator_variable=var_yU[r.id], active_when=1)
                indicator_constraintU2[r.id] = prob.Constraint(var_wU[r.id], lb=0, ub=0,
                                                                indicator_variable=var_yU[r.id], active_when=0)
        model.add_cons_vars([*indicator_constraintL1.values(), *indicator_constraintL2.values(),
                             *indicator_constraintU1.values(), *indicator_constraintU2.values()])
    indicator_constraintLU = dict()
    for r in model.reactions:
        if r.id in subset_reactions:
            indicator_constraintLU[r.id] = prob.Constraint(var_yL[r.id] + var_yU[r.id], ub=1)
    model.add_cons_vars([*indicator_constraintLU.values()])
    # Set the outer problem objective
    model.objective = prob.Objective(add([v_min_WT[r]*var_yL[r] for r in subset_reactions] +
                                         [-var_wL[r] for r in subset_reactions]) +
                                     add([-v_max_WT[r]*var_yU[r] for r in subset_reactions] + 
                                         [var_wU[r] for r in subset_reactions]))
    # Set solver specific parameters
    if solver_str == 'gurobi':
        model.solver.problem.params.FeasibilityTol = 1e-9
        model.solver.problem.params.OptimalityTol = 1e-9
        model.solver.problem.params.IntFeasTol = 1e-9
        model.solver.problem.params.MIPgapAbs = 0
        model.solver.problem.params.MIPgap = 0
    elif solver_str == 'cplex':
        model.solver.problem.parameters.simplex.tolerances.feasibility.set(1e-9)
        model.solver.problem.parameters.simplex.tolerances.optimality.set(1e-9)
        model.solver.problem.parameters.mip.tolerances.integrality.set(0)
        model.solver.problem.parameters.mip.tolerances.absmipgap.set(0)
        model.solver.problem.parameters.mip.tolerances.mipgap.set(0)
    model.solver.update()

def find_force(model, biomass, target, maximum_knockout, minimum_growth,
               fva_solution_MT, MustL_set, MustU_set, maximum_L, maximum_U,
               excluded_reactions=[], knockout_penalty=0.0,
               number_of_solutions=10, time_limit=60, threads=None, verbose=False):
    """
    Identify force set
   
    Parameters
    ----------
    model : cobra.Model
        The model to compute an OptForce solution for
    biomass : str
        Biomass reaction id
    target : str
        Target reaction id
    maximum_knockout : int
        Maximum knockout limit
    minimum_growth : float
        Minimum growth requirement
    fva_solution_MT : ouput from cobra.flux_analysis.flux_variability_analysis
        Minimum and maximum possible flux value for each reaction in MT
    MustL_set : set
        A set of mustL reactions
    MustU_set : set
        A set of mustU reactions
    maximum_L : int
        Maximum number of reactions with decreased flux
    maximum_U : int
        Maximum number of reactions with increased flux
    exlcluded_reactions : list, optional
        List of reactions excluded from must set (default [])
    knockout_penalty : float, optional
        Minimum increase in target production per knockout (default 0.0)
    number_of_solutions : int, optional
        Maximum number of OptForce solutions (default 10)
    time_limit : any nonnegative value, optional
        Maximum time in seconds for each optimization (default 60)
    threads : int, optional
        Maximum number of parallel threads (default None)
    verbose : bool, optional
        Whether to print the solver log (default False)

    Returns
    -------
    solution
        A dictionary containing the force set solution
    See Also
    --------
    add_force : add variables, constraints, and objective for force set identification
    """
    # Make a copy of the model and add variables, constraints, and objective 
    optmodel = model.copy()
    solver_str = interface_to_str(optmodel.problem)
    add_force(optmodel, biomass, target, maximum_knockout, minimum_growth,
              fva_solution_MT, MustL_set, MustU_set, maximum_L, maximum_U,
              excluded_reactions, knockout_penalty)
    # Set solver specific parameters
    if solver_str == 'gurobi':
        if time_limit:
            optmodel.solver.problem.params.TimeLimit = time_limit
        if threads:
            optmodel.solver.problem.params.Threads = threads
        if verbose:
            optmodel.solver.problem.params.OutputFlag = 1
    elif solver_str == 'cplex':
        if time_limit:
            optmodel.solver.problem.parameters.timelimit.set(time_limit)
        if threads:
            optmodel.solver.problem.parameters.threads.set(threads)
        if verbose:
            optmodel.solver.problem.set_results_stream(sys.stdout)
    # Solve OptForce and read solution
    # optlang's optimize function involves extra steps
    # Gurobi interface resets if the existing status is not optimal
    solution = []
    temp = 0
    while temp < number_of_solutions:
        if solver_str == 'gurobi':
            optmodel.solver.problem.optimize()
            status = optmodel.solver.problem.getAttr("Status")
            objval = optmodel.solver.problem.getAttr("ObjVal")
            if verbose:
                print(temp, 'status:', status, 'objval:', objval)
            if (status == 2 or status == 9) and objval > 0.01:
                variable_names = [v.VarName for v in optmodel.solver.problem.getVars()]
                variable_values = [v.x for v in optmodel.solver.problem.getVars()]
            else:
                break
        elif solver_str == 'cplex':
            optmodel.solver.problem.solve()
            status = optmodel.solver.problem.solution.get_status()
            objval = optmodel.solver.problem.solution.get_objective_value()
            if verbose:
                print(temp, 'status:', status, 'objval:', objval)
            if (status == 101 or status == 107) and objval > 0.01:
                variable_names = optmodel.solver.problem.variables.get_names()
                variable_values = optmodel.solver.problem.solution.get_values()
            else:
                break
        sol_yU = []
        sol_yL = []
        sol_yK = []
        for j in range(len(variable_names)):
            if variable_names[j].startswith("yU_") and variable_values[j]:
                if verbose:
                    print('U', variable_names[j][3:], variable_values[j])
                sol_yU.append(variable_names[j][3:])
            if variable_names[j].startswith("yL_") and variable_values[j]:
                if verbose:
                    print('L', variable_names[j][3:], variable_values[j])
                sol_yL.append(variable_names[j][3:])
            if variable_names[j].startswith("yK_") and variable_values[j]:
                if verbose:
                    print('K', variable_names[j][3:], variable_values[j])
                sol_yK.append(variable_names[j][3:])
        solution.append({'status': status, 'objective value': objval, 'solution': (sol_yU,sol_yL,sol_yK)})
        intcut = optmodel.problem.Constraint(add([optmodel.variables.get('yU_'+r) for r in sol_yU]) +
                                             add([optmodel.variables.get('yL_'+r) for r in sol_yL]) +
                                             add([optmodel.variables.get('yK_'+r) for r in sol_yK]),
                                             ub=len(sol_yU+sol_yL+sol_yK)-1, name='intcut'+str(temp))
        optmodel.add_cons_vars(intcut)
        optmodel.solver.update()
        temp = temp + 1
    return solution
    
def add_force(model, biomass, target, maximum_knockout, minimum_growth,
              fva_solution_MT, MustL_set, MustU_set, maximum_L, maximum_U, 
              excluded_reactions=[], knockout_penalty=0.0, primal_bound=True):
    """
    Add variables, constraints, and objective for OptForce

    Parameters
    ----------
    model : cobra.Model
        The model to add OptForce constraints and objective to
    biomass : str
        Biomass reaction id
    target : str
        Target reaction id
    maximum_knockout : int
        Maximum knockout limit
    minimum_growth : float
        Minimum growth requirement
    fva_solution_MT : ouput from cobra.flux_analysis.flux_variability_analysis
        Minimum and maximum possible flux value for each reaction in MT
    MustL_set : set
        A set of mustL reactions
    MustU_set : set
        A set of mustU reactions
    maximum_L : int
        Maximum number of reactions with decreased flux
    maximum_U : int
        Maximum number of reactions with increased flux
    exlcluded_reactions : list, optional
        List of reactions excluded from must set (default [])
    knockout_penalty : float, optional
        Minimum increase in target production per knockout (default 0.0)
    primal_bound : bool, optional
        Whether to place [-1000.0, 1000.0] bounds on the primal variables (default True)
    """
    # Get the solver name
    prob = model.problem
    solver_str = interface_to_str(prob)
    model.tolerance = 1e-9
    # Remove arbitrary bounds (usually -1000.0 and 1000.0)
    if not primal_bound:
        for r in model.reactions:
            if r.lower_bound <= -1000.0:
                r.lower_bound = -np.inf
            if r.upper_bound >= 1000.0:
                r.upper_bound = np.inf
    model.reactions.get_by_id(biomass).lower_bound = minimum_growth
    v_min_MT = fva_solution_MT['minimum'].to_dict()
    v_max_MT = fva_solution_MT['maximum'].to_dict()
    subset_reactions = [r.id for r in model.reactions if r.id not in excluded_reactions]
    # Define variables
    var_mu = dict()
    var_l = dict()
    var_u = dict()
    var_uL = dict()
    var_yL = dict()
    var_lU = dict()
    var_yU = dict()
    var_k = dict()
    var_yK = dict()
    for m in model.metabolites:
        var_mu[m.id] = prob.Variable('mu_'+m.id)
    for r in model.reactions:
        if r.lower_bound > -np.inf:
            var_l[r.id] = prob.Variable('l_'+r.id, lb=0)
        if r.upper_bound < np.inf:
            var_u[r.id] = prob.Variable('u_'+r.id, lb=0)
        if r.id in MustL_set:
            var_uL[r.id] = prob.Variable('uL_'+r.id, lb=0)
            var_yL[r.id] = prob.Variable('yL_'+r.id, type='binary')
        if r.id in MustU_set:
            var_lU[r.id] = prob.Variable('lU_'+r.id, lb=0)
            var_yU[r.id] = prob.Variable('yU_'+r.id, type='binary')
        if r.id in subset_reactions:
            var_k[r.id] = prob.Variable('k_'+r.id)
            var_yK[r.id] = prob.Variable('yK_'+r.id, type='binary')
        model.solver.update()
    # Add constraints
    maximum_L_constraint = prob.Constraint(add([var_yL[r.id] for r in model.reactions if r.id in MustL_set]),
                                           ub=maximum_L, name='maximum_L_constraint')
    maximum_U_constraint = prob.Constraint(add([var_yU[r.id] for r in model.reactions if r.id in MustU_set]),
                                           ub=maximum_U, name='maximum_U_constraint')
    
    print(f'Inside add_force function. The value of maximum_knockout is {maximum_knockout}')
    maximum_K_constraint = prob.Constraint(add([var_yK[r.id] for r in model.reactions if r.id in subset_reactions]),
                                           ub=maximum_knockout, name='maximum_K_constraint')
    strong_duality = prob.Constraint(add([r.flux_expression for r in model.reactions
                                         if r.id == target]) -
                                     add([r.lower_bound*var_l[r.id] for r in model.reactions 
                                         if r.lower_bound > -np.inf]) +
                                     add([r.upper_bound*var_u[r.id] for r in model.reactions
                                         if r.upper_bound < np.inf]) -
                                     add([v_min_MT[r.id]*var_lU[r.id] for r in model.reactions 
                                         if r.id in MustU_set]) +
                                     add([v_max_MT[r.id]*var_uL[r.id] for r in model.reactions 
                                         if r.id in MustL_set]),
                                     lb=0, ub=0)
    dual_constraint = dict()
    for r in model.reactions:
        dual_constraint[r.id] = prob.Constraint(add([s*var_mu[m.id] for m, s in r.metabolites.items()]) +
                                                (var_l[r.id] if r.lower_bound > -np.inf else 0) -
                                                (var_u[r.id] if r.upper_bound < np.inf else 0) +
                                                (var_lU[r.id] if r.id in MustU_set else 0) -
                                                (var_uL[r.id] if r.id in MustL_set else 0) +
                                                (var_k[r.id] if r.id in subset_reactions else 0) -
                                                (1 if r.id == target else 0),
                                                lb=0, ub=0)
    model.add_cons_vars([maximum_L_constraint, maximum_U_constraint, maximum_K_constraint, 
                         strong_duality, *dual_constraint.values()])
    # optlang gurobi interface does not support indicator
    if solver_str == 'gurobi':
        model.solver.update()
        for r in model.reactions:
            if r.id in MustL_set:
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yL[r.id].name), True,
                                                           model.solver.problem.getVarByName(r.forward_variable.name) - 
                                                           model.solver.problem.getVarByName(r.reverse_variable.name) <= v_max_MT[r.id])
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yL[r.id].name), False,
                                                           model.solver.problem.getVarByName(var_uL[r.id].name) == 0)
                if r.upper_bound < np.inf:
                    model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yL[r.id].name), True,
                                                               model.solver.problem.getVarByName(var_u[r.id].name) == 0)
            if r.id in MustU_set:
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yU[r.id].name), True,
                                                           model.solver.problem.getVarByName(r.forward_variable.name) - 
                                                           model.solver.problem.getVarByName(r.reverse_variable.name) >= v_min_MT[r.id])
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yU[r.id].name), False,
                                                           model.solver.problem.getVarByName(var_lU[r.id].name) == 0)
                if r.lower_bound > -np.inf:
                    model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yU[r.id].name), True,
                                                               model.solver.problem.getVarByName(var_l[r.id].name) == 0)
            if r.id in subset_reactions:
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yK[r.id].name), True,
                                                           model.solver.problem.getVarByName(r.forward_variable.name) - 
                                                           model.solver.problem.getVarByName(r.reverse_variable.name) == 0)
                model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yK[r.id].name), False,
                                                           model.solver.problem.getVarByName(var_k[r.id].name) == 0)
                if r.lower_bound > -np.inf:
                    model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yK[r.id].name), True,
                                                               model.solver.problem.getVarByName(var_l[r.id].name) == 0)
                if r.upper_bound < np.inf:
                    model.solver.problem.addGenConstrIndicator(model.solver.problem.getVarByName(var_yK[r.id].name), True,
                                                               model.solver.problem.getVarByName(var_u[r.id].name) == 0)
    else:
        indicator_constraintL1 = dict()
        indicator_constraintL2 = dict()
        indicator_constraintL3 = dict()
        indicator_constraintU1 = dict()
        indicator_constraintU2 = dict()
        indicator_constraintU3 = dict()
        indicator_constraintK1 = dict()
        indicator_constraintK2 = dict()
        indicator_constraintK3 = dict()
        indicator_constraintK4 = dict()
        for r in model.reactions:
            if r.id in MustL_set:
                indicator_constraintL1[r.id] = prob.Constraint(r.flux_expression, ub=v_max_MT[r.id],
                                                               indicator_variable=var_yL[r.id], active_when=1)
                indicator_constraintL2[r.id] = prob.Constraint(var_uL[r.id], lb=0, ub=0,
                                                               indicator_variable=var_yL[r.id], active_when=0)
                if r.upper_bound < np.inf:
                    indicator_constraintL3[r.id] = prob.Constraint(var_u[r.id], lb=0, ub=0,
                                                                   indicator_variable=var_yL[r.id], active_when=1)
            if r.id in MustU_set:
                indicator_constraintU1[r.id] = prob.Constraint(r.flux_expression, lb=v_min_MT[r.id],
                                                               indicator_variable=var_yU[r.id], active_when=1)
                indicator_constraintU2[r.id] = prob.Constraint(var_lU[r.id], lb=0, ub=0,
                                                               indicator_variable=var_yU[r.id], active_when=0)
                if r.lower_bound > -np.inf:
                    indicator_constraintU3[r.id] = prob.Constraint(var_l[r.id], lb=0, ub=0,
                                                                   indicator_variable=var_yU[r.id], active_when=1)
            if r.id in subset_reactions:
                indicator_constraintK1[r.id] = prob.Constraint(r.flux_expression, lb=0, ub=0,
                                                               indicator_variable=var_yK[r.id], active_when=1)
                indicator_constraintK2[r.id] = prob.Constraint(var_k[r.id], lb=0, ub=0,
                                                               indicator_variable=var_yK[r.id], active_when=0)
                if r.lower_bound > -np.inf:
                    indicator_constraintK3[r.id] = prob.Constraint(var_l[r.id], lb=0, ub=0,
                                                                   indicator_variable=var_yK[r.id], active_when=1)
                if r.upper_bound < np.inf:
                    indicator_constraintK4[r.id] = prob.Constraint(var_u[r.id], lb=0, ub=0,
                                                                   indicator_variable=var_yK[r.id], active_when=1)
        model.add_cons_vars([*indicator_constraintL1.values(), *indicator_constraintL2.values(), *indicator_constraintL3.values(),
                             *indicator_constraintU1.values(), *indicator_constraintU2.values(), *indicator_constraintU3.values(),
                             *indicator_constraintK1.values(), *indicator_constraintK2.values(),
                             *indicator_constraintK3.values(), *indicator_constraintK4.values()])
    # Set the outer problem objective
    model.objective = prob.Objective(add([r.flux_expression for r in model.reactions if r.id == target]) -
                                     knockout_penalty*(add([var_yL[r] for r in MustL_set] +
                                                           [var_yU[r] for r in MustU_set] +
                                                           [var_yK[r] for r in subset_reactions])))
    # The objective offset is commented out in the optlang interface
    if solver_str == 'gurobi':
        model.solver.problem.params.FeasibilityTol = 1e-9
        model.solver.problem.params.OptimalityTol = 1e-9
        model.solver.problem.params.IntFeasTol = 1e-9
        model.solver.problem.params.MIPgapAbs = 0
        model.solver.problem.params.MIPgap = 0
    elif solver_str == 'cplex':
        model.solver.problem.parameters.simplex.tolerances.feasibility.set(1e-9)
        model.solver.problem.parameters.simplex.tolerances.optimality.set(1e-9)
        model.solver.problem.parameters.mip.tolerances.integrality.set(0)
        model.solver.problem.parameters.mip.tolerances.absmipgap.set(0)
        model.solver.problem.parameters.mip.tolerances.mipgap.set(0)
    model.solver.update()