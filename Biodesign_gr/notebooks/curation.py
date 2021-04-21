# This function uses the fuctions found in notebooks/curation_gr.ipynb
# That file is signicantly more reader friendly

from matplotlib import pyplot as plt
from IPython.display import IFrame
import numpy as np
import pandas as pd
import json
import urllib
import cobra
import cplex
import os
import requests
import collections
import itertools

# Main Function First functions first

def curate(model):
    
    # Check model for metabolites with equivalent formulas
    for m in [m for m in model.metabolites if ';' in m.formula and not all_formulas_equivalent(m)]:
        m.formula = m.formula.split(';')[0]
        
    # Assign metabolite formulas by checking if one formula makes all reactions the where its the only undefined metabolite balanced
    metabolites_that_can_be_defined = 1
    while metabolites_that_can_be_defined > 0:
        metabolites_that_can_be_defined = 0
        for m in [m for m in model.metabolites if ';' in m.formula]:
            for f in m.formula.split(';'):
                if fraction_of_reactions_formula_balances(m, f, reactions_where_m_is_only_undefined_metabolite(m)) == 1:
                    metabolites_that_can_be_defined += 1
                    m.formula = f
                    
    # Repeat but allow for imperfect fitting
    highest_fraction = 1
    while highest_fraction > 0:
        highest_fraction = 0
        metabolites_that_can_be_defined = 0

        # find highest fraction of reactions that are solved by a given formula
        for m in [m for m in model.metabolites if ';' in m.formula]:
            for f in m.formula.split(';'):
                if fraction_of_reactions_formula_balances(m, f, reactions_where_m_is_only_undefined_metabolite(m)) > highest_fraction:
                    highest_fraction = fraction_of_reactions_formula_balances(m, f, reactions_where_m_is_only_undefined_metabolite(m))

        # assign formulas to metabolites with formula that gives a score equal to the best fraction
        if highest_fraction > 0:
            for m in [m for m in model.metabolites if ';' in m.formula]:
                for f in m.formula.split(';'):
                    if fraction_of_reactions_formula_balances(m, f, reactions_where_m_is_only_undefined_metabolite(m)) == highest_fraction:
                        m.formula = f
                        metabolites_that_can_be_defined += 1
                        
    # Assign remaining metabolites by checking if combinations of formulas are balanaced or only off by hydrogen
    undefined_metabolites = [m for m in model.metabolites if ';' in m.formula]
    possible_formulas = [m.formula.split(';') for m in model.metabolites if ';' in m.formula]

    # only do this step if there is a reasonable number of undefined metbolites due to exponential growth of formula combinations
    if len(undefined_metabolites) < 10:

        # inital best formulas is their original values
        best_formulas = [m.formula for m in undefined_metabolites]
        best_score = len(reactions_off_by_more_than_hydrogen(model))
        best_length = 0

        # goes through all permutations of formulas for undefined metabolites
        for formulas in list(itertools.product(*possible_formulas)):
            # assign the formulas to the metabolites
            for count, m in enumerate(undefined_metabolites):
                model.metabolites.get_by_id(m.id).formula = formulas[count]

            # get the number of reactions that are off by more than hydrogen
            unacceptable_reactions = reactions_off_by_more_than_hydrogen(model)

            # if its the best fit replace the best formulas
            if len(unacceptable_reactions) <= best_score:
                if chars_in_string_list(formulas) > best_length:
                    best_formulas = formulas
                    best_score = len(reactions_off_by_more_than_hydrogen(model))
                    best_length = chars_in_string_list(formulas)
                    
    for count, m in enumerate(undefined_metabolites):
        model.metabolites.get_by_id(m.id).formula = best_formulas[count]
              
    # Balance Hydrogen
    for r in [r for r in model.reactions if only_hydrogen_unbalanced(r)]:
        fix_unbalanced_hydrogen(model, r)
        
    return model

        ####################
        #                  #
        # Helper Functions #
        #                  #
        ####################
    

def get_initial_number_string(substring):
    initial_string = ''
    for char in substring:
        if char.isdigit():
            initial_string += char
        else:
            return initial_string
    return initial_string

def formula_dict_from_string(formula_string):
    formula_dict = {}
    elements = [char for char in formula_string if char.isalpha()]
    for element in elements:
        string_after_element = formula_string.split(element, 1)[1]
        coefficient = get_initial_number_string(string_after_element)
        if coefficient == '':
            coefficient = '1'
        formula_dict[element] = int(coefficient)
    return formula_dict

def all_formulas_equivalent(m):
    first_formula = m.formula.split(';')[0]
    return len([f for f in m.formula.split(';') if formula_dict_from_string(f) != formula_dict_from_string(first_formula)]) > 0

def m_only_undefined_metabolite(m1, r):
    return ';' in m1.formula and len([m2 for m2 in r.metabolites if ';' in m2.formula and m1 != m2]) == 0

def reactions_where_m_is_only_undefined_metabolite(m):
    return [r for r in m.reactions if m_only_undefined_metabolite(m,r)]

def fraction_of_reactions_formula_balances(m, formula, rxn_list):
    original_formula = m.formula
    m.formula = formula
    balanced_reactions   = [r for r in rxn_list if r.check_mass_balance() == {}]
    unbalanced_reactions = [r for r in rxn_list if r.check_mass_balance() != {}]
    m.formula = original_formula
    
    # avoid divide by zero
    if len(balanced_reactions) + len(unbalanced_reactions) == 0:
        return 0
    return len(balanced_reactions) / (len(balanced_reactions) + len(unbalanced_reactions))            
                    
def is_balanced(r):
    return abs(sum(list(r.check_mass_balance().values()))) < 1e-5

# This one needs to be fixed
def should_be_balanced(r):
    return not (r.id.startswith('EX_') or r.id.startswith('sink_') or r.id.startswith('Growth'))

def balanced_or_only_hydrogen_unbalanced(r):
    return should_be_balanced(r) and (is_balanced(r) or list(r.check_mass_balance().keys()) == ['H'])

def reactions_off_by_more_than_hydrogen(model):
    return [r for r in model.reactions if should_be_balanced(r) and not balanced_or_only_hydrogen_unbalanced(r)]

def chars_in_string_list(string_list):
    total_chars = 0
    for string in string_list:
        total_chars += len(string)
    return total_chars
    
def only_hydrogen_unbalanced(r):
    return list(r.check_mass_balance().keys()) == ['H'] and should_be_balanced(r)

def fix_unbalanced_hydrogen(model, r):
    hydrogen_error = int(r.check_mass_balance()['H'])
    r.subtract_metabolites({model.metabolites.get_by_id("h_c"): hydrogen_error})