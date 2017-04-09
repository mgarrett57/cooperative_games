# -*- coding: utf-8 -*-
"""
coop_games.py

@author: mgarrett
"""
from itertools import combinations
from scipy.special import comb

SOLUTION_FUNCTIONS = ('Shapley','Banzhaf')
        
class CooperativeGame:
    '''
    This class finds solutions for fully connected cooperative games,
    i.e. undirected graph with an edge connecting each player.
    
    Parameters
    ----------
    inputted_player_set : list-like or set-like
        Total players in game. Player set gets converted into a 'frozenset' data type.
        
    cost_function : dictionary with tuple/frozenset keys
        Tuple keys are required for directed graphs
        
    minimize : boolean
        True if cost minimization game, false if gain maximization game.
        Default is True.

    Attribute
    ----------
    player_set : list-like or set-like
        Specified player set
    
    cost_function : dictionary with tuple/frozenset keys
        Specified cost function, with minor adjustments (null coalition be set
        to zero, if not explicitly submitted)
        
    minimize : boolean
        Specified 'minimize' boolean
    
    a_priori_unions : list-like
        A list of a priori unions that are used for calculating the Owen value
    
    solutions : dictionary
        Dictionary of all the solutions requested. Example: after running 
        'calculate_ shapley' method, Shapley values for 'Player1' can be found
        in self.solutions['Shapley']['Player1'] and all Shapley values are in 
        self.solutions['Shapley']
    '''
    
    def __init__(self, inputted_player_set, cost_function, minimize = True):
        self.player_set = frozenset(inputted_player_set)
        self.n_players = len(inputted_player_set)
        self.minimize = minimize
        self._validate_cost_function(cost_function)

        #create empty dictionary for solutions for each player
        self.solutions = {sltns:{plyrs:None for plyrs in self.player_set}
                          for sltns in SOLUTION_FUNCTIONS}
        
    def _validate_cost_function(self, cost_function):
        n = self.n_players
        
        #verify size of cost_function is correct
        expected_size = 2.**n
        if len(cost_function) == expected_size-1:
            assert () not in cost_function, "Cost function is not a complete graph"
            cost_function[()]=0
        else:
            assert len(cost_function) == expected_size, "Cost function is not a complete graph"
            
        #convert coalitions to immutable frozensets
        cost_function = {frozenset(coalition):value for coalition,value in cost_function.items()}        
        self.cost_function = cost_function
        return self
    
    def get_grand_coalition_value(self):
        return self.cost_function[self.player_set]
        
    def normalize_solutions(self, solutions = ['Shapley','Banzhaf']):
        try:
            banzhaf_total = sum(x for x in self.solutions['Banzhaf'].itervalues()
                                if x is not None)
        except AttributeError:
            banzhaf_total = sum(x for x in self.solutions['Banzhaf'].values()
                                if x is not None)
        for p1 in self.player_set:
            if self.solutions['Banzhaf'][p1] is not None:
                self.solutions['Banzhaf'][p1] /= banzhaf_total
            for sol_type in ['Shapley']:
                if self.solutions[sol_type][p1] is not None:
                    self.solutions[sol_type][p1] /= self.get_grand_coalition_value()

    def marginal_contribution(self,player,combination):
        if player in combination:
            sans_player = frozenset({x for x in combination if x != player})
            if self.minimize:
                marginal = self.cost_function[combination] - self.cost_function[sans_player]
            else:
                marginal = self.cost_function[sans_player] - self.cost_function[combination]
        else:
            marginal = 0
        return marginal

    def Shapley(self):
        #calculate coefficents for marginal contribs: [(S-1)!(N-S)!]/N!
        n = self.n_players
        shapley_coefficients = [1/(comb(n,s)*s) for s in range(1,n+1)]
        
        #initialize Shapley value to 0 for each player           
        for p0 in self.player_set:
            self.solutions['Shapley'][p0] = 0
            
        #loop through different coalition sizes to add marginal contributions for each player
        for i in range(self.n_players):
            set_i = [frozenset(j) for j in combinations(self.player_set, i+1)]
            for p1 in self.player_set:
                temp_contribution = 0
                n_combos = 0
                for combo in set_i:
                    if p1 in combo:
                        coalition_size = len(combo)
                        n_combos += 1
                        temp_marginal = self.marginal_contribution(p1,combo)
                        temp_contribution += shapley_coefficients[coalition_size-1]*temp_marginal
                self.solutions['Shapley'][p1] += temp_contribution
        return self.solutions['Shapley']
 
    def Banzhaf(self, standardize = False):
        #http://www.lamsade.dauphine.fr/~airiau/Teaching/CoopGames/2011/lecture8.pdf
        #calculate the banzhaf coefficient
        n = self.n_players
        banzhaf_coefficient = 1.0/(2**(n-1))
        
        #initialize Banzhaf index value to 0 for each player           
        for p0 in self.player_set:
            self.solutions['Banzhaf'][p0] = 0
            
        #loop through different coalition sizes to add marginal contributions for each player
        total_contribution = 0
        for i in range(self.n_players):
            set_i = [frozenset(j) for j in combinations(self.player_set, i+1)]
            for p1 in self.player_set:
                temp_contribution = 0
                for combo in set_i:
                    if p1 in combo:
                        temp_marginal = self.marginal_contribution(p1,combo)
                        temp_contribution += temp_marginal
                self.solutions['Banzhaf'][p1] += temp_contribution
        for p2 in self.player_set:
            self.solutions['Banzhaf'][p2] *= banzhaf_coefficient
            
        #if standardization is requested
        if standardize:    
            try:
                total_contribution = sum(x for x in self.solutions['Banzhaf'].itervalues())
            except AttributeError:
                total_contribution = sum(x for x in self.solutions['Banzhaf'].values())
            standardization_factor = self.get_grand_coalition_value()/total_contribution
            for p3 in self.player_set:
                self.solutions['Banzhaf'][p3] *= standardization_factor                      
        return self.solutions['Banzhaf']
    
