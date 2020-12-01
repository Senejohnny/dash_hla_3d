""" This file contains all the methods consumed by Epitope Data Base """

import pandas as pd
from collections import defaultdict
from app.analysis_utils import flatten2set

class Epitope_DB:
    
    def __init__(self, path:str='../data/20201123_EpitopevsHLA.pickle'):
        self.df = pd.read_pickle(path)
        
    def __repr__(self):
        return f""" Epitope_DB(records={len(self.df)}, columns={self.df.columns}) """
    
    def __str__(self):
        return __repr__()
    
    def epitope(self, value):
        if isinstance(value, str):
            ind = self.df.Epitope == value 
        else: 
            ind = self.df.Epitope.apply(lambda x: x in value) 
        self.df = self.df[ind]
        return self
    
    def ellipro(self, value):
        if isinstance(value, str):
            ind = self.df.Epitope == value 
        else: 
            ind = self.df['ElliPro Score'].apply(lambda x: x in value) 
        self.df = self.df[ind]
        return self
    
    def hlavsep(self, hla_allel:str='Luminex Alleles'):
        
        hlas = flatten2set(self.df[hla_allel].values)
        hlavsep_dict = defaultdict(list)
        for hla in hlas:
            ind = self.df[hla_allel].apply(lambda x: hla in x)
            epitopes = flatten2set(self.df[ind]['Epitope'].values)
            hlavsep_dict['HLA'].append(hla)
            hlavsep_dict['Epitope'].append(epitopes)
        self._hlavsep = pd.DataFrame(hlavsep_dict)
        return self._hlavsep

    def min_hlavsep(self, epitopes):
        _epitopes = epitopes.copy()
        hlavsep_df = self.hlavsep()
        hla_ep = defaultdict(set)
        while len(_epitopes) != 0:
            ind_max = hlavsep_df.Epitope.apply(lambda x: len(x.intersection(_epitopes))).sort_values().index[-1]
            hla = hlavsep_df.iloc[ind_max].HLA
            ep = hlavsep_df.iloc[ind_max].Epitope.intersection(_epitopes)
            hla_ep[hla] = ep
            _epitopes.difference_update(ep)
        return dict(hla_ep)
    