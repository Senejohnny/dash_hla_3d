""" This file contains all the methods consumed by Epitope Data Base """
from typing import Union, List
from collections import defaultdict
import pandas as pd
from app.analysis_utils import flatten2set

class Epitope:
    """ This is a class that entails the data base [Pandas DataFrame] of all epitopes and
        all the related methods tha can be applied to this data base  """

    def __init__(self, path:str='./data/20201123_EpitopevsHLA.pickle'):
        # the path is consistent if dash_hla_3d/app.py is ran
        self.df = pd.read_pickle(path) # pylint: disable=invalid-name
        self._hlavsep = None
        self._hlavsep_df = None

    def __repr__(self):
        return f""" Epitope_DB(records={len(self.df)}, columns={self.df.columns}) """

    def get_epitopes(self, value:Union[str, List[str]]):
        """ get epitope info from the df
        value: can be str or a list of strings """

        if isinstance(value, str):
            ind = self.df.Epitope == value
        else:
            ind = self.df.Epitope.apply(lambda x: x in value)
        self.df = self.df[ind]
        return self

    def ellipro(self, value):
        """ filter EpitopeDB based on desired ellipro score """
        if isinstance(value, str):
            ind = self.df.Epitope == value
        else:
            ind = self.df['ElliPro Score'].apply(lambda x: x in value)
        self.df = self.df[ind]
        return self
    
    def hlavsep(self, hla_allel:str='Luminex Alleles') -> pd.DataFrame:
        """ returns a DataFrame { 'HLA' : {'epitopes'}} """
        
        hlas = flatten2set(self.df[hla_allel].values)
        hlavsep_dict = defaultdict(list)
        for hla in hlas:
            ind = self.df[hla_allel].apply(lambda x: hla in x)
            epitopes = flatten2set(self.df[ind]['Epitope'].values)
            hlavsep_dict['HLA'].append(hla)
            hlavsep_dict['Epitope'].append(epitopes)
        self._hlavsep_df = pd.DataFrame(hlavsep_dict)
        return self._hlavsep_df

    def min_hlavsep(self, epitopes:set) -> dict:
        """ Returns the HLA vs epitope dictionary
            based on minimum number of HLA possible
            format { 'HLA' : {'epitopes'} }
        """
        # Deep copy of epitopes set for later epitope removal
        _epitopes = epitopes.copy()
        hlavsep_df = self.hlavsep()
        hla_ep = defaultdict(set)
        while len(_epitopes) != 0:
            ind_max = hlavsep_df.Epitope.apply(
                lambda x: len(x.intersection(_epitopes))
            ).sort_values().index[-1]
            hla = hlavsep_df.iloc[ind_max].HLA
            set_of_ep = hlavsep_df.iloc[ind_max].Epitope.intersection(_epitopes)
            hla_ep[hla] = set_of_ep
            _epitopes.difference_update(set_of_ep)
        return dict(hla_ep)

    def epvshla2hlavsep(self, epvshla:dict) -> dict:
        """ Transform an ep vs hla dict 2 hla vs ep dict """
        hlavsep = defaultdict(set)
        for epitope, hla in epvshla.items():
            hlavsep[hla].add(epitope)
        return hlavsep
