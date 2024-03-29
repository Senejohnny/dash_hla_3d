""" This file contains all the methods consumed by Epitope Data Base """
import os
from typing import Union, List
from collections import defaultdict
import pandas as pd
from app.common.logger import logging, get_logger
from app.common.utilities import (
    get_inventory_hlas,
    flatten2list,
    flatten2set,
    flatten_dict_values,
)

class Epitope:
    """ This is a class that entails the data base [Pandas DataFrame] of all epitopes and
        all the related methods tha can be applied to this data base  """

    def __init__(self, path:str='./data/EpitopevsHLA.pickle'):
        # the path is consistent if dash_hla_3d/app.py is ran
        self.path = os.path.expanduser(path)
        # Get hlas with pdb files
        self.pdb_inventory = flatten_dict_values(get_inventory_hlas('./data/HLAMolecule'))
        self.df = pd.read_pickle(self.path) # pylint: disable=invalid-name
        self._hlavsep = None
        self._hlavsep_df = None
        self.log = get_logger('Epitope', logging.INFO)

    def __repr__(self):
        return f""" Epitope_DB(records={len(self.df)}, columns={self.df.columns}) """

    @staticmethod
    def epvshla2hlavsep(epvshla:dict) -> dict:
        """ Transform an ep vs hla dict 2 hla vs ep dict """
        hlavsep = defaultdict(set)
        for epitope, hla in epvshla.items():
            hlavsep[hla].add(epitope)
        return hlavsep

    def filter_mAb(self):
        ind = self.df.mAb == 'Yes'
        self.df = self.df[ind]
        return self

    def is_IgG(self):
        self.isotype()
        if len(self.df != 0):
            return True
        else:
            return False

    def isotype(self, isotype:str='IgG'):
        ind = self.df.isotype.apply(lambda x: isotype in x)
        self.df = self.df[ind]
        return self

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

    def hlavsep(self,
                hla_allel:str='Luminex Alleles',
                only_with_pdb:bool=False,
                ignore_hla:set =set()) -> pd.DataFrame:
        """ returns a DataFrame of HLA vs epitoes
        hla_allel [default is 'Luminex Alleles']: determines the allel type
        only_with_pdb [default is False]: If True, includes only luminex allels that pdb file
        is available:
        only_with_pdb: Include only Luminex Alleles that pdb file is available
        { 'HLA' : {'epitopes'}} """

        if only_with_pdb:
            # Luminex Alleles with available pdb files
            hlas = flatten2set(self.df[hla_allel].values).intersection(self.pdb_inventory) - ignore_hla
        else:
            hlas = flatten2set(self.df[hla_allel].values)

        hlavsep_dict = defaultdict(list)
        for hla in hlas:
            ind = self.df[hla_allel].apply(lambda x: hla in x)
            epitopes = flatten2set(self.df[ind]['Epitope'].values)
            hlavsep_dict['HLA'].append(hla)
            hlavsep_dict['Epitope'].append(epitopes)
        self._hlavsep_df = pd.DataFrame(hlavsep_dict)
        return self._hlavsep_df

    def polymorphic_residues(self, epitope:str) -> set:
        """ Gets the aminoacide sequence of one epitope from
        the polymorphic residue column """
        return self.df[self.df.Epitope == epitope].PolymorphicResidues.values[0]

    def min_hlavsep(self, epitopes:set, ignore_hla:set=set(), most_freq_hla:bool=False) -> dict:
        """ Returns the HLA vs epitope dictionary
            based on minimum number of HLA possible
            ignore_hla: ignores some hla
            format { 'HLA' : {'epitopes'} }
        """
        # Deep copy of epitopes set for later epitope removal
        _epitopes = epitopes.copy()
        hlavsep_df = self.hlavsep(only_with_pdb=True, ignore_hla=ignore_hla)
        if most_freq_hla:
            freq_hla = {
                'A*02:01', 'A*01:01', 'A*03:01', #'A*24:02', 'A*11:01', 'A*68:01', 'A*31:01', 'A*32:01', 'A*26:01', 'A*29:02',
                'B*07:02', 'B*08:01', 'B*15:01', #'B*44:02', 'B*40:01', 'B*35:01', 'B*51:01', 'B*44:03', 'B*57:01', 'B*18:01',
                'DRB1*15:01', 'DRB1*03:01', 'DRB1*01:01', #'DRB1*07:01', 'DRB1*04:01', 'DRB1*13:01', 'DRB1*11:01', 'DRB1*13:02', 'DRB1*04:04', 'DRB1*14:54',
                'DQB1*03:01', 'DQB1*06:02', 'DQB1*02:01', #'DQB1*05:01', 'DQB1*03:02', 'DQB1*02:02', 'DQB1*06:03', 'DQB1*06:04', 'DQB1*03:03', 'DQB1*05:03',
            }
            hlavsep_df = hlavsep_df[hlavsep_df['HLA'].apply(lambda x: x in freq_hla)]
        hla_ep = defaultdict(set)
        for _ in range(10):
            intersect = hlavsep_df.Epitope.apply(
                lambda x: len(x.intersection(_epitopes))
            )
            # find the indexes with maximum value
            value_max = intersect.max()
            if value_max:
                ind_maxes = intersect == value_max
                hlavsep_max_df = hlavsep_df[ind_maxes]
                max_hla = hlavsep_max_df.HLA.values.tolist()[0]
                ind_max = hlavsep_max_df.HLA == max_hla
                set_of_ep = hlavsep_max_df[ind_max].Epitope.values[0].intersection(_epitopes)
                hla_ep[max_hla] = set_of_ep
                _epitopes.difference_update(set_of_ep)
        if _epitopes:
            self.log.info(
                f'Epitopes :{_epitopes} could not be assigned',
                extra={'messagePrefix': 'Epitope.min_hlavsep'}
            )
        return dict(hla_ep)
