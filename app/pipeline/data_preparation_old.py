""" This file contains the data preparation pipeline for Transplants with DESA """
from collections import Counter
import pandas as pd
from app.common.utilities import get_exp

transplant = pd.read_csv('~/UMCUtrecht/RawData/NOTR/PRO_NOTR_TRANSPLANT.csv', sep=';')
desa_path = '/Users/Danial/UMCUtrecht/ProcessedData/DSAandDESA/20200916_DESA_new_All.pickle'
desa = pd.read_pickle(desa_path)
desa = desa.merge(transplant[['TransplantID', 'GraftFunction_NOTR', 'FailureCode_NOTR', 'GraftSurvivalYears_NOTR']], on='TransplantID')
ind_desa = desa.Status == 'DESA'
desa_yes = desa[ind_desa]

def get_hla_class(hla:str):
    """ get hla class from high resolution typing """
    if hla.split('*')[0] in ['A', 'B', 'C']:
        return 'I'
    return 'II'

def get_class(x):
    """ get hla class """
    _class = {get_hla_class(hla) for hla in set(x.values())}
    return ','.join(list(_class))

desa_yes = desa_yes.assign(Donor_HLAvsnumDESA= desa_yes.EpvsHLA_Donor\
                .apply(lambda x: Counter([value for key,  value in x.items()]) if x else 0),
                Donor_HLA= desa_yes.EpvsHLA_Donor\
                .apply(lambda x: {_.split(':')[0] for _ in set(x.values())}),
                Donor_HLA_Class= desa_yes.EpvsHLA_Donor.apply(get_class),
                Donor_HLA_exp = desa_yes.EpvsHLA_Donor.apply(
                    lambda x: {hla:get_exp(hla) for hla in set(x.values())}
                )
            )
desa_yes.rename(columns={
                    'FailureCode_NOTR': 'Failure',
                    'GraftSurvivalYears_NOTR': 'Survival[Y]',
                    'GraftFunction_NOTR': 'Graft_Func',
                }, inplace=True)
print(desa_yes.columns)
desa_yes.to_pickle('/Users/Danial/Repos/STRIDE/dash_hla_3d/data/desa_3d_view2.pickle')
