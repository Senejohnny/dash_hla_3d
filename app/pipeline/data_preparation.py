""" This file contains the data preparation pipeline for Transplants with DESA """
from collections import Counter
import pandas as pd
from app.common.utilities import get_hla_exp, get_class


class DataPreparation:
    """
    A class to prepare any analysis data set

    Feature:
    --------
        - Include the antibody analysis from the pipeline
        - Include other required data such as demographical and clinical
        -
    """
    def __init__(self):
        self.donor, self.km, self.tx, self.desa = self.load_datasets()

    @staticmethod
    def load_datasets():
        """
        Important columns:
        - NOTR transplant dataset:
            'TransplantID', 'DeelnemerID', 'DonorID', 'RecipientAge_NOTR','DonorAge_NOTR',
            'CauseGraftFailure_NOTR', 'GraftSurvivalYears_NOTR',
            'TransplantOrderNumber_NOTR': n-th transplanted kidney for the recipient,

            GraftFunction_NOTR: Function of graft directly after transplantation
                - DIR: Directly functioned grafts
                - NEV: Never functioned grafts
                - DEL: delayed functioned grafts
                - NDI:
                - UNK: Unknown

            'ColdIschaemicPeriod_NOTR' (minute, 5940:unknown),
            'FailureCode_NOTR':
                0 : no failure,
                1 : failure,
                2 : censoring event [death with functioning graft]),

        - NOTR donor dataset
            'DonorID', 'DonorSex_NOTR', 'TypeOfDonor_NOTR', 'TypeCadaveric_NOTR',
            'DonorCauseDeath_NOTR', 'DonorHLAMergedType_NOTR'

        - NOTR recipient dataset
            'DeelnemerID', 'RecipientSex_NOTR', RecipientHLAMergedType_NOTR

        'Retransplant', 'TxYear', 'Center',
        IL2rMoAb, InductionTherapy: this is on induction therapy
        'should' transplant

        DialysisYears

        - Pipeline antibody analysis:

        """
        desa_km = pd.read_csv('~/UMCUtrecht/KaplanMeier/DESAsurvival_original.csv', sep=';')
        notr_tx_path = '~/UMCUtrecht/RawData/NOTR/PRO_NOTR_TRANSPLANT.csv'
        transplant = pd.read_csv(notr_tx_path, sep=';')
        notr_donor_path = '~/UMCUtrecht/RawData/NOTR/PRO_NOTR_DONOR.csv'
        donor = pd.read_csv(notr_donor_path, sep=';')
        transplant = pd.read_csv(notr_tx_path, sep=';')
        desa_path = '~/UMCUtrecht/ProcessedData/DSAandDESA/20200916_DESA_new_All.pickle'
        desa = pd.read_pickle(desa_path)
        return donor, desa_km, transplant, desa

    def engineering(self):
        """
        change cold schemia time from minute to hour
        find the period patient was on dialysis
        """
        {
            'desa': 'pass',
            'tx': 'pass'
        }

        ind_desa = self.desa.Status == 'DESA'
        self.desa = self.desa[ind_desa]
        self.desa = self.desa.assign(
            Donor_HLAvsnumDESA = self.desa.EpvsHLA_Donor.apply(
                    lambda x: Counter([value for key,  value in x.items()]) if x else 0
                ),
            Donor_HLA = self.desa.EpvsHLA_Donor.apply(
                    lambda x: {_.split(':')[0] for _ in set(x.values())}
                ),
            Donor_HLA_Class= self.desa.EpvsHLA_Donor.apply(get_class),
            Donor_HLA_exp = self.desa.EpvsHLA_Donor.apply(
                    lambda x: {hla:get_hla_exp(hla) for hla in set(x.values())}
                ),
        )
        return self

    def merge(self):
        self.desa = self.desa\
            .merge(
                self.tx[['TransplantID', 'GraftFunction_NOTR', 'FailureCode_NOTR',
                        'GraftSurvivalYears_NOTR']], on='TransplantID')\
            .merge(
                self.km[['TransplantID', 'TypeOfDonor_NOTR',]], on='TransplantID')
        return self

    def modify_columns(self):

        self.desa.rename(
            columns={
                'FailureCode_NOTR': 'Failure',
                'GraftSurvivalYears_NOTR': 'Survival[Y]',
                'TypeOfDonor_NOTR': 'Donor_Type',
                'GraftFunction_NOTR': 'Graft_Func',
                # 'ColdIschaemicPeriod_NOTR': 'CIP_hour',
            }, inplace=True)

        return self

    def save_dataset(self, path=None, dataset='desa', ):
        if path is None:
            path = '~/UMCUtrecht/ProcessedData/20210310_prepared_dataset.pickle'
        {
            'desa': self.desa,
            'tx': self.tx,
        }.get(dataset).to_pickle(path)


if __name__ == '__main__':
    data = DataPreparation()
    path = '/Users/Danial/Repos/STRIDE/dash_hla_3d/data/desa_3d_view2.pickle'
    data.engineering().merge().modify_columns().save_dataset(path)
    print(data.desa.columns)
    # print(data.desa.loc[data.desa.Graft_Func == 'NEV', ('Graft_Func', 'Survival[Y]')])
