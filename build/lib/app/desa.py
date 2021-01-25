""" DESA class """
import pandas as pd

class DESA:
    """ This is a class that entails the data base [Pandas DataFrame] of all transplants
        with DESA and all the related methods that can be applied or update this DataFrame """
    
    def __init__(self, path:str='../data/desa_3d_view.pickle'):
        self.df = pd.read_pickle(path)

    def __repr__(self):
        return f""" DESA_DB(records={len(self.df)}, columns={self.df.columns}) """
    
    def __str__(self):
        return __repr__()
    
    def donor_type(self, donor_type:str='Deceased'):
        if donor_type not in ['Living', 'Deceased']:
            raise KeyError(f'{donor_type} does not exist in the df values, accepted values: {self.df.Donor_Type.unique()}')
        ind = self.df.Donor_Type.apply(lambda x: x == donor_type)                          
        self.df = self.df[ind]
        return self
    
    def get_tx(self, TxID:int) -> pd.DataFrame:
        ind = self.df.TransplantID == TxID
        if sum(ind) == 0:
            raise ValueError(f'Transplant ID {TxID} does not exist in the datat set')
        else:
            return self.df[ind]
        
    def hla_class(self, hla_class):
        if hla_class not in ['I', 'II', 'I,II']:
            raise KeyError(f'{hla_class} does not exist in the df values, accepted values: {self.df.Donor_HLA_Class.unique()}')
        ind = self.df.Donor_HLA_Class.apply(lambda x: x == hla_class)                          
        self.df = self.df[ind]
        return self
        
    def early_failed(self, threshold):
        ind_t = self.df['Survival[Y]'].apply(lambda x: x < threshold)
        ind_e = self.df.Failure.apply(lambda x: x == 1)
        self.df = self.df[ind_t & ind_e]
        return self

    def late_failed(self, threshold):
        ind_t = self.df['Survival[Y]'].apply(lambda x: x > threshold)
        ind_e = self.df.Failure.apply(lambda x: x != 1 )
        self.df = self.df[ind_t & ind_e]
        return self

if __name__ == '__main__':
    desa = DESA()
    print(desa.df.columns)
    print(desa.df.EpvsHLA_Donor)
# print('getcwd:      ', os.getcwd())