""" """
import os
import regex as re
from typing import Set
from app.epitope import Epitope_DB
from app.analysis_utils import flatten2list, flatten2set

def get_hla_locus(hla:str) -> str:
    """ get the long locus (max 3 letters of gene) of hla """

    gene = hla.split('*')[0]
    return gene if len(gene) == 1 else gene[0:3]

def polymorphic_residues(epitope_set:set, epitope_db) -> set:
    """ Find the aminoacide sequence of each epitope """
    ind = epitope_db.Epitope.apply(lambda x: x in epitope_set)
    poly_residues = epitope_db[ind].PolymorphicResidues.values
    return flatten2list(poly_residues)
    
def hla_to_filename(hla:str):
    """ Translates the hla to the pertinant .pdb file name in the repo """
    locus, specificity = hla.split('*')
    filename = '_'.join([locus, *specificity.split(':')]) + '_V1.pdb'
    return re.split('\d', locus)[0], filename


def find_molecule_path(locus:str, filename:str) -> str:
    """This function makes use of the locus and filename resulted from 'hla_to_filename' function 
        to find the path to the relevant file .pdb file """
    
    path = os.path.expanduser(f'./data/HLAMolecule/{locus[0:2]}') # get until the first 2 character of locus if exist
    pdb_files = [file for file in os.listdir(path) if filename.split('_V1.pdb')[0] in file ]
    if len(pdb_files) != 0:
        return  True, os.path.join(path, f'{pdb_files[0]}')
    else:
        return  False, 'No path exists'



class HLA_Vis:
    def __init__(self, hla:str, epitopes:Set[str]):
        self.hla = hla
        self.epitopes = epitopes
        self.ep_df = Epitope_DB().df 

    # def 
    def hla(epitope_db, 
                    desa_db:pd.DataFrame, 
                    TxIDs:List[int], 
                    rAb=False,
                    mAb=False) -> dict:
    """ This function receives the Tx and desa database and returns HLA vs DESA 
    the output of this function is a nested dictionary. The underscore keys 
    are needed for style parser analysis, while the other keys are useful for presentation
    {TxID:
        {   
            hla:{   
                desa: list of values,
                _desa: list of values,
                desa_rAb: list of tuples
                _desa_rAb: list of tuples
                }
        }
    }
    """
    hlavsdesa = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for TxID in TxIDs:
        desa_db_per_Tx = desa_db[desa_db.TransplantID == TxID]
        if desa_db_per_Tx.shape[0] == 0:
            raise ValueError(f'Transplant ID {TxID} does not exist in the datat set')
        if 'EpvsHLA_Donor_filt' in desa_db_per_Tx.columns:
            desavshla = desa_db_per_Tx.EpvsHLA_Donor_filt.values[0]
        else:
            desavshla = desa_db_per_Tx.EpvsHLA_Donor.values[0]
        for desa, hla in desavshla.items():
            hlavsdesa[TxID][hla]['desa'].append(desa)
            hlavsdesa[TxID][hla]['_desa'].extend(polymorphic_residues(desa, epitope_db))
            if rAb | mAb:
                try:
                    ind = epitope_db.Epitope == desa
                    if epitope_db[ind]['AntibodyReactivity'].values[0] == 'Yes':
                        hlavsdesa[TxID][hla]['desa_rAb'].append(desa)
                        hlavsdesa[TxID][hla]['_desa_rAb'].extend(polymorphic_residues(desa, epitope_db))                    
                    if epitope_db[ind]['mAb'].values[0] == 'Yes':
                        hlavsdesa[TxID][hla]['desa_mAb'].append(desa)
                        hlavsdesa[TxID][hla]['_desa_mAb'].extend(polymorphic_residues(desa, epitope_db))
                except IndexError:
                    logger.info(f'desa {desa} does not exist in the epitope_db')
    return hlavsdesa