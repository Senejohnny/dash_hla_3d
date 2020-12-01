import pandas as pd
import regex as re
import os 
import json
import sys
import logging

from typing import Dict, List
from collections import defaultdict

from app import styles_parser as sparser, pdb_parser as parser
# from dash_bio_utils import pdb_parser as parser
import dash_bio as dashbio

######### logging ########

def setup_logger(name, log_file, level=logging.DEBUG):
    """To setup as many loggers as you want"""
    formatter = logging.Formatter('%(name)s - %(asctime)s - %(levelname)s - %(message)s')
    handler = logging.FileHandler(log_file)        
    handler.setFormatter(formatter)

    # Create & configure logger
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

logger = setup_logger('main_loger', 'data/log.log')



###################################################################
# Auxiliary functions for analysis
###################################################################
def load_data():
    """ This function loads DESA and EpitopeDB data frames"""

    return load_epitope_db(), load_desa_db()

def load_desa_db():
    """ This function loads DESA data frame"""

    desa_path = os.path.expanduser('./data/desa_3d_view.pickle')
    desa_db = pd.read_pickle(desa_path)
    # desa_db = desa_db.drop('Status', axis=1)
    return desa_db

def load_epitope_db():
    """ This function loads EpitopeDB data frame"""

    ep_path = os.path.expanduser('./data/20201123_EpitopevsHLA.pickle')
    epitope_db = pd.read_pickle(ep_path)
    return epitope_db

def flatten2list(object):
    """ This function flattens objects in a nested structure and returns a list"""
    gather = []
    for item in object:
        if isinstance(item, (list, set)):
            gather.extend(flatten2list(item))            
        else:
            gather.append(item)
    return gather

def flatten2set(object) -> set:
    """ This function flattens objects in a nested structure and returns a set"""

    return set(flatten2list(object))

###################################################################
# Main functions for analysis
###################################################################

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


def hlavsdesa_donor(epitope_db, 
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

    
def get_hla_polychain(hla:str)-> str:
    """ get the ploymorphic chain of hla locus """

    polychain =  {'A': 'A', 'B': 'A', 'C': 'A',
                'DRB': 'B', 'DQA': 'A', 'DQB':'B'}
    return polychain.get(get_hla_locus(hla), 'hla is invalid')

def _3dview_data_preparation(hlavsdesa:dict, style)-> Dict:
    """ Data Preparation (model & style) for the 3d viewer """

    _3d_models_data = defaultdict(lambda: defaultdict())
    for TxID in hlavsdesa.keys():
        for hla in hlavsdesa[TxID].keys():
            desa_info = {
                'chain': get_hla_polychain(hla),
                'desa': {_[0]:_[1] for _ in hlavsdesa[TxID][hla]['_desa']},
                'desa_rAb': set([int(_[0]) for _ in hlavsdesa[TxID][hla]['_desa_rAb']]),
                'desa_mAb': set([int(_[0]) for _ in hlavsdesa[TxID][hla]['_desa_mAb']])
            }

            locus, filename = hla_to_filename(hla)
            pdb_exist, pdb_path = find_molecule_path(locus, filename)
            # Create the model data from the pdb files 
            if pdb_exist:
                model_data = parser.create_data(pdb_path)
            else:
                logger.info(f'{hla}:{pdb_path}: IOError: No such file or directory')
                continue
            # Create the style data from the decoded contents
            # style_data = sparser.create_style(pdb_path, style, desa_info=desa_info)
            style_data = sparser.create_style(pdb_path, style, mol_color='chain', desa_info=desa_info)
            _3d_models_data[TxID][hla] = {'model':model_data, 'style':style_data}
    return _3d_models_data


def _3d_dashbio(model, style, opacity=0):
    """ Dashbio output """

    return dashbio.Molecule3dViewer(
                                    selectionType='atom',
                                    modelData=model,
                                    styles=style,
                                    selectedAtomIds=[],
                                    backgroundOpacity=str(opacity),
                                    atomLabelsShown=False,
                                    )


def div_3dviewer(hla:str, TxID, _3d_data):
    """ This function returns the hla name and 3d structure components 
        that are needed by dash component for visualisztion  
        _hlavsdesa: dict can added later for testing 
        """
    
    model_data = json.loads(_3d_data[TxID][hla]['model'])
    style_data = json.loads(_3d_data[TxID][hla]['style'])
    component = _3d_dashbio(model_data, style_data, opacity=0.6)
    return component


def data_3dviewer(desa_db, 
                epitope_db, 
                TxIDs:List[int], 
                style:str='sphere',
                rAb=False,
                mAb=False):
    """ This function provides the require data consumed by 'div_3dviewer' """

    hlavsdesa = hlavsdesa_donor(epitope_db, desa_db, TxIDs, rAb, mAb)
    _3d_data = _3dview_data_preparation(hlavsdesa, style)
    return _3d_data, hlavsdesa