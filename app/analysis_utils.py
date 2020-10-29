import pandas as pd
import regex as re
import os 
import json

from typing import Dict
from collections import defaultdict

from app import styles_parser as sparser
from dash_bio_utils import pdb_parser as parser
import dash_bio as dashbio


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
    return desa_db

def load_epitope_db():
    """ This function loads EpitopeDB data frame"""

    ep_path = os.path.expanduser('./data/20201027_EpitopevsHLA_distance.pickle')
    epitope_db = pd.read_pickle(ep_path)
    return epitope_db

def flatten2list(object):
    """ This function flattens objects in a nested structure """
    gather = []
    for item in object:
        if isinstance(item, (list, set)):
            gather.extend(flatten2list(item))            
        else:
            gather.append(item)
    return gather

###################################################################
# Main functions for analysis
###################################################################

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
        return  False, ''


def get_hla_locus(hla:str) -> str:
    """ return the locus of hla"""

    star_location = hla.find('*')
    return hla[:star_location]


def hlavsdesa_donor(epitope_db, desa_db:pd.DataFrame, TxID:int, rAb) -> dict:
    """ This function receives the Tx and desa database and returns HLA vs DESA 
    the output of this function is a nested dictionary. The underscore keys 
    are needed for style parser analysis, while the other keys are useful for presentation
    {
        hla:{   desa: list of values,
                _desa: list of values,
                desa_rAb: list of tuples
                _desa_rAb: list of tuples
            }
    }
    """

    desa_3d_view_db = desa_db[desa_db.TransplantID == TxID]
    if desa_3d_view_db.shape[0] == 0:
        raise ValueError(f'Transplant ID {TxID} does not exist in the datat set')
    desavshla = desa_3d_view_db.EpvsHLA_Donor.values[0]
    hlavsdesa = defaultdict(lambda: defaultdict(list))
    for desa, hla in desavshla.items():
        hlavsdesa[hla]['desa'].append(desa)
        hlavsdesa[hla]['_desa'].extend(polymorphic_residues(desa, epitope_db))
        if rAb:
            try:
                if epitope_db[epitope_db.Epitope == desa]['AntibodyReactivity'].values[0] == 'Yes':
                    hlavsdesa[hla]['desa_rAb'].append(desa)
                    hlavsdesa[hla]['_desa_rAb'].extend(polymorphic_residues(desa, epitope_db))
            except IndexError:
                hlavsdesa[hla]['desa_rAb'] = []
                hlavsdesa[hla]['_desa_rAb'] = []
    return hlavsdesa

def get_hla_locus(hla:str) -> str:
    """ get the long locus (max 3 letters of gene) of hla """

    gene = hla.split('*')[0]
    return gene if len(gene) == 1 else gene[0:3]
    
def get_hla_polychain(hla:str)-> str:
    """ get the ploymorphic chain of hla locus """

    polychain =  {'A': 'A', 'B': 'A', 'C': 'A',
                'DRB': 'B', 'DQA': 'A', 'DQB':'B'}
    return polychain.get(get_hla_locus(hla), 'hla is invalid')

def _3dview_data_preparation(hlavsdesa:dict, style, rAb)-> Dict:
    """ Data Preparation (model & style) for the 3d viewer """

    _3d_models_data = defaultdict()
    print(dict(hlavsdesa))
    for hla in hlavsdesa.keys():
        desa_info = {
            'chain': get_hla_polychain(hla),
            'desa': set([int(_[0]) for _ in hlavsdesa[hla]['_desa'] ]),
            'desa_rAb': set([int(_[0]) for _ in hlavsdesa[hla]['_desa_rAb'] ])
        }

        print('desa_info', desa_info)
        locus, filename = hla_to_filename(hla)
        pdb_exist, pdb_path = find_molecule_path(locus, filename)
        # Create the model data from the pdb files 
        if pdb_exist:
            model_data = parser.create_data(pdb_path)
        else:
            print(f'{hla}:{pdb_path}: IOError: No such file or directory')
            continue
        # Create the cartoon style from the decoded contents
        style_data = sparser.create_style(pdb_path, style, mol_color='chain', desa_info=desa_info)
        _3d_models_data[hla] = {'model':model_data, 'style':style_data}
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


def div_3dviewer(hla:str, _3d_data):
    """ This function returns the hla name and 3d structure components 
        that are needed by dash component for visualisztion  
        _hlavsdesa: dict can added later for testing 
        """
    
    model_data = json.loads(_3d_data[hla]['model'])
    style_data = json.loads(_3d_data[hla]['style'])
    component = _3d_dashbio(model_data, style_data, opacity=0.6)
    return component


def data_3dviewer(desa_db, epitope_db, TxID:int, style:str='sphere', rAb=False):
    """ This function provides the require data consumed by 'div_3dviewer' """

    hlavsdesa = hlavsdesa_donor(epitope_db, desa_db, TxID, rAb)
    # _hlavsdesa = {hla:{key:polymorphic_residues(value, epitope_db)} for hla, val_dict in hlavsdesa.items() for key, value in val_dict.items()}
    # print(_hlavsdesa )
    _3d_data = _3dview_data_preparation(hlavsdesa, style, rAb)
    return _3d_data, hlavsdesa