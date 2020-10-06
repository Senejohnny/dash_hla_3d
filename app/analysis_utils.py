import pandas as pd
import regex as re
import os 
import json

from typing import Dict
from collections import defaultdict

from app import styles_parser as sparser
from dash_bio_utils import pdb_parser as parser
import dash_bio as dashbio

def load_data():
    """ This function loads DESA and EpitopeDB data frames"""

    # Load Epitope Database
    ep_path = os.path.expanduser('./data/20200804_EpitopevsHLA_distance.pickle')
    epitope_db = pd.read_pickle(ep_path)

    # Load Epitope Database
    # desa_path = '/Users/Danial/UMCUtrecht/ProcessedData/DSAandDESA/20200916_DESA_new_All.pickle'
    desa_path = os.path.expanduser('./data/desa_db.pickle')
    desa_db = pd.read_pickle(desa_path)
    return epitope_db, desa_db

def polymorphic_residues(epitope_set:set, epitope_db) -> set:
    """ Find the aminoacide sequence of each epitope """
    ind = epitope_db.Epitope.apply(lambda x: x in epitope_set)
    poly_residues = epitope_db[ind].PolymorphicResidues.values
    return [__ for _ in poly_residues for __ in _]
    
def hla_to_filename(hla:str):
    """ Translates the hla to the pertinant .pdb file name in the repo"""
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


def hlavsdesa_donor(desa_db:pd.DataFrame, TxID:int) -> dict:
    """ This function receives the Tx and desa database and returns HLA vs DESA """

    desa_3d_view = desa_db[desa_db.TransplantID == TxID]
    desavshla = desa_3d_view.EpvsHLA_Donor.values[0]
    hlavsdesa = defaultdict(list)
    for key, value in desavshla.items():
        hlavsdesa[value].append(key)
    return hlavsdesa

def get_hla_locus(hla:str) -> str:
    """ get the long locus (max 3 letters of gene) of hla """

    gene = hla.split('*')[0]
    return gene if len(gene) == 1 else gene[0:3]
    
def get_hla_polychain(hla:str)-> str:
    """ get the ploymorphic chain of hla """

    polychain =  {'A': 'A', 'B': 'A', 'C': 'A',
                'DRB': 'B', 'DQA': 'A', 'DQB':'B'}
    return polychain.get(get_hla_locus(hla), 'hla is invalid')

def _3dview_data_preparation(hlavsdesa:dict)-> Dict:
    """ Data Preparation (model & style) for the 3d viewer """

    _3d_models_data = defaultdict()
    for hla in hlavsdesa.keys():
        desa_info = {'chain': get_hla_polychain(hla),
        'relevant_desa': set([int(_[0]) for _ in hlavsdesa.get(hla)]),
        'irrelevant_desa': set(),
        }

        locus, filename = hla_to_filename(hla)
        pdb_exist, pdb_path = find_molecule_path(locus, filename)
        # print(pdb_path)
        # Create the model data from the pdb files 
        if pdb_exist:
            model = parser.create_data(pdb_path)
        else:
            print(f'{hla}:{pdb_path}: IOError: No such file or directory')
        # Create the cartoon style from the decoded contents
        style = sparser.create_style(pdb_path, style='sphere', mol_color='chain', desa_info=desa_info)
        _3d_models_data[hla] = {'model':model, 'style':style}
    return _3d_models_data


def _3d_dashbio(id:str, model, style, opacity=0):
    """ Dashbio output """

    return dashbio.Molecule3dViewer(
                                    id=id,
                                    selectionType='atom',
                                    modelData=model,
                                    styles=style,
                                    selectedAtomIds=[],
                                    backgroundOpacity=str(opacity),
                                    atomLabelsShown=False,
                                    )


def div_3dviewer(hlas:list, hla_ind:int, _3d_data, _hlavsdesa: dict):
    """ This function returns the hla name and 3d structure components 
        that are needed by dash component for visualisztion  """

    if hla_ind < len(hlas):
        hla = hlas[hla_ind]
        model = json.loads(_3d_data[hla]['model'])
        style =  json.loads(_3d_data[hla]['style'])
        component = _3d_dashbio(hla, model, style, opacity=0.6)
        # id  = 'mol3d-viewer' + '-' + str(hla_ind + 1)
        return hla, component
    else:
        return hla_ind + 1, None


def data_3dviewer(desa_db, epitope_db, TxID:int):
    """ This function provides the require data consumed by 'div_3dviewer' """

    hlavsdesa = hlavsdesa_donor(desa_db, TxID)
    _hlavsdesa = {key:polymorphic_residues(value, epitope_db) for key, value in hlavsdesa.items()}
    _3d_data = _3dview_data_preparation(hlavsdesa)
    hlas = list(hlavsdesa.keys())
    return hlas, _3d_data, _hlavsdesa, hlavsdesa