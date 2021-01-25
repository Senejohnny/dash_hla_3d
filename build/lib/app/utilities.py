import os
import regex as re
from typing import Set

def flatten2list(object):
    """ This function flattens objects in a nested structure and retu"""
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

def flatten_dict_values(dictionary:dict) -> set:
    """ This function flattens objects in a nested structure and returns a set"""

    return flatten2set(dictionary.values())

def get_hla_locus(hla:str) -> str:
    """ get the long locus (max 3 letters of gene) of hla """

    gene = hla.split('*')[0]
    return gene if len(gene) == 1 else gene[0:3]

def get_hla_polychain(hla:str)-> str:
    """ get the ploymorphic chain of hla locus """

    polychain =  {'A': 'A', 'B': 'A', 'C': 'A',
                'DRB': 'B', 'DQA': 'A', 'DQB':'B'}
    return polychain.get(get_hla_locus(hla), 'hla is invalid')

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

def dict_depth(my_dict):
    """ function to find the depth of dictionary """ 
    
    if isinstance(my_dict, dict): 
        return 1 + (max(map(dict_depth, my_dict.values() )) if my_dict else 0 ) 
    return 0