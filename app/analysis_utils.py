import pandas as pd
import regex as re
import os 


def polymorphic_residues(epitope_set:set, epitope_db) -> set:
    ind = epitope_db.Epitope.apply(lambda x: x in epitope_set)
    poly_residues = epitope_db[ind].PolymorphicResidues.values
    return [__ for _ in poly_residues for __ in _]
    

# HLA to file name
def hla_to_filename(hla:str):
    """ """
    locus, specificity = hla.split('*')
    filename = '_'.join([locus, *specificity.split(':')]) + '_V1.pdb'
    return re.split('\d', locus)[0], filename

#######################################################
# Find HLA molecule path
#######################################################
def find_molecule_path(locus:str, filename:str) -> str:
    """This function makes use of the locus and filename resulted from 'hla_to_filename' function """
    
    path = os.path.expanduser(f'./data/HLAMolecule/{locus[0:2]}') # get until the first 2 character of locus if exist
    pdb_files = [file for file in os.listdir(path) if filename.split('_V1.pdb')[0] in file ]
    if len(pdb_files) != 0:
        return  True, os.path.join(path, f'{pdb_files[0]}')
    else:
        return  False, ''
