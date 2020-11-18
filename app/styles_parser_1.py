"""Styles parser
This module contains the functions that can specify the style of a
molecule that is rendered using the Molecule3dViewer component in Dash
Bio.
Modified October 2020 by Danial Senejohnny, d.senejohnny@gmail.com """

import json
import os
import logging
# from app.analysis_utils import setup_logger
# logger = setup_logger('aminoacid', 'data/aminoacid.log')

logging.basicConfig(filename= 'data/aminoacid.log',
                    filemode = 'w',
                    format= '%(name)s - %(asctime)s - %(levelname)s - %(message)s',
                    level=logging.DEBUG,
)
# Create & configure logger
logger = logging.getLogger('aminoacid')

Aminoacid_conversion = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 
                        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 
                        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',  
                        'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 
                        'LYS': 'K', 'ASN': 'N', 'TRP': 'W', 'MET': 'M'}

chain_colors = {
    "A": "#65A5E2", # Blue
    "B": "#CA7FE5", # purple
    "C": "#65E2AB", # green
}
desa_color = {
    'desa': '#FFFF00', # Yellow
    'desa_mAb': '#FFA500', # Orange
    'desa_rAb': '#FF0000', # Red
}

# class StyleParser:
#     def __init__(self,  pdb_path:str):
#         self.path = pdb_path
    
def create_style(pdb_path:str, 
                style:str,
                desa_info:dict):
    """Function to create the different styles (stick, cartoon, sphere)
    using the protein data bank (PDB) file as input. This function outputs
    the styles as a JSON file
    @param pdb_path
    Name of the biomolecular structure file in PDB format
    @param style
    Type of representation of the biomolecule (options: stick, cartoon, sphere)
    @param mol_color
    Coloring scheme for depicting biomolecules (options: residue_type, atom, residue, chain)
    @param custom_dict
    optional parameter to specify the color scheme for different chains
    in JSON format
    @param atm_color
    optional parameter to specify the color scheme for different atoms
    in JSON format
    """

    # Read input file
    with open(pdb_path, 'r') as pdb_file:
        # store only non-empty and relevant lines (comments start not with ATOM or HETATM)
        lines = [line.strip() for line in pdb_file if line.strip() if line.split()[0] in ["ATOM", "HETATM"]]

   
    # Initialize variables
    data = {}
    chains = []
    atm_types = []
    res_names = []
    res_indxs = []

    # Variables that store the character positions of different
    # parameters from the molecule PDB file
    res_name_ind = slice(17, 20)
    chain_ind = slice(21, 22)
    res_index_ind = slice(23, 26)
    atom_type_ind = slice(77, 78)

    # print(len(lines)),
    for line in lines:
        res_name = line[res_name_ind].strip()
        chain = line[chain_ind]
        res_indx = line[res_index_ind].strip()
        atm_type = line[atom_type_ind]
        # print(res_name, chain, res_indx, atm_type)

        chains.append(chain)
        atm_types.append(atm_type)
        res_names.append(res_name)
        res_indxs.append(res_indx)
        index = len(chains) - 1

        data[index] = {
                    'color': chain_colors.get(chain, '#BEA06E'),
                    'visualization_type': style
                }
        if chain == desa_info['chain']:
            if res_indx in desa_info['desa'].keys():
                data[index] = {
                        "color": desa_color['desa'],
                        "visualization_type": style
                    }
                try:
                    found_res_name, exp_res_name = Aminoacid_conversion.get(res_name), desa_info['desa'][res_indx]
                    assert found_res_name == exp_res_name
                except AssertionError:
                    logger.info(f"""In the .pdb file {pdb_path}, chain:{chain}, 
                        Expected {res_indx}{exp_res_name}, found: {res_indx}{found_res_name}""")
            if int(res_indx) in desa_info['desa_rAb']: # we put if to recolor the epitopes in the former step
                data[index] = {
                        "color": desa_color['desa_rAb'],
                        "visualization_type": style
                    }
            if int(res_indx) in desa_info['desa_mAb']: # we put if to recolor the epitopes in the former step
                data[index] = {
                        "color": desa_color['desa_mAb'],
                        "visualization_type": style
                    }
        # print(json.dumps(data))
    return json.dumps(data)

if __name__ == '__main__':
    
    path = os.path.expanduser('./data/HLAMolecule/A/A_01_01_V1.pdb')
    create_style(path, None, None, 'sphere')

    def test_style_parser():
        pass