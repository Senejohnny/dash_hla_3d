"""PDB parser
This module contains functions that can read PDB files and return a
JSON representation of the structural data."""

import re
import json
import os
from shutil import copy2

import parmed as pmd


def create_data(pdb_path):
    """
    Parse the protein data bank (PDB) file to generate
    input modelData
    @param pdb_path
    Name of the biomolecular structure file in PDB format
    """

    # Create local copy of temp file
    copy2(pdb_path, './tmp.pdb')

    # Use parmed to read the bond information from temp file
    # top = pmd.load_file('tmp.pdb')
    top = pmd.load_file(pdb_path)

    # Remove the created temp file
    os.remove('tmp.pdb')

    # Read PDB file to create atom/bond information
    with open(pdb_path, 'r') as infile:
        # store only non-empty lines
        lines = [l.strip() for l in infile if l.strip()]

    # Initialize all variables
    var_nchains = []
    serial = []
    atm_name = []
    res_name = []
    chain = []
    res_id = []
    positions = []
    occupancy = []
    temp_factor = []
    atom_type = []
    ct = 0

    datb = {
        'atoms': [],
        'bonds': []
    }

    # Variables that store the character positions of different
    # parameters from the molecule PDB file
    serialpos = slice(6,11) #[6, 11]
    atm_namepos = slice(12, 16) #[12, 16]
    r_namepos = slice(17, 20) #[17, 20]
    chainpos = slice(21, 22) #[21, 22]
    r_idpos = slice(22, 26) #[22, 26]
    xpos = slice(30, 38) #[30, 38]
    ypos = slice(38, 46) #[38, 46]
    zpos = slice(46, 54) #[46, 54]
    occupos = slice(54, 60) #[54, 60]
    bfacpos = slice(60, 66) #[60, 66]
    atm_typepos = slice(77,79) #[77, 79]

    for l in lines:
        line = l.split()
        if line[0] in ["ATOM", "HETATM"]:
            serial.append(int(l[serialpos]))
            atm_name.append(l[atm_namepos].strip())
            val_r_name = l[r_namepos].strip()
            res_name.append(val_r_name)
            chain_val = l[chainpos].strip()
            chain.append(chain_val)
            if chain_val not in var_nchains:
                var_nchains.append(chain_val)
            val_r_id = int(l[r_idpos])
            res_id.append(val_r_id)
            x = float(l[xpos])
            y = float(l[ypos])
            z = float(l[zpos])
            positions.append([x, y, z])
            occupancy.append(l[occupos].strip())
            temp_factor.append(l[bfacpos].strip())
            atom_type.append(l[atm_typepos].strip())
            ct += 1

    # Create list of atoms
    tmp_res = res_id[0]
    resct = 1
    for i in range(len(chain)):  # pylint: disable=consider-using-enumerate
        if tmp_res != res_id[i]:
            tmp_res = res_id[i]
            resct += 1
        datb['atoms'].append({
            "name": atm_name[i],
            "chain": chain[i],
            "positions": positions[i],
            "residue_index": resct,
            "element": atom_type[i],
            "residue_name": res_name[i] + str(res_id[i]),
            "serial": i,
        })

    # Create list of bonds using the parmed module
    for i in range(len(top.bonds)):
        bondpair = top.bonds[i].__dict__
        atom1 = re.findall(r"\[(\d+)\]", str(bondpair['atom1']))
        atom2 = re.findall(r"\[(\d+)\]", str(bondpair['atom2']))
        datb['bonds'].append({
            'atom2_index': int(atom1[0]),
            'atom1_index': int(atom2[0])
        })

    return json.dumps(datb)