""" """
import os
import json
import logging
import regex as re
from collections import defaultdict
from typing import Set, List
from app.epitope import Epitope_DB
from dash_bio import Molecule3dViewer
from app import styles_parser as sparser, pdb_parser as parser
from app.utilities import (
    flatten2list,
    flatten2set,
    dict_depth,
    get_hla_locus,
    get_hla_polychain,
    polymorphic_residues,
    hla_to_filename,
    find_molecule_path,
)

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


class HLA_Epitope_Vis:

    def __init__(self):
        pass

    def from_epitopes(self, 
                    epitopes, 
                    style:str='sphere', 
                    rAb:bool=False,
                    mAb:bool=False,
        ) -> dict:
        """ 3D visualisation from epitopes as input. This script finds the 
        minimum number of hla on which all the epitopes can be located """

        epitope_obj = Epitope_DB()
        _min_hlavsep = epitope_obj.min_hlavsep(epitopes)
        hlavsep = self._hlavsep_poly(_min_hlavsep, rAb, mAb)
        vis_data = self._get_pdb_data(hlavsep, style)
        self.vis_data = vis_data
        return self
        
    def from_transplant(self, 
                        TxID:int,
                        style:str='sphere', 
                        rAb:bool=False,
                        mAb:bool=False,
        ) -> dict:
        desa_obj = DESA_DB()
        hlavsep = desa_obj.transplant(TxID).get_hlavsdesa()
        hlavsep = self._hlavsep_poly(self.min_hlavsep, rAb, mAb)
        vis_data = self._get_pdb_data(hlavsep, style)
        self.vis_data = vis_data
        return self

    def visualise3D(self, hla:str, TxID:int=None):
        
        model_data = json.loads(self.vis_data[hla]['model'])
        style_data = json.loads(self.vis_data[hla]['style'])
        return self._molecule_viewer(model_data, style_data, opacity=0.6)

    def _hlavsep_poly(self,
                hlavsep:dict, 
                rAb:bool=False,
                mAb:bool=False,
                ) -> dict:
        hlavsep_ = defaultdict(lambda: defaultdict(list))
        for hla, ep in hlavsep.items():
            hlavsep_[hla]['desa'].append(ep)
            hlavsep_[hla]['_desa'].extend(polymorphic_residues(ep, self.ep_db))
            try:
                ind = self.ep_db.Epitope == ep
                if rAb & (self.ep_db[ind]['AntibodyReactivity'].values[0] == 'Yes'):
                    hlavsep_[hla]['desa_rAb'].append(ep)
                    hlavsep_[hla]['_desa_rAb'].extend(polymorphic_residues(ep, self.ep_db))                    
                if mAb & (self.ep_db[ind]['mAb'].values[0] == 'Yes'):
                    hlavsep_[hla]['desa_mAb'].append(ep)
                    hlavsep_[hla]['_desa_mAb'].extend(polymorphic_residues(ep, self.ep_db))
            except IndexError:
                logger.info(f'desa {ep} does not exist in the epitope db')
        return hlavsep_

    def _get_pdb_data(self, hlavsep:dict, style)-> Dict:
        """ Preparing 3D model & style data for visualisation """

        vis_data = defaultdict(lambda: defaultdict())
        desa_info = defaultdict(list)
        for hla in hlavsep.keys():
            desa_info['chain'] = get_hla_polychain(hla)
            desa_info['desa'] = {_[0]:_[1] for _ in hlavsep[hla]['_desa']}
            try:
                desa_info['desa_rAb'] = set( map(lambda x: int(x[0]), hlavsep[hla]['_desa_rAb']) )
            except KeyError:
                desa_info['desa_rAb'] = set()
            try:
                desa_info['desa_mAb'] = set( map(lambda x: int(x[0]), hlavsep[hla]['_desa_mAb']) )
            except KeyError:
                desa_info['desa_mAb'] = set()

            locus, filename = hla_to_filename(hla)
            pdb_exist, pdb_path = find_molecule_path(locus, filename)
            if pdb_exist:
                model_data = parser.create_data(pdb_path)
                style_data = sparser.create_style(pdb_path, style, mol_color='chain', desa_info=desa_info)
                vis_data[hla] = {'model':model_data, 'style':style_data}
            else:
                logger.info(f'{hla}:{pdb_path}: IOError: No such file or directory')
                continue            
        return vis_data

    def _molecule_viewer(self, model, style, opacity=0):
            """ Dashbio output """

            return Molecule3dViewer(
                                    selectionType='atom',
                                    modelData=model,
                                    styles=style,
                                    selectedAtomIds=[],
                                    backgroundOpacity=str(opacity),
                                    atomLabelsShown=False,
                    )