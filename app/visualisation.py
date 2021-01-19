""" This file encompass all the 3D visualisation scripts"""
import json
import logging
from typing import Set, Dict
from collections import defaultdict
from dash_bio import Molecule3dViewer
import dash_bootstrap_components as dbc
import dash_html_components as html
from app.common.logger import get_logger

from app.epitope import Epitope
from app.desa import DESA
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



SERVICE_NAME = 'VisualiseHLA'

class VisualiseHLA:

    def __init__(self):
        """ uses DESA and Epitope Classes """

        self.DESA = DESA()
        self.Epitope = Epitope()
        self.hlavsep = None
        self.txvshlavsep = None
        self.log = get_logger(SERVICE_NAME, logging.INFO)

    @staticmethod
    def _molecule_viewer(model:dict, style:dict, opacity=0):
        """ Dashbio front-end visualisation function """
        return Molecule3dViewer(            # pylint: disable=not-callable
            selectionType='atom',
            modelData=model,
            styles=style,
            selectedAtomIds=[],
            backgroundOpacity=str(opacity),
            atomLabelsShown=False,
            )
        
    def _get_vis_data_from_pdb(self, hlavsep:dict, style, called_by)-> Dict[str, Dict[str, list]]:
        """
        Extract visualisation data like model & style from pdb file and construct a
        nested dictionary in the form of
            {'hla': {'model': , 'style':} }
        the consumer of this intermediary nested dict is visualise3D mehtod.
        """
        vis_data = defaultdict(lambda: defaultdict())
        desa_info = defaultdict(list)
        for hla in hlavsep.keys():
            self.log.info(f'Start getting .pdb files for HLA {hla}', extra={'messagePrefix': called_by})
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
                try:
                    model_data = parser.create_data(pdb_path)
                    style_data = sparser.create_style(pdb_path, style, mol_color='chain', desa_info=desa_info)
                    vis_data[hla] = {'model':model_data, 'style':style_data}
                    self.log.info(f'Successfully loaded & parsed this file', extra={'messagePrefix': called_by})
                except:
                    self.log.info(f'Can not parse .pdb for HLA {hla}', extra={'messagePrefix': called_by})
                    raise
            else:
                self.log.warning(                                                # pylint: disable=logging-fstring-interpolation
                    f'No .pdb file for HLA {hla} not found', extra={'messagePrefix': called_by}
                )
                continue
        return vis_data

    def _get_min_hlavsep(self, epitopes:set) -> dict:
        if not isinstance(epitopes, set):
            raise TypeError('Epitopes should be given as a set')
        return self.Epitope.min_hlavsep(epitopes)

    def _hlavsep_poly(
                self,
                hlavsep:dict,
                rAb:bool=False,
                mAb:bool=False,
                ) -> dict:
        """ An internal method to get hla vs polymorphic epitopes dictionary from
        hla vs epitopes dictionary  """

        ep_db = self.Epitope.df
        _hlavsep = defaultdict(lambda: defaultdict(list))
        for hla in hlavsep.keys():
            for ep in hlavsep[hla]:
                _hlavsep[hla]['desa'].append(ep)
                _hlavsep[hla]['_desa'].extend(polymorphic_residues(ep, ep_db))
                try:
                    ind = ep_db.Epitope == ep
                    if rAb & (ep_db[ind]['AntibodyReactivity'].values[0] == 'Yes'):
                        _hlavsep[hla]['desa_rAb'].append(ep)
                        _hlavsep[hla]['_desa_rAb'].extend(polymorphic_residues(ep, ep_db))
                    if mAb & (ep_db[ind]['mAb'].values[0] == 'Yes'):
                        _hlavsep[hla]['desa_mAb'].append(ep)
                        _hlavsep[hla]['_desa_mAb'].extend(polymorphic_residues(ep, ep_db))
                except IndexError:
                    self.log.warning(
                        f'Epitope {ep} from HLA {hla} does not exist in the epitope database',
                        extra={'messagePrefix': 'Polymorphic resisdues'}
                    )
        return _hlavsep

    def from_epitopes(self,
                    epitopes:set,
                    style:str='sphere',
                    rAb:bool=False,
                    mAb:bool=False,
        ) -> Dict[str, set]:
        """
        Prepare 3D visualisation data from epitopes as input. This method finds the
        minimum number of hla on which all the epitopes can be located and return a dict
        {'hla': {'ep'} }
        """

        self.log.info(f'Start with Epitopes {epitopes}', extra={'messagePrefix': 'from_epitopes'})
        _min_hlavsep = self._get_min_hlavsep(epitopes)
        self.log.info(f'min HLA vs Epitope is:{_min_hlavsep.keys()}', extra={'messagePrefix': 'from_epitopes'})
        hlavspolyep = self._hlavsep_poly(_min_hlavsep, rAb, mAb)
        vis_data = self._get_vis_data_from_pdb(hlavspolyep, style, 'from_epitopes')
        self.vis_data_from_epitopes = vis_data
        self.hlavsep = _min_hlavsep
        return self
    
    def from_transplant(self,
                        TxIDs:Set[int],
                        style:str='sphere',
                        rAb:bool=False,
                        mAb:bool=False,
        ) -> dict:
        """
        Prepare 3D visualisation data from transplants as input. This method return a dict
        {
            'TxID': {'hla': {'ep'} }
        }
        """

        self.log.info(f'Start with TxID {TxIDs}', extra={'messagePrefix': 'from_transplant'})
        if not isinstance(TxIDs, set):
            raise TypeError('Transplants should be given as a set')
        vis_data = defaultdict(dict)
        txvshlavsep = defaultdict(dict)

        for TxID in TxIDs:
            df_tx = self.DESA.get_tx(TxID)
            epvshla = df_tx.EpvsHLA_Donor.values[0]
            hlavsep = self.Epitope.epvshla2hlavsep(epvshla)  # call a method in Epitope class
            txvshlavsep[TxID] = hlavsep
            hlavspolyep = self._hlavsep_poly(hlavsep, rAb, mAb)
            vis_data[TxID] = self._get_vis_data_from_pdb(hlavspolyep, style, 'from_transplants')

        self.vis_data_from_transplants = vis_data
        self.txvshlavsep = txvshlavsep
        return self

    def visualise3D(self, hla:str, TxID:int=None):
        """ This function passes model and style data to the front-end viewer """
        self.log.info(
            'Preparing the visualisation payload', extra={'messagePrefix': f'TXID:{TxID}, HLA:{hla}'}
        )
        if self.hlavsep: # this attribute is updated by from_epitopes method
            model_data = json.loads(self.vis_data_from_epitopes.get(hla).get('model'))
            style_data = json.loads(self.vis_data_from_epitopes.get(hla).get('style'))

        elif self.txvshlavsep: # this attribute is updated by from_transplants method
            model_data = json.loads(self.vis_data_from_transplants.get(TxID).get(hla).get('model'))
            style_data = json.loads(self.vis_data_from_transplants.get(TxID).get(hla).get('style'))
        return self._molecule_viewer(model_data, style_data, opacity=0.6)

####################################################
# Front-end visulaisation healper functions
####################################################

def vis_cards(vis_object):
    """
    This is the visualisation method used in the app. It uses the
    vis_object from VisualiseHLA and based on the object nature
    from_transplants or from_epitopes it calls vis_card method.
    In the later case, each transplant belongs to one card and in
    case of the former, the all the hla's are visualised in a card.
    """

    if vis_object.txvshlavsep: # Object is carrying vis data from transplants
        return [_vis_card(vis_object, TxID) for TxID in vis_object.txvshlavsep.keys() ]
    if vis_object.hlavsep:
        return _vis_card(vis_object)


def _vis_card(vis_object, TxID:int=None):
    """ Front-end card embelishing the visualisation payload in one card"""

    if TxID:
        card_header = dbc.CardHeader(
                                    html.H4(f""" Transplant ID: {TxID},
                                            Survival Time [Y]: {vis_object.DESA.get_tx(TxID)['Survival[Y]'].values[0]: .2f}
                                            """
                                            , className="card-title"
                                    )
                        )
        visualisation_payload = [ dbc.Col(vis_payload(vis_object, i, hla, TxID), width=6)
                                                for i, hla in enumerate(vis_object.txvshlavsep[TxID].keys()) ]
    else:
        # If visualisation is from_epitopes there is no header
        card_header = dbc.CardHeader( html.H4(''), className="card-title")
        visualisation_payload = [ dbc.Col(vis_payload(vis_object, i, hla), width=6)
                                                for i, hla in enumerate(vis_object.hlavsep.keys()) ]
    card_body =  dbc.CardBody(dbc.Row(visualisation_payload))
    return dbc.Card(
                [card_header, card_body]
                , color='secondary', style={'padding': 10}
            )

def vis_payload(vis_object, i, hla, TxID:int=None):
    """
    This methods embelishes and wraps up the visualise3D method in
    VisualiseHLA class and deliver it to v-s_card method
    """

    hlavsep = vis_object.txvshlavsep[TxID] if TxID else vis_object.hlavsep
    return [
            html.H6(
                    html.Span(f'{hla}',
                        id=f"tooltip-target-{TxID}-{i}",
                        style={"textDecoration": "underline", "cursor": "pointer"},
                        className="card-subtitle",
                    )
            ),
            dbc.Tooltip(
                f"DESA #{len(hlavsep[hla])}: {hlavsep[hla]}",
                target=f"tooltip-target-{TxID}-{i}",
                placement='top'
            ),
            vis_object.visualise3D(hla, TxID)
        ]
