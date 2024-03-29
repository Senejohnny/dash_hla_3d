""" This file encompass all the 3D visualisation scripts"""
import json
import logging
from typing import Set, Dict, Union
from collections import defaultdict
from dash_bio import Molecule3dViewer
import dash_bootstrap_components as dbc
import dash_html_components as html
from app.common.logger import get_logger

from app.epitope import Epitope
from app.desa import DESA
from app import styles_parser as sparser, pdb_parser as parser
from app.common.utilities import (
    get_hla_exp,
    get_hla_polychain,
    hla_to_filename,
    find_molecule_path,
)

SERVICE_NAME = 'VisualiseHLA'

class VisualiseHLA:

    def __init__(self, ignore_hla:set=set(),  path_desa=None, path_epitope=None):
        """ uses DESA and Epitope Classes """
        self.ignore_hla = ignore_hla
        self.desa = DESA(path_desa) if path_desa else DESA()
        self.epitope = Epitope(path_epitope) if path_epitope else Epitope()
        self.log = get_logger(SERVICE_NAME, logging.INFO)
        self.hlavsep = None
        self.txvshlavsep = None
        self.vis_data_from_transplants = None


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
            # try:
            #     desa_info['desa_rAb'] = set( map(lambda x: int(x[0]), hlavsep[hla]['_desa_rAb']) )
            # except KeyError:
            #     desa_info['desa_rAb'] = set()
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

    def _get_min_hlavsep(self, epitopes:set, most_freq_hla:bool=True) -> dict:
        if not isinstance(epitopes, set):
            raise TypeError('Epitopes should be given as a set')
        return self.epitope.min_hlavsep(epitopes, ignore_hla=self.ignore_hla, most_freq_hla=most_freq_hla)

    def _keep_hla_with_pdb(self, hlavsep:dict) -> dict:
        """
        Keeps the HLA's that a pdb file is available in the inventory.
        Gets and returns HLA vs Epitope dictionary
        """
        for hla in list(hlavsep):
            if hla not in self.epitope.pdb_inventory:
                self.log.warning(
                    f'HLA {hla} is removed: No relevant pdb file found',
                    extra={'messagePrefix': 'Filtering HLA not in inventory'}
                )
                del hlavsep[hla]
        return hlavsep

    def _hlavsep_poly(
                self,
                hlavsep:dict,
                mAb:bool=False,
                ) -> dict:
        """ An internal method to get hla vs polymorphic epitopes dictionary from
        hla vs epitopes dictionary  """

        ep_db = self.epitope.df
        _hlavsep = defaultdict(lambda: defaultdict(list))
        for hla in hlavsep.keys():
            for ep in hlavsep[hla]:
                _hlavsep[hla]['desa'].append(ep)
                _hlavsep[hla]['_desa'].extend(self.epitope.polymorphic_residues(ep))
                try:
                    ind = ep_db.Epitope == ep
                    if mAb & (ep_db[ind]['isotype'].values[0] == 'IgG'):
                        _hlavsep[hla]['desa_mAb'].append(ep)
                        _hlavsep[hla]['_desa_mAb'].extend(self.epitope.polymorphic_residues(ep))
                except IndexError:
                    self.log.warning(
                        f'Epitope {ep} from HLA {hla} does not exist in the epitope database',
                        extra={'messagePrefix': 'Polymorphic resisdues'}
                    )
        return _hlavsep

    def from_epitopes(self,
                    epitopes:set,
                    style:str='sphere',
                    mAb:bool=False,
                    elliproscore:Union[set, list, str]=None,
        ) -> Dict[str, set]:
        """
        Prepare 3D visualisation data from epitopes as input. This method finds the
        minimum number of hla on which all the epitopes can be located and return a dict
        {'hla': {'ep'} }
        """

        self.log.info(f'Start with Epitopes {epitopes}', extra={'messagePrefix': 'from_epitopes'})
        if not isinstance(epitopes, set):
            raise TypeError('Epitopes should be given as a set')

        # filter epitope df (under the hood) by ellipro score
        if elliproscore:
            self.epitope.ellipro(elliproscore)

        _min_hlavsep = self._get_min_hlavsep(epitopes)
        self.log.info(f'min HLA vs Epitope is:{list(_min_hlavsep.keys())}', extra={'messagePrefix': 'from_epitopes'})
        hlavspolyep = self._hlavsep_poly(_min_hlavsep, mAb)
        vis_data = self._get_vis_data_from_pdb(hlavspolyep, style, 'from_epitopes')
        self.vis_data_from_epitopes = vis_data
        self.hlavsep = _min_hlavsep
        return self

    def from_transplant(self,
                        TxIDs:Set[int],
                        style:str='sphere',
                        mAb:bool=False,
                        elliproscore:Union[set, list, str]=None,
        ) -> dict:
        """
        Prepare 3D visualisation data from transplants as input. This method return a dict
        {
            'TxID': {'hla': {'ep'} }
        }
        """
        self.log.info(f'Start with TxID {TxIDs}', extra={'messagePrefix': 'from_transplant'})
         # filter epitope df (under the hood) by ellipro score
        if elliproscore:
            self.epitope.ellipro(elliproscore)


        if not isinstance(TxIDs, set):
            raise TypeError('Transplants should be given as a set')

        vis_data = defaultdict(dict)
        txvshlavsep = defaultdict(dict)

        for TxID in TxIDs:
            epvshla = self.desa.get_tx(TxID).EpvsHLA_Donor.values[0]
            # Apply ellipro filter to epitopes
            epvshla = {ep:epvshla[ep] for ep in epvshla.keys()
                                        if ep in set(self.epitope.df.Epitope.values.tolist()) }
            hlavsep = self.epitope.epvshla2hlavsep(epvshla)
            # Keep HLA's for which a pdb file exist
            hlavsep = self._keep_hla_with_pdb(hlavsep)
            txvshlavsep[TxID] = hlavsep
            hlavspolyep = self._hlavsep_poly(hlavsep, mAb)
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
# Front-end visualisation helper functions
####################################################

def vis_cards(vis_object: VisualiseHLA):
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


def _vis_card(vis_object: VisualiseHLA, TxID:int=None):
    """ Front-end card embellishing the visualisation payload in one card """

    if TxID:
        card_header = dbc.CardHeader(
                                    html.H4(f""" Transplant ID: {TxID},
                                            Survival Time [Y]: {vis_object.desa.get_tx(TxID)['Survival[Y]'].values[0]: .2f}
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

def vis_payload(vis_object: VisualiseHLA, i, hla, TxID:int=None):
    """
    This methods embellishes and wraps up the visualise3D method in
    VisualiseHLA class and deliver it to vis_card method
    """
    hlavsep = vis_object.txvshlavsep[TxID] if TxID else vis_object.hlavsep
    return [
            html.H6(
                    [
                        html.Span(f'{hla}',
                            id=f"tooltip-target-{TxID}-{i}",
                            style={"textDecoration": "underline", "cursor": "pointer"},
                            className="card-subtitle",
                        ),
                        # f" [Exp:{get_hla_exp(hla)}]"
                    ]
            ),
            dbc.Tooltip(
                f"DESA #{len(hlavsep[hla])}: {hlavsep[hla]}",
                target=f"tooltip-target-{TxID}-{i}",
                placement='top'
            ),
            vis_object.visualise3D(hla, TxID)
        ]
