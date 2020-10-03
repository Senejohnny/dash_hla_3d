import json
import tempfile
import os
import pandas as pd
import six.moves.urllib.request as urlreq

import dash
from dash import no_update
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bio as dashbio
import dash_daq as daq
from dash_bio_utils import pdb_parser as parser
from app import styles_parser as sparser

from app.analysis_utils import polymorphic_residues, hla_to_filename, find_molecule_path
from app.dash_utils import parse_contents, upload_file
from pygit2 import Repository

# Load Epitope Database
ep_path = os.path.expanduser('./data/20200804_EpitopevsHLA_distance.pickle')
epitope_db = pd.read_pickle(ep_path)

des_path = '/Users/Danial/UMCUtrecht/ProcessedData/DSAandDESA/20200916_DESA_new_All.pickle'
desa_db = pd.read_pickle(des_path)
print(desa_db)

#EpvsHLA
epvshla_donor = {'44RMA': 'B*57:01',
                '62GE': 'B*57:01',
                '94I': 'B*57:01',
                '71SA': 'B*57:01',
                '74Y': 'B*57:01',
                '62GRN': 'B*57:01',
                '97V': 'B*57:01'} 

from collections import defaultdict
reversed_dict = defaultdict(list)
for key, value in epvshla_donor.items():
    reversed_dict[value].append(key)

# print({hla:[].append(ep) for ep, hla in epvshla_donor.items()})
poly_res = polymorphic_residues(set(epvshla_donor.keys()), epitope_db)
# print(poly_res)
# Function to create the modelData and style files for molecule visualization
# hla = set(epvshla_donor.values()) 
# print(list(hla))
locus, filename = hla_to_filename('B*57:01')
pdb_exist, pdb_path = find_molecule_path(locus, filename)



app = dash.Dash(__name__, title='HLA 3D Viewer')

# Create the model data from the decoded contents
if pdb_exist:
    model = parser.create_data(pdb_path)
else:
    print(f'{pdb_path}: IOError: No such file or directory')
    pass 

desa_info = {'chain': 'A',
             'relevant_desa': set([int(_[0]) for _ in poly_res]),
             'irrelevant_desa': set([101, 102, 103, 104, 105]),
            }
# print(desa_info)
# Create the cartoon style from the decoded contents
style = sparser.create_style(pdb_path, style='sphere', mol_color='chain', desa_info=desa_info)

model = json.loads(model)
style =  json.loads(style)

component = dashbio.Molecule3dViewer(
                                    id='mol-3d',
                                    selectionType='atom ',
                                    modelData=model,
                                    styles=style,
                                    selectedAtomIds=[],
                                    backgroundOpacity='0',
                                    atomLabelsShown=False,
                                    )

# ######################################################################
# App Layout
# ######################################################################
app.layout = html.Div([
    html.Div(
        className="app-header",
        children=[
            html.Div('HLA 3D view: V_0.0.1', className="app-header--title")
        ]
    ),
    
    html.Div([
    dcc.Tabs(id='tabs', value='tab-1', children=[
        dcc.Tab(label='Transplants', 
                value='tab-1',
                children=html.Div([
                html.H5(['Upload DESA Analysis [.csv & .xlsx]']),
                    html.Div([
                        upload_file('DESA', id='upload-desa-data'),
                    ], style={'padding': 0, 'columnCount': 1},
                    ),
                html.Div(id='output-desa-upload'),
                ]),
        ), 
        dcc.Tab(label='3D View', 
                value='tab-2',
                children=html.Div([
                    html.H5('Overview'),
                    html.Div('''
                        This is a 3D view of an HLA with DESA color coded  
                    '''),
                    html.Div(
                        component
                    ),
                ]),
        ),
    ]),
    # Hidden div inside the app that stores the intermediate value
    html.Div(id='hidden-desa-data', style={'display': 'none'}),
    ]),   

])
# ######################################################################
# App Layout Pieces
# ######################################################################



# ######################################################################
# Callbacks
# ######################################################################

# 
@app.callback([Output('hidden-desa-data', 'children'),
            Output('output-desa-upload', 'children'),], 
            [Input('upload-desa-data', 'contents'),],
            [State('upload-desa-data', 'filename'),])
def update_output_load_data(contents, filename):
    if contents is not None:
        df_json, children = parse_contents(contents, filename, 'DESA')
        return df_json, children
    else:
        return no_update,  html.Div(['No data is uploaded'])


if __name__ == '__main__':
    if Repository('.').head.shorthand == 'master':
        debug = False
    else:
        debug = True
    app.run_server(debug=debug, host = '0.0.0.0', port=8080)