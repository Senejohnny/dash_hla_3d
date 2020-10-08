import json
import tempfile
import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

import six.moves.urllib.request as urlreq

import dash
from dash import no_update
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State

import dash_daq as daq
import dash_table
from app.analysis_utils import load_data, data_3dviewer, div_3dviewer

from app.dash_utils import parse_contents, upload_file
from pygit2 import Repository


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, title='HLA 3D Viewer', external_stylesheets=external_stylesheets)

# app.css.append_css({
#     "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
# })

epitope_db, desa_db = load_data()

desa_db = desa_db[desa_db.Status=='DESA']
Tx = 4606

hlas, _3d_data, _hlavsdesa, hlavsdesa = data_3dviewer(desa_db, epitope_db, Tx)

hla_1, component_1 = div_3dviewer(hlas, 0, _3d_data, _hlavsdesa)
hla_2, component_2 = div_3dviewer(hlas, 1, _3d_data, _hlavsdesa)
hla_3, component_3 = div_3dviewer(hlas, 2, _3d_data, _hlavsdesa)
hla_4, component_4 = div_3dviewer(hlas, 3, _3d_data, _hlavsdesa)
hla_5, component_5 = div_3dviewer(hlas, 4, _3d_data, _hlavsdesa)
hla_6, component_6 = div_3dviewer(hlas, 5, _3d_data, _hlavsdesa)


desa_db_table = desa_db[['TransplantID', 'Status', '#DESA', 'FailureCode_NOTR', 'GraftSurvivalYears_NOTR']]
desa_db_table['DESA->Donor_HLA'] = desa_db['DESA->Donor_HLA'].apply(lambda x: str(dict(x)))
desa_db_table = desa_db_table.sort_values(by='#DESA', ascending=False)
# ######################################################################
# App Layout
# ######################################################################
app.layout = html.Div([
    html.Div([
        html.Div(
            html.Div(
                children=[
                html.H2('HLA 3D Viewer: V 0.0.2'),

            ]),
            className ="nine columns",
        ),
        html.Div(
            html.Img(
                src="assets/logo.svg",
                className ="logo-header",
                style={
                    'height': '12%',
                    'width': '12%', 
                    'float': 'right',
                    'position': 'relative',
                    'padding-top': 2, 
                    'padding-right': 4,
                }
            )
        )
    ], className="row"),
    html.Div([
        dcc.Tabs(id='tabs', value='tab-1', children=[
            dcc.Tab(
                label='Transplants',
                className='control-upload',
                value='tab-1',
                children=[
                    html.Div([
                        html.Div([
                            html.H5('Upload DESA Analysis [.csv & .xlsx]'),
                            html.Div(
                                upload_file('DESA', id='upload-desa-data'),
                                style={'padding': 0, 'columnCount': 1},
                            ),
                            html.Div(id='output-desa-upload'),
                        ]),
                    html.H5("Choose Tx for 3d View:"),
                    dcc.Input(
                        id="transplant_3d_view", 
                        placeholder="Enter Tx number'",
                        type='text',
                        value=''
                    )
                    ], className="four columns"),
                    html.H1(""),
                    html.Div(id='output', 
                            children=dash_table.DataTable(
                                id='table',
                                columns=[{"name": i, "id": i} for i in desa_db_table.columns],   
                                data=desa_db_table.to_dict('records'),
                                page_size= 20,
                                editable=True,
                                style_table={'height': '400px', 
                                             'overflowY': 'auto'},
                                style_cell_conditional=[
                                                        {
                                                            'if': {'column_id': c},
                                                            'textAlign': 'left'
                                                        } for c in ['Date', 'Region']
                                                    ],
                                style_data_conditional=[
                                                        {
                                                            'if': {'row_index': 'odd'},
                                                            'backgroundColor': 'rgb(248, 248, 248)'
                                                        }
                                                    ],
                                style_header={
                                                'backgroundColor': 'rgb(230, 230, 230)',
                                                'fontWeight': 'bold'
                                            }
                            ), className="eight columns"),  # end of table
            ]), 
            dcc.Tab(label='3D View', 
                    value='tab-2',
                    children=html.Div([
                        html.H5('Overview'),
                        html.Div(f'''
                            This is a 3D view of HLA molecules in Transplant ID: {Tx}, with DESA color depicted in yellow  
                        '''),
                        dcc.Loading(
                            type="circle",
                            children=[
                            html.Div(      # HLA 1
                                id='mol3d-viewer-1',
                                children=[
                                html.H5(f'HLA: {hla_1}'),
                                # html.P(f'DESA: {hlavsdesa_short[hla_1]}'),
                                component_1,
                            ], className="four columns"),
                            html.Div(      # HLA 2
                                id='mol3d-viewer-2',
                                children=[
                                html.H5(f'HLA: {hla_2}'),
                                # html.P(f'DESA: {hlavsdesa_short[hla_2]}'),
                                component_2,
                            ], className="four columns"),
                            html.Div(      # HLA 3
                                id='mol3d-viewer-3',
                                children=[
                                html.H5(f'HLA: {hla_3}'),
                                # html.P(f'DESA: {hlavsdesa_short[hla_3]}'),
                                component_3
                            ], className="four columns"),
                            html.Div(      # HLA 4
                                id='mol3d-viewer-4',
                                children=[
                                html.H5(f'HLA: {hla_4}'),
                                component_4
                            ], className="four columns"),
                            html.Div(      # HLA 5
                                id='mol3d-viewer-5',
                                children=[
                                html.H5(f'HLA: {hla_5}'),
                                component_5
                            ], className="four columns"),
                            html.Div(      # HLA 6
                                id='mol3d-viewer-6',
                                children=[
                                html.H4(f'HLA: {hla_6}'),
                                None
                            ], className="four columns"),
                        ])
                    ]),
            ),
            dcc.Tab(
                label='3D Replication',
                className='control-upload',
                children=[
                    html.Div([
                        html.H5("Which HLA should be repeated?"),
                        dcc.Input(
                            id="rep-hla", 
                            placeholder='Enter a HLA',
                            type='text',
                            value=''
                        ), 
                        html.H5("Which DESA Should be marked?"),
                        dcc.Input(
                            id="rep-desa", 
                            placeholder="Enter DESA's",
                            type='text',
                            value=''
                        )
                    ], className= "three columns"),
                    html.Div(
                        html.H4("The 3D HLA structure would be represented below"),
                        className= "nine columns",
                    )
                ]
            )
        ]),
        # Hidden div inside the app that stores the intermediate value
        html.Div(id='text-output', children='')#, style={'display': 'none'}),
    ]),
],)
# ######################################################################
# App Layout Pieces
# ######################################################################



# ######################################################################
# Callbacks
# ######################################################################

# 
# @app.callback([Output('hidden-desa-data', 'children'),
#             Output('output-desa-upload', 'children'),], 
#             [Input('upload-desa-data', 'contents'),],
#             [State('upload-desa-data', 'filename'),])
# def update_output_load_data(contents, filename):
#     if contents is not None:
#         df_json, children = parse_contents(contents, filename, 'DESA')
#         return df_json, children
#     else:
#         return no_update,  html.Div(['No data is uploaded'])

@app.callback(Output('text-output', 'children'),
            [Input('transplant_3d_view', 'value')],)
def get_active_cell_value (value):
    if value:
        return value
    else:
        return no_update

def _output(hla, component):
    return html.H4(f'HLA: {hla_6}'), component

# @app.callback(Output('mol3d-viewer-6', 'children'),
#             [Input('transplant_3d_view', 'value')],)
# def get_tx_3d(Tx):
#     print(Tx)
#     if Tx:
#         print('hi')
#         hlas, _3d_data, _hlavsdesa, hlavsdesa = data_3dviewer(desa_db, epitope_db, Tx)
#         print(Tx, hlavsdesa)
#         hla_1, component_1 = div_3dviewer(hlas, 0, _3d_data, _hlavsdesa)
#         # hla_2, component_2 = div_3dviewer(hlas, 1, _3d_data, _hlavsdesa)
#         # hla_3, component_3 = div_3dviewer(hlas, 2, _3d_data, _hlavsdesa)
#         # hla_4, component_4 = div_3dviewer(hlas, 3, _3d_data, _hlavsdesa)
#         # hla_5, component_5 = div_3dviewer(hlas, 4, _3d_data, _hlavsdesa)
#         # hla_6, component_6 = div_3dviewer(hlas, 5, _3d_data, _hlavsdesa)
#         return _output(hla_1, component_1)
#     else:
#         return no_update

if __name__ == '__main__':
    if Repository('.').head.shorthand == 'master':
        debug = False
    else:
        debug = True
    app.run_server(debug=debug, host = '0.0.0.0', port=8080)
