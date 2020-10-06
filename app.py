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
from app.analysis_utils import load_data, data_3dviewer, div_3dviewer

from app.dash_utils import parse_contents, upload_file
from pygit2 import Repository




app = dash.Dash(__name__, title='HLA 3D Viewer')

app.css.append_css({
    "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
})

epitope_db, desa_db = load_data()
hlas, _3d_data, _hlavsdesa, hlavsdesa = data_3dviewer(desa_db, epitope_db, 4080)

hla_1, component_1 = div_3dviewer(hlas, 0, _3d_data, _hlavsdesa)
hla_2, component_2 = div_3dviewer(hlas, 1, _3d_data, _hlavsdesa)
hla_3, component_3 = div_3dviewer(hlas, 2, _3d_data, _hlavsdesa)
hla_4, component_4 = div_3dviewer(hlas, 3, _3d_data, _hlavsdesa)
hla_5, component_5 = div_3dviewer(hlas, 4, _3d_data, _hlavsdesa)
hla_6, component_6 = div_3dviewer(hlas, 5, _3d_data, _hlavsdesa)


# ######################################################################
# App Layout
# ######################################################################
app.layout = html.Div([
    html.Div([
        html.Div(
            html.Div(
                className="app-header", 
                children=
                html.Div('HLA 3D View: V 0.0.1', className="app-header--title"),
            ),
            className ="nine columns",
        ),
        html.Div(

            className ="nine columns",
        )],
    className="row"),
    html.Div([
        dcc.Tabs(id='tabs', value='tab-1', children=[
            dcc.Tab(label='Transplants',
                    className='control-upload',
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
                        html.Div(      # HLA 1
                            id='mol3d-viewer-1',
                            children=[
                            html.P(f'HLA: {hla_1}'),
                            # html.P(f'DESA: {hlavsdesa_short[hla_1]}'),
                            component_1,
                        ], className="four columns"),
                        html.Div(      # HLA 2
                            id='mol3d-viewer-2',
                            children=[
                            html.P(f'HLA: {hla_2}'),
                            # html.P(f'DESA: {hlavsdesa_short[hla_2]}'),
                            component_2,
                        ], className="four columns"),
                        html.Div(      # HLA 3
                            id='mol3d-viewer-3',
                            children=[
                            html.P(f'HLA: {hla_3}'),
                            # html.P(f'DESA: {hlavsdesa_short[hla_3]}'),
                            component_3
                        ], className="four columns"),
                        html.Div(      # HLA 4
                            id='mol3d-viewer-4',
                            children=[
                            html.P(f'HLA: {hla_4}'),
                            component_4
                        ], className="four columns"),
                        html.Div(      # HLA 5
                            id='mol3d-viewer-5',
                            children=[
                            html.P(f'HLA: {hla_5}'),
                            component_5
                        ], className="four columns"),
                        html.Div(      # HLA 6
                            id='mol3d-viewer-6',
                            children=[
                            html.P(f'HLA: {hla_6}'),
                            component_6
                        ], className="four columns"),
                    ]),
            ),
        ]),
        # Hidden div inside the app that stores the intermediate value
        html.Div(id='hidden-desa-data', style={'display': 'none'}),
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


if __name__ == '__main__':
    if Repository('.').head.shorthand == 'master':
        debug = False
    else:
        debug = True
    app.run_server(debug=debug, host = '0.0.0.0', port=8080)
