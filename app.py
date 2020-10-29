import json
import tempfile
import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

import six.moves.urllib.request as urlreq

import dash
from dash import no_update
import dash_daq as daq
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State


from app.analysis_utils import load_desa_db, data_3dviewer, div_3dviewer, load_epitope_db

from app.dash_utils import filtering_logic, logo_img, dashtable_data_compatibility
from pygit2 import Repository


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, title='HLA 3D Viewer', external_stylesheets=external_stylesheets)
app.config['suppress_callback_exceptions'] = True

# app.css.append_css({
#     "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
# })




# ######################################################################
# App Layout
# ######################################################################
app.layout = html.Div([
    html.Div([
        html.Div(
            html.Div(
                children=[
                html.H2('HLA 3D Viewer: V 0.0.3'),
            ]),
            className ="nine columns",
        ),
        html.Div(
            logo_img()
        )
    ], className="row"),
    html.Div([
        dcc.Tabs(id='tabs-mother', value='tab-1', children=[
            dcc.Tab(
                label='Transplants',
                className='control-upload',
                value='tab-1',
                children=[
                    html.Div([
                        dcc.Tabs(id='tabs-children', value='tab-1-1', 
                            children=[
                                dcc.Tab(label='About', value='tab-1-1'),
                                dcc.Tab(label='Data', value='tab-1-2'),
                                dcc.Tab(label='3D View', value='tab-1-3'),
                        ], className="four columns"),
                        html.Div(id='tabs-children-content')
                    ], style={'padding': 0 }),
                    html.Div(id='tx-table',
                            children=[],
                            className="eight columns"),
            ]), 
            dcc.Tab(label='3D View', 
                    value='tab-2',
                    children=html.Div([
                        html.H5('HLA Molecule 3D View'),
                        html.Div(f'''
                            This is a 3D view of HLA molecules in Transplant ID: , with DESA color depicted in yellow  
                        '''),
                        dcc.Loading(
                            id='3d-view-loading', 
                            type="circle",
                        )
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
        html.Div(id='hidden-df', style={'display': 'none'}),
        html.Div(id='hidden-data', style={'display': 'none'}),
    ]),
],)
# ######################################################################
# App Layout Pieces
# ######################################################################  
# tab_2_layout = html.Div([
#                     html.Div([
#                         html.H5('Upload DESA Analysis [.csv & .xlsx & .pkl]'),
#                         html.Div(
#                             upload_file('DESA', id='upload-desa-data'),
#                             style={'padding': 0, 'columnCount': 1},
#                         ),
#                         html.Div(id='output-desa-upload'),
#                     ]),
#                 ], className="four columns"),
# table_layout = 
tab_2_layout = html.Div([
                        html.H3('Loading & Filtering:'),
                        html.Div([
                        html.Button('Show Table', id='show-table', n_clicks=0),
                        html.H6('Sort By'),
                        html.Div(
                            dcc.Dropdown(
                                id= 'dropdown_sortby',
                                options=[
                                    {'label': 'Early Failure', 'value': 'early_failure'},
                                    {'label': 'Late Surviving', 'value': 'late_surviving'},
                                    {'label': 'DESA', 'value': 'desa'}
                                ],
                                placeholder="Select HLA Class",
                            ),
                        )
                    ]), 
                    html.Div([
                        html.H6('HLA Class'),
                        html.Div(
                            dcc.Dropdown(
                                id= 'dropdown_class',
                                options=[
                                    {'label': 'I', 'value': 'I'},
                                    {'label': 'II', 'value': 'II'},
                                    {'label': 'I & II', 'value': 'III'}
                                ],
                                placeholder="Select HLA Class",
                            ),
                        )
                    ]), 
                    html.Div([
                        html.H6('HLA Molecule'),
                        html.Div(dcc.Input(id='input_hla', type='text', placeholder="HLA")),
                    ]),
                ], className="four columns", style={'padding': 10}),

tab_3_layout = html.Div([
                    html.Div([
                            daq.BooleanSwitch(
                                id='rAb-switch',
                                on=False,
                                label="Monoclonal Antobodies",
                                labelPosition="top",
                            ),
                    ], style={'padding': 10}),
                    html.Div([
                        html.H6('Transplant ID'),
                        html.Div(dcc.Input(id='input-tx', type='text', placeholder="Tx ID")),
                        html.H6('Style'),
                        dcc.Dropdown(
                                id= 'dropdown_style',
                                options=[
                                    {'label': 'Stick', 'value': 'stick'},
                                    {'label': 'Cartoon', 'value': 'cartoon'},
                                    {'label': 'Sphere', 'value': 'sphere'}
                                ],
                                placeholder="Select HLA Class",
                        ),
                        html.Div(
                            html.Button('Show', id='submit-tx-show', n_clicks=0),
                        style={'padding': 10})
                    ], style={'padding': 10})
                ], 
                className="four columns",)


                
# ######################################################################
# Callbacks
# ######################################################################


@app.callback(Output('tabs-children-content', 'children'),
              [Input('tabs-children', 'value')])
def render_content(tab):
    if tab == 'tab-1-1':
        return html.Div([
            html.H3('What is 3D HLA Viewer'),
            html.P("""
            Molecule3D is a visualizer that allows you to view biomolecules in multiple representations: 
            sticks, spheres, and cartoons.You can select a preloaded structure, or upload your own, in the
            "Data" tab. A sample structure is also available to download. In the "View" tab, you can change 
            the style and coloring of the various components of your molecule.
            """)
        ], className="four columns", style={'padding': 10})
    elif tab == 'tab-1-2':
        return tab_2_layout
    elif tab == 'tab-1-3':
        return tab_3_layout



@app.callback(Output('tx-table', 'children'),
              [Input('show-table','n_clicks'),],
              [State('dropdown_sortby', 'value'),
               State('dropdown_class', 'value'),
               State('input_hla', 'value')])
def update_output_data(n_clicks, sort_failure, sort_class, hla):
    if n_clicks == 0:
        return html.P('Click on Show button to see the table')
    else:
        desa_db = load_desa_db()
        desa_db = filtering_logic(desa_db, sort_failure, sort_class, hla)
        df_print = dashtable_data_compatibility(desa_db)
        return dash_table.DataTable(
                                columns=[{"name": i, "id": i} for i in df_print.columns],   
                                data=df_print.to_dict('records'),
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
                        )         


@app.callback(Output('3d-view-loading', 'children'),
              Input('submit-tx-show','n_clicks'),
              [State('input-tx', 'value'),
               State('dropdown_style', 'value'),
               State('rAb-switch', 'on')])
def show_3d_tx(n_clicks, TxID, style, rAb_switch):
    if n_clicks:
        desa_df = load_desa_db()
        epitope_db = load_epitope_db()
        _3d_data, hlavsdesa = data_3dviewer(desa_df, epitope_db, int(TxID), style=style, rAb=rAb_switch)
        return [html.Div([
                    html.H4(f'Transplant ID: {TxID}'),
                    html.Div([
                        html.Div([
                            html.Div([
                                html.H6([
                                        html.Span(f'HLA: {hla}',
                                            id=f"tooltip-target-{i}",
                                            style={"textDecoration": "underline", "cursor": "pointer"},
                                        )
                                ]),
                                dbc.Tooltip(
                                    f"DESA #{len(hlavsdesa[hla]['desa'])}: {hlavsdesa[hla]['desa']}",
                                    # "DESA rAb: {hlavsdesa[hla].get('desa_rAb')}"
                                    target=f"tooltip-target-{i}",
                                    placement='top'
                                )
                            ]),
                            div_3dviewer(hla, _3d_data)], 
                            className="four columns") for i, hla in enumerate(hlavsdesa.keys())
                        ])
                    ])
                ]
    else:
        return 'No Transplant ID is selected for visualisation'



if __name__ == '__main__':
    if Repository('.').head.shorthand == 'master':
        debug = False
    else:
        debug = True
    app.run_server(debug=debug, host = '0.0.0.0', port=8080)

# Error 1272
