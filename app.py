import json
import tempfile
import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# import pymongo
# from pymongo import MongoClient

import dash
from dash import no_update
import dash_daq as daq
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from flask_caching import Cache
from app.epitope import Epitope_DB

from app.analysis_utils import (
    load_desa_db, 
    data_3dviewer, 
    div_3dviewer,
    load_epitope_db
    )

from app.dash_utils import filtering_logic, Header, dashtable_data_compatibility
from pygit2 import Repository


# mongo_client = MongoClient('localhost', 27017) # build a new client instance of MongoClient
# db = mongo_client.desa_database # create new database
# desa_col = db['desa_db'] # create new collection instance
# desa_df = load_desa_db()
# data_dict = desa_df.to_dict("records") # convert to dictionary
# desa_col.insertMany({"index":"desa_db","data":data_dict}) # inesrt into DB

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, title='HLA-Epitope 3D', external_stylesheets=[dbc.themes.CERULEAN])
app.config['suppress_callback_exceptions'] = True

# app.css.append_css({
#     "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
# })

cache = Cache(app.server, config={
    'CACHE_TYPE': 'redis',
    # Note that filesystem cache doesn't work on systems with ephemeral
    # filesystems like Heroku.
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': 'cache-directory',

    # should be equal to maximum number of users on the app at a single time
    # higher numbers will store more data in the filesystem / redis cache
    'CACHE_THRESHOLD': 200
})




# ######################################################################
# App Layout
# ######################################################################
app.layout = dbc.Container([
    Header('HLA-Epitope 3D'),
    html.Hr(),
    dbc.Row([
        dbc.Col(dbc.Tabs(id='tabs-mother', children=[
            dbc.Tab(
                label='Transplants',
                children=dbc.Row([
                    dbc.Col([
                        dbc.Tabs(id='tabs-children', 
                            children=[
                                dbc.Tab(label='About', tab_id='tab-1-1'),
                                dbc.Tab(label='Data', tab_id='tab-1-2'),
                        ]),
                        html.Div(id='tabs-children-content')
                    ], width={"size": 4, "order": 1, "offset": 0}),
                    dbc.Col(
                        [
                            html.H5(id='tx-table-records', children=[], style={'padding':5}),
                            html.Div(
                                id='tx-table',
                                children=[],
                            )          
                        ],
                        width={"size": 8, "order": 2, "offset": 0}
                    ),
            ])), 
            dbc.Tab(label='3D View', 
                    children=html.Div([
                        dcc.Loading(
                            id='3d-view-loading', 
                            type="circle",
                        )
                    ]),
            ),
            dbc.Tab(label='HLA-Epitope 3D View', 
                    children=html.Div([
                        dcc.Loading(
                            id='hla-epitope-3d-view-loading', 
                            type="circle",
                        )
                    ]),
            ),
        ])),
    ]),
], fluid=False,)
# ######################################################################
# App Layout Pieces
# ######################################################################  
####################################################################
show_button =  dbc.Button('Show Table', id='show-table', n_clicks=0)
sortby_dropdown = [
    html.H6('Sort By'),
    dcc.Dropdown(
                id= 'dropdown_sortby',
                options=[
                    {'label': 'Early Failure', 'value': 'early_failure'},
                    {'label': 'Late Surviving', 'value': 'late_surviving'},
                    {'label': 'DESA', 'value': 'desa'}
                ],
                placeholder="",
    ),
]
hla_class = [
            html.H6('HLA Class'),
            dcc.Dropdown(
                id= 'dropdown_class',
                options=[
                    {'label': 'I', 'value': 'I'},
                    {'label': 'II', 'value': 'II'},
                    {'label': 'I & II', 'value': 'III'}
                ],
                placeholder="HLA Class",
            ),
]
hla_molecule = [
    html.H6('HLA Molecule'),
    dbc.Input(id='input_hla', type='text', placeholder="HLA"),
]
donor_type = [
    html.H6('Donor Type'),
    dcc.Dropdown(
                id= 'dropdown_donor_type',
                options=[
                    {'label': 'Living', 'value': 'Living'},
                    {'label': 'Deceased', 'value': 'Deceased'},
                ],
                placeholder="",
    ),
]
ellipro_score = [
    html.H6('ElliPro Score'),
    dcc.Dropdown(
                id= 'dropdown_elli_pro',
                options=[
                    {'label': 'High', 'value': 'High'},
                    {'label': 'Intermediate', 'value': 'Intermediate'},
                    {'label': 'Low', 'value': 'Low'},
                    {'label': 'Very Low', 'value': 'Very Low'},
                ],
                placeholder="",
                multi=True,
    ),
]
filter_card = dbc.Card(
    [
        dbc.CardHeader(html.H5('Filtering',  className="card-title")),
        dbc.CardBody(
            [
                dbc.Row(
                    [
                        dbc.Col(show_button, style={'padding':5}, align="center")
                    ],
                ),
                dbc.Row(
                    [
                        dbc.Col(donor_type,  style={'padding':5}),
                        dbc.Col(sortby_dropdown,  style={'padding':5})
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(ellipro_score,  style={'padding':5})
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(hla_class,  style={'padding':5}),
                        dbc.Col(hla_molecule,  style={'padding':5})
                    ]
                )
            ]
        )
    ]
)
mAb_switch = daq.BooleanSwitch(
                id='mAb-switch',
                on=False,
                label="Monoclonal Abs",
                labelPosition="top",
            )
rAb_switch = daq.BooleanSwitch(
                id='rAb-switch',
                on=False,
                label="Reactive Abs",
                labelPosition="top",
            )
transplant_id = [
        html.H6('Transplant ID'),
        dbc.Input(id='input-tx', type='text', placeholder="Tx ID")
    ]
style_dropdown = [
    html.H6('Style'),
    dcc.Dropdown(
                id= 'dropdown_style',
                options=[
                    {'label': 'Stick', 'value': 'stick'},
                    {'label': 'Cartoon', 'value': 'cartoon'},
                    {'label': 'Sphere', 'value': 'sphere'}
                ],
                placeholder="Style",
    )
]
vis_buttion = dbc.Button('Visualise', id='submit-tx-show', n_clicks=0)
Tx_vis_card = dbc.Card(
    [
        dbc.CardHeader(html.H5('Visualise Transplants:',  className="card-title")),
        dbc.CardBody(
            [
                dbc.Row(
                    [
                        dbc.Col(style_dropdown, style={'padding':5}),
                        dbc.Col(transplant_id, style={'padding':5})
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(mAb_switch, style={'padding':5}),
                        dbc.Col(rAb_switch, style={'padding':5})
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(vis_buttion, style={'padding':5})
                    ]
                )
            ]
        )
    ]
)
ep_vis_buttion = dbc.Button('Visualise', id='submit-hla-ep-show', n_clicks=0)
ep_vis_card = dbc.Card(
    [
        dbc.CardHeader(html.H5('Visualise HLA-Epitopes',  className="card-title")),
        dbc.CardBody(
            [
                dbc.Row(
                    dbc.Textarea(
                        id='input-textarea',
                        bs_size='md', 
                        className="mb-3", 
                        placeholder="Epitopes for visualisation"
                    )
                ), 
                dbc.Row(
                        ep_vis_buttion, style={'padding':5}
                )
            ]
        )
    ]
)

tab_2_layout = html.Div(
    [
        # dbc.CardHeader(html.H3('Filtering & Visualization')),
        dbc.Row(filter_card, className="w-75 mb-4"),
        dbc.Row(Tx_vis_card, className="w-75 mb-4"),
        dbc.Row(ep_vis_card, className="w-75 mb-4")
    ]    
)

# ######################################################################
# Callbacks
# ######################################################################


@app.callback(Output('tabs-children-content', 'children'),
              [Input('tabs-children', 'active_tab')])
def render_content(tab):
    if tab == 'tab-1-1':
        return dbc.Card([
                    dbc.CardHeader(html.H4('What is HLA-Epitope 3D ?')),
                    dbc.CardBody(
                        html.P("""
                        HLA-Epitope 3D is a visualizer that allows you to view the epitopes on HLA molecules
                        with multiple representations: sticks, spheres, and cartoons and color code for 
                        epitopes that are recognised by reactive and human monoclonal antibodies. 
                        You can search and filter the Procare database for all the HLA epitopes to which 
                        antibody can bind.  
                        """)
                    )
                ], style={"width": "4"})
    elif tab == 'tab-1-2':
        return tab_2_layout
    # elif tab == 'tab-1-3':
    #     return None



@app.callback([Output('tx-table', 'children'),
               Output('tx-table-records', 'children')],
              [Input('show-table','n_clicks'),],
              [State('dropdown_sortby', 'value'),
               State('dropdown_class', 'value'),
               State('input_hla', 'value'),
               State('dropdown_donor_type', 'value'),
               State('dropdown_elli_pro', 'value')])
def update_output_data(n_clicks, sort_failure, sort_class, hla, donor_type, ellipro_score):
    if n_clicks == 0:
        return html.P('Click on "Show Table" button to see the table'), None
    else:
        desa_db = load_desa_db()
        ep_db = load_epitope_db()
        desa_db = filtering_logic(desa_db, ep_db, sort_failure, sort_class, hla, donor_type, ellipro_score)
        df_print = dashtable_data_compatibility(desa_db)
        return dbc.Table.from_dataframe(
            df_print, 
            bordered=True,
            dark=True,
            hover=True,
            striped=True,
            size='sm',
            ), f'Records: {len(desa_db)}'

def vis_payload(i, hla, TxID, hlavsdesa, _3d_data):
    """
    This wraps the visualisation payload into a function
    """
    
    return [
            html.H6(
                    html.Span(f'{hla}',
                        id=f"tooltip-target-{TxID}-{i}",
                        style={"textDecoration": "underline", "cursor": "pointer"},
                        className="card-subtitle",
                    )
            ),
            dbc.Tooltip(
                f"DESA #{len(hlavsdesa[TxID][hla]['desa'])}: {hlavsdesa[TxID][hla]['desa']}",
                target=f"tooltip-target-{TxID}-{i}",
                placement='top'
            ),
            div_3dviewer(hla, TxID, _3d_data)
        ]
def get_survivaltime_from_df(desa_df, TxID):
    return desa_df[desa_df['TransplantID'] == TxID]['Survival[Y]'].values[0] 

@app.callback(Output('3d-view-loading', 'children'),
              Input('submit-tx-show','n_clicks'),
              [State('input-tx', 'value'),
               State('dropdown_style', 'value'),
               State('rAb-switch', 'on')])
def show_3d_tx(n_clicks, TxIDs, style, rAb_switch):
    if n_clicks:
        if TxIDs is None:
            return no_update
        desa_df = load_desa_db()
        
        epitope_db = load_epitope_db()
        TxIDs = list(map(int, TxIDs.split(',')))
        _3d_data, hlavsdesa = data_3dviewer(desa_df, epitope_db, TxIDs, style=style, rAb=rAb_switch)
        return [
                dbc.Card(
                    [
                        dbc.CardHeader(
                            html.H4(f'Transplant ID: {TxID}, Survival Time [Y]: {get_survivaltime_from_df(desa_df, TxID): .2f}', 
                                    className="card-title")
                        ),
                        dbc.CardBody(
                            [
                                dbc.Row(
                                    [
                                        dbc.Col(
                                                vis_payload(i, hla, TxID, hlavsdesa, _3d_data), width=6
                                        ) for i, hla in enumerate(hlavsdesa[TxID].keys())
                                    ]
                                )
                            ]
                        ),
                    ],
                color='secondary', style={'padding': 10}) for TxID in TxIDs 
            ]
    else:
        return 'No Transplant ID is selected for visualisation'


@app.callback(Output('hla-epitope-3d-view-loading', 'children'),
              Input('submit-hla-ep-show','n_clicks'),
              [State('input-textarea', 'value'),])
def show_hla_epitope_3d(n_clicks, epitopes):
    if n_clicks:
        if epitopes is None:
            return no_update
        epitope_db =  Epitope_DB()
        hlavsep = epitope_db.min_hlavsep(epitopes)
        return [
                dbc.Card(
                    [
                        # dbc.CardHeader(
                        #     html.H4(f'', 
                        #             className="card-title")
                        # ),
                        dbc.CardBody(
                            [
                                dbc.Row(
                                    [
                                        dbc.Col(
                                                vis_payload(i, hla, TxID, hlavsdesa, _3d_data), width=6
                                        ) for i, hla in enumerate(hlavsdesa[TxID].keys())
                                    ]
                                )
                            ]
                        ),
                    ],
                color='secondary', style={'padding': 10}) for TxID in TxIDs 
            ]
    else:
        return 'No HLA Epitope is given for visualisation'


if __name__ == '__main__':
    if Repository('.').head.shorthand in ['master']:
        debug = False
    else:
        debug = True
    app.run_server(debug=debug, host = '0.0.0.0', port=8080)

# Error 1272
