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

from app.dash_utils import filtering_logic, Header, dashtable_data_compatibility
from pygit2 import Repository


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, title='HLA3D Epitopes', external_stylesheets=[dbc.themes.CERULEAN])
app.config['suppress_callback_exceptions'] = True

# app.css.append_css({
#     "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
# })





# ######################################################################
# App Layout
# ######################################################################
app.layout = dbc.Container([
    Header('HLA3D Epitopes'),
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
                    dbc.Col(id='tx-table',
                            children=[],
                            width={"size": 8, "order": 2, "offset": 0}
                    ),
            ])), 
            dbc.Tab(label='HLA 3D View', 
                    children=html.Div([
                        # html.H5('HLA Molecule 3D View'),
                        # html.Div(f'''
                        #     This is a 3D view of HLA molecules in Transplant ID: , with DESA color depicted in yellow  
                        # '''),
                        dcc.Loading(
                            id='3d-view-loading', 
                            type="circle",
                        )
                    ]),
            ),
            # dbc.Tab(
            #     label='3D Replication',
            #     className='control-upload',
            #     children=[
            #         html.Div([
            #             html.H5("Which HLA should be repeated?"),
            #             dcc.Input(
            #                 id="rep-hla", 
            #                 placeholder='Enter a HLA',
            #                 type='text',
            #                 value=''
            #             ), 
            #             html.H5("Which DESA Should be marked?"),
            #             dcc.Input(
            #                 id="rep-desa", 
            #                 placeholder="Enter DESA's",
            #                 type='text',
            #                 value=''
            #             )
            #         ],),
            #         html.Div(
            #             html.H4("The 3D HLA structure would be represented below")
            #         )
            #     ]
            # )
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
filter_card = dbc.Card(
    [
        dbc.CardHeader(html.H5('Filtering',  className="card-title")),
        dbc.CardBody(
            [
                dbc.Row(
                    [
                        dbc.Col(show_button),
                        dbc.Col(sortby_dropdown)
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(hla_class),
                        dbc.Col(hla_molecule)
                    ]
                )
            ]
        )
    ]
)
mAb_switch = daq.BooleanSwitch(
                id='mAb-switch',
                on=False,
                label="Monoclonal Antibody",
                labelPosition="top",
            )
rAb_switch = daq.BooleanSwitch(
                id='rAb-switch',
                on=False,
                label="Reactive Antibody",
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
show_buttion = dbc.Button('Show', id='submit-tx-show', n_clicks=0)
vis_card = dbc.Card(
    [
        dbc.CardHeader(html.H5('Visualisation',  className="card-title")),
        dbc.CardBody(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            dbc.Row(
                                [
                                    mAb_switch,
                                    rAb_switch,
                                ],  align="start",
                            ), 
                        ),
                        dbc.Col(transplant_id, )
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(style_dropdown),
                        dbc.Col(show_buttion)
                    ]
                )
            ]
        )
    ]
)

tab_2_layout = html.Div(
    [
        # dbc.CardHeader(html.H3('Filtering & Visualization')),
        dbc.Row(filter_card, className="mb-4"),
        dbc.Row(vis_card, className="mb-4")
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
                    dbc.CardHeader(html.H3('What is HLA3D ?')),
                    dbc.CardBody(
                        html.P("""
                        HLA3D is a visualizer that allows you to view the epitopes on HLA molecules
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
        return dbc.Table.from_dataframe(
            df_print, 
            bordered=True,
            dark=True,
            hover=True,
            striped=True,
            size='sm',
            )

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



if __name__ == '__main__':
    if Repository('.').head.shorthand in ['master']:
        debug = False
    else:
        debug = True
    app.run_server(debug=debug, host = '0.0.0.0', port=8080)

# Error 1272
