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
from app.analysis_utils import load_desa_db, data_3dviewer, div_3dviewer, load_epitope_db

from app.dash_utils import parse_contents, upload_file
from pygit2 import Repository


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, title='HLA 3D Viewer', external_stylesheets=external_stylesheets)
app.config['suppress_callback_exceptions'] = True

# app.css.append_css({
#     "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
# })

# epitope_db, desa_db = load_data()

# desa_db = desa_db[desa_db.Status=='DESA']
# Tx = 4606



def data_to_str(df):
    df['Survival[Y]'] = df['Survival[Y]'].apply(lambda x: round(x,3))
    df['DESA->Donor_HLA'] = df['DESA->Donor_HLA'].apply(lambda x: str(dict(x)))
    df['Donor_HLA'] = df['Donor_HLA'].apply(lambda x: str(x)) 
    df['Donor_HLA_Class'] = df['Donor_HLA_Class'].apply(lambda x: str(x))
    return df[['TransplantID', 'Status', '#DESA', 'Failure', 'Survival[Y]', 'Donor_HLA', 'Donor_HLA_Class']]

def filtering_logic(df, sort_failure, sort_class, hla):
    if sort_failure:
        if sort_failure == 'desa':  
            df = df.sort_values(by='#DESA', ascending=False)

        if sort_failure == 'early_failure': 
            ind_T_early = df['Survival[Y]'].apply(lambda x: x < 1/4)
            ind_E_early = df['Failure'].apply(lambda x: x == 1)
            df = df[ind_T_early & ind_E_early]
        print(sort_failure)
        if sort_failure == 'late_surviving': 
            ind_T_late = df['Survival[Y]'].apply(lambda x: x > 10)
            ind_E_late = df['Failure'].apply(lambda x: x == 0 | x == 2)
            df = df[ind_T_late & ind_E_late]
    if sort_class:
        ind = df['Donor_HLA_Class'].apply(lambda x: x == {'I':'I', 'II': 'II', 'III': 'I,II'}.get(sort_class))
        df = df[ind]

    if hla:
        ind = df['Donor_HLA'].apply(lambda x: hla in x)
        df = df[ind]

    return df


def logo_img():
    return html.Img(
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
                    ]),
                    html.Div(id='table',
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
                ], className="four columns"),

tab_3_layout = html.Div([
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
                        html.Button('Show', id='submit-tx-show', n_clicks=0),
                    ])
                ], 
                className="four columns")


                
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
        ], className="four columns")
    elif tab == 'tab-1-2':
        return tab_2_layout
    elif tab == 'tab-1-3':
        return tab_3_layout



@app.callback(Output('table', 'children'),
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
        df_print = data_to_str(desa_db)
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

# @app.callback([Output('output-desa-upload', 'children'),
#                Output('hidden-df', 'children'),
#                Output('table', 'children'),],
#               [Input('upload-desa-data', 'contents')],
#             #    Input('submit-filter-data','n_clicks'),],
#               [State('upload-desa-data', 'filename'),
#             #    State('dropdown_sortby', 'value'),
#             #    State('dropdown_class', 'value'),
#             #    State('input_hla', 'value')           
#                ])
# def update_output_data(contents, filename):
#     if contents is None:
#         return html.Div(['No desa DB is uploaded']), no_update, no_update 
#     else:
#         df, message = parse_contents(contents, filename, 'Transplants with DESA')
#         # Any filtering should be before as data_eng stringify values
#         df_eng = data_eng(df)
#         return message, df.to_json(orient='split', date_format='iso'), \
#             dash_table.DataTable(
#                                 columns=[{"name": i, "id": i} for i in df_eng.columns],   
#                                 data=df_eng.to_dict('records'),
#                                 page_size= 20,
#                                 editable=True,
#                                 style_table={'height': '400px', 
#                                              'overflowY': 'auto'},
#                                 style_cell_conditional=[
#                                                         {
#                                                             'if': {'column_id': c},
#                                                             'textAlign': 'left'
#                                                         } for c in ['Date', 'Region']
#                                                     ],
#                                 style_data_conditional=[
#                                                         {
#                                                             'if': {'row_index': 'odd'},
#                                                             'backgroundColor': 'rgb(248, 248, 248)'
#                                                         }
#                                                     ],
#                                 style_header={
#                                                 'backgroundColor': 'rgb(230, 230, 230)',
#                                                 'fontWeight': 'bold'
#                                             }
#                         )




# @app.callback([Output('hidden-data', 'children'),],
#               [Input('submit-filter-data','n_clicks'),],
#               [State('dropdown_sortby', 'value'),
#                State('dropdown_class', 'value'),
#                State('input_hla', 'value')           
#                ])
# def update_output_data(n_clicks, sort_failure, sort_class, hla):
#     if n_clicks:
#         desa_db = load_desa_db()
#         desa_db = filtering_logic(desa_db, sort_failure, sort_class, hla)
#         return desa_db.to_json(date_format='iso', orient='split')
#     else:
#         return no_update


    
@app.callback(Output('3d-view-loading', 'children'),
              Input('submit-tx-show','n_clicks'),
              [State('input-tx', 'value'),
               State('dropdown_style', 'value')])
def show_3d_tx(n_clicks, TxID, style):
    if n_clicks:
        desa_df = load_desa_db()
        epitope_db = load_epitope_db()
        hlas, _3d_data, _hlavsdesa, hlavsdesa = data_3dviewer(desa_df, epitope_db, int(TxID), style=style)
        return [html.Div(
                    children=[
                        html.H5(f'HLA: {hla}'),
                        div_3dviewer(hla, _3d_data)], 
                        className="four columns") for hla in hlas]
    else:
        return 'No Transplant ID is selected for visualisation'



if __name__ == '__main__':
    if Repository('.').head.shorthand == 'master':
        debug = False
    else:
        debug = True
    app.run_server(debug=debug, host = '0.0.0.0', port=8080)
