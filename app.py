""" The core functionality of the app is here
    written by Danial Senejohnny """
# from app.epitope import Epitope
import warnings
import dash
from dash import no_update
from dash.dependencies import Input, Output, State
import dash_daq as daq
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from pygit2 import Repository

from app.visualisation import VisualiseHLA, vis_cards
from app.desa import DESA
from app.dash_utils import filtering_logic, Header, dashtable_data_compatibility

warnings.filterwarnings("ignore")

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

# cache = Cache(app.server, config={
#     'CACHE_TYPE': 'redis',
#     # Note that filesystem cache doesn't work on systems with ephemeral
#     # filesystems like Heroku.
#     'CACHE_TYPE': 'filesystem',
#     'CACHE_DIR': 'cache-directory',

#     # should be equal to maximum number of users on the app at a single time
#     # higher numbers will store more data in the filesystem / redis cache
#     'CACHE_THRESHOLD': 200
# })




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
            dbc.Tab(label='Visualise Transplants',
                    children=html.Div([
                        dcc.Loading(
                            id='3d-view-loading',
                            type="circle",
                        )
                    ]),
            ),
            dbc.Tab(label='Visualise Epitopes',
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
                    {'label': 'I & II', 'value': 'I,II'}
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
# exclude_nev = daq.BooleanSwitch(        # pylint: disable=not-callable
#                 id='exclude-nev-switch',
#                 on=False,
#                 label="Exclude Never Functioning Grafts",
#                 labelPosition="top",
#             )

filter_card = dbc.Card(
    [
        dbc.CardHeader(html.H5('Data Base',  className="card-title")),
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
                        dbc.Col(hla_class,  style={'padding':5}),
                        dbc.Col(hla_molecule,  style={'padding':5})
                    ]
                ),
                # dbc.Row(
                #     [
                #         dbc.Col(exclude_nev,  style={'padding':5}),

                #     ]
                # )
            ]
        )
    ]
)
mAb_switch = daq.BooleanSwitch(        # pylint: disable=not-callable
                id='mAb-switch',
                on=False,
                label="Monoclonal Abs",
                labelPosition="top",
            )
# rAb_switch = daq.BooleanSwitch(        # pylint: disable=not-callable
#                 id='rAb-switch',
#                 on=False,
#                 label="Reactive Abs",
#                 labelPosition="top",
#             )
transplant_id = [
        html.H6('By Transplant'),
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
vis_buttion_tx = dbc.Button('Visualise Transplant', id='submit-tx-show', n_clicks=0)
vis_buttion_epitope = dbc.Button('Visualise Epitopes', id='submit-ep-show', n_clicks=0)
text_area = [
        html.H6('By Epitopes'),
        dbc.Textarea(
                    id='input-textarea',
                    bs_size='md',
                    className="mb-3",
                    placeholder="Epitopes for visualisation"
                )
]
Tx_vis_card = dbc.Card(
    [
        dbc.CardHeader(html.H5('3D Visualisation:',  className="card-title")),
        dbc.CardBody(
            [
                dbc.Row(
                    [
                        dbc.Col(ellipro_score,  style={'padding':5})
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(style_dropdown, style={'padding':5}),
                        dbc.Col(transplant_id, style={'padding':5})
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(text_area, style={'padding':5})
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(mAb_switch, style={'padding':5}),
                    ]
                ),
                dbc.Row(
                    [
                        dbc.Col(vis_buttion_tx, style={'padding':5}),
                        dbc.Col(vis_buttion_epitope, style={'padding':5})
                    ]
                )
            ]
        )
    ]
)
ep_vis_buttion = dbc.Button('Visualise', id='submit-hla-ep-show', n_clicks=0)


tab_2_layout = html.Div(
    [
        # dbc.CardHeader(html.H3('Filtering & Visualization')),
        dbc.Row(filter_card, className="w-75 mb-4"),
        dbc.Row(Tx_vis_card, className="w-75 mb-4"),
        # dbc.Row(ep_vis_card, className="w-75 mb-4")
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

###############################################################
# Callback Table
###############################################################
@app.callback([Output('tx-table', 'children'),
               Output('tx-table-records', 'children')],
              [Input('show-table','n_clicks'),],
              [State('dropdown_sortby', 'value'),
               State('dropdown_class', 'value'),
               State('input_hla', 'value'),
               State('dropdown_donor_type', 'value'),])
def table(n_clicks, sort_by, hla_class, hla, donor_type):
    if n_clicks == 0:
        return html.P('Click on "Show Table" button to see the table'), None
    desa = DESA()
    desa_df = filtering_logic(desa, sort_by, hla_class, hla, donor_type, None)
    desa_df = dashtable_data_compatibility(desa_df)
    return dbc.Table.from_dataframe(
        desa_df,
        bordered=True,
        dark=True,
        hover=True,
        striped=True,
        size='sm',
        ), f'Records: {len(desa_df)}'

###############################################################
# Callback Visualisation from Transplants
###############################################################
@app.callback(Output('3d-view-loading', 'children'),
              Input('submit-tx-show','n_clicks'),
              [State('input-tx', 'value'),
               State('dropdown_style', 'value'),
               State('mAb-switch', 'on'),
               State('dropdown_elli_pro', 'value')])
def visualise_from_transplants(n_clicks, TxIDs, style, mAb_switch, elliproscore):
    if n_clicks:
        if TxIDs == None :
            return no_update
        style = 'sphere' if style == None else style
        TxIDs = set(map(int, TxIDs.split(',')))
        vis = VisualiseHLA(ignore_hla={'B*13:01'})
        vis_object = vis.from_transplant(TxIDs, style, mAb=mAb_switch, elliproscore=elliproscore)
        return vis_cards(vis_object)
    else:
        return 'No Transplant ID is selected for visualisation'

###############################################################
# Callback Visualisation from Epitopes
###############################################################
@app.callback(Output('hla-epitope-3d-view-loading', 'children'),
              Input('submit-ep-show','n_clicks'),
              [State('input-textarea', 'value'),
               State('dropdown_style', 'value'),
               State('mAb-switch', 'on'),
               State('dropdown_elli_pro', 'value')])
def visualise_from_epitopes(n_clicks, epitopes, style, mAb_switch, elliproscore):
    if n_clicks:
        if epitopes is None:
            return no_update
        style = 'sphere' if style == None else style
        epitopes = epitopes.replace("'", "").replace("\n", "")
        epitopes = set(map(str.strip, epitopes.split(',')))
        vis = VisualiseHLA(ignore_hla={'B*13:01', 'DQB1*06:03'})
        vis_object = vis.from_epitopes(epitopes, style, mAb=mAb_switch, elliproscore=elliproscore)
        return vis_cards(vis_object)
        # return no_update
    else:
        return 'No HLA Epitope is given for visualisation'

###############################################################
if __name__ == '__main__':
    if Repository('.').head.shorthand in ['master']:
        debug = False
    else:
        debug = True
    app.run_server(debug=debug, host = '0.0.0.0', port=8080)

