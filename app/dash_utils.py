""" some dash utility functions """
import io 
import base64
import pandas as pd
from dash import no_update
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from app.desa import DESA

# def filtering_logic(df,
#                     ep_db,
#                     sort_failure, 
#                     sort_class, 
#                     hla, 
#                     donor_type,
#                     ellipro_score):
#     """ This function wraps all the logic pertain to
#     filtering the table """
#     if sort_failure:
#         if sort_failure == 'desa':  
#             df = df.sort_values(by='#DESA', ascending=False)

#         if sort_failure == 'early_failure': 
#             ind_T_early = df['Survival[Y]'].apply(lambda x: x < 1/4)
#             ind_E_early = df['Failure'].apply(lambda x: x == 1)
#             df = df[ind_T_early & ind_E_early]
            
#         if sort_failure == 'late_surviving': 
#             ind_T_late = df['Survival[Y]'].apply(lambda x: x > 10) # Still have Kidney after 10 years
#             ind_E_late = df['Failure'].apply(lambda x: x == 0 | x == 2)
#             df = df[ind_T_late & ind_E_late]
#     if sort_class:
#         ind = df['Donor_HLA_Class'].apply(lambda x: x == {'I':'I', 'II': 'II', 'III': 'I,II'}.get(sort_class))
#         df = df[ind]

#     if hla:
#         ind = df['Donor_HLA'].apply(lambda x: hla in x)
#         df = df[ind]
#     if donor_type:
#         ind = df['Donor_Type'].apply(lambda x: x == donor_type)
#         df = df[ind]
#     if ellipro_score:
#         try:
#             df.drop('EpvsHLA_Donor_filt', axis=1)
#         except:
#             pass
#         ind = ep_db['ElliPro Score'].apply(lambda x: x in ellipro_score)
#         ep_ellipro = ep_db[ind].Epitope.values.tolist()
#         df['EpvsHLA_Donor_filt'] = df.EpvsHLA_Donor.apply(lambda x: 
#                         {key:value for key, value in x.items() 
#                                     if key in ep_ellipro})
#         df.to_pickle('./data/desa_3d_view.pickle')
#     return df

def filtering_logic(desa: DESA, sort_by:str, hla_class:str, hla:str, donor_type:str):
    """ This function wraps all the logic pertain to filtering the table """
    if donor_type:
        desa.donor_type(donor_type)
    if hla_class:
        desa.hla_class(hla_class)
    if hla:
        desa.tx_with_hla(hla)
    if sort_by:
        if sort_by == 'early_failure':
            desa.early_failed(1/4)
        if sort_by == 'late_surviving':
            desa.late_failed(10)
        if sort_by == 'desa':
            return desa.df.sort_values(by='#DESA', ascending=False)
    return desa.df


def dashtable_data_compatibility(df):
    """ This function helps to make the data compatible to dash table
    by turning all unstructured values into a string """

    df['Survival[Y]'] = df['Survival[Y]'].apply(lambda x: round(x,3))
    df['DESA->Donor_HLA'] = df['DESA->Donor_HLA'].apply(lambda x: str(dict(x)))
    df['Donor_HLA'] = df['Donor_HLA'].apply(lambda x: str(x)) 
    df['Donor_HLA_Class'] = df['Donor_HLA_Class'].apply(lambda x: str(x))
    return df[['TransplantID', '#DESA', 'Failure', 'Survival[Y]', 'Donor_HLA', 'Donor_HLA_Class']]

def Header(name):
    title = html.H1(name, style={"margin-top": 5})
    logo = html.Img(
                src="assets/logo.svg",
                style={
                        "float": "right", 
                        "height": 60,
                        'float': 'right',
                        'position': 'relative',
                        'padding-top': 4, 
                        'padding-right': 4,
                }
            )
    link = html.A(logo, href="https://plotly.com/dash/")

    return dbc.Row([dbc.Col(title, md=8), dbc.Col(link, md=4)])
def vis_payload(hlavsdesa, _3d_data, hla, i):
    return [
                                html.H6(
                                        html.Span(f'{hla}',
                                            id=f"tooltip-target-{i}",
                                            style={"textDecoration": "underline", "cursor": "pointer"},
                                        )
                                ),
                                dbc.Tooltip(
                                    f"DESA #{len(hlavsdesa[hla]['desa'])}: {hlavsdesa[hla]['desa']}",
                                    target=f"tooltip-target-{i}",
                                    placement='top'
                                ),
                            div_3dviewer(hla, _3d_data)]
#####################################################################################
# File Upload
#####################################################################################
def upload_file(Data_Type:str, id:str = ''):
    return dcc.Upload(id=id,
                children=html.Div(
                [
                    # f'Drag and Drop {Data_Type} Data ',
                    html.A(f'Upload {Data_Type} File')
                ]),
                style={
                    'width': '60%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '20px'
                },
                # Allow multiple files to be uploaded
                multiple=False,
           )


#####################################################################################
# Parse Uploaded file
#####################################################################################
def parse_contents(contents, filename, Data_type:str)-> pd.DataFrame: # the df output will not be stored so no json format.  

    if '.' not in filename:
        return html.Div(f"""Error: filename {filename} does not contain extension. Make sure filename contains
        supported file extension e.g. .csv, .xlsx.""")

    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string) # decode byte-like object or ASCII into bytes 
    try:
        if filename.rsplit('.', 1)[1].lower() == 'csv':
            # Assume that the user uploaded a CSV file
            string = decoded.decode('utf-8') # decode bytes into string 
            file = io.StringIO(string)    # StringIO is an in-memory stream for text I/O
            delimiter = detectdelimiter(file) # find the delimiter of the csv file
            df = pd.read_csv(file, sep=delimiter)
            file.close()                       # close StringIO object
        elif filename.rsplit('.', 1)[1].lower() in ['xls', 'xlsx']:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))     #BytesIO is an in-memory stream for text I/O
        elif filename.rsplit('.', 1)[1].lower() in ['pkl', 'pickle']:
            # Assume that the user uploaded an pickle file
            df = pd.read_pickle(io.BytesIO(decoded))     #BytesIO is an in-memory stream for text I/O
    except Exception as e:
        print(e)
        return no_update, html.Div(["""
                                There was an error processing this file. 
                                File extension might not supported. 
                                Supported file extensions .csv, .xls, .xlsx, .pickle
                                """])
    return df, html.Div(html.P(f'{Data_type} data is successfully uploaded as {filename}!'))


def detectdelimiter(file): 
    """
    this will detect the famous delimiters [',', ';', ':', '|', '\t'] 
    The input is _io.StringIO object. 
    """

    first_line = file.getvalue().split('\n')[0] # get the first line of the StringIO object
    return detect(first_line)   # apply the detect method