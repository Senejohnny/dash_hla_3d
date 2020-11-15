import json
import pandas as pd 
from app.analysis_utils import load_data, hla_to_filename, get_hla_locus, find_molecule_path, flatten2list, hlavsdesa_donor, _3dview_data_preparation
import pytest 
import warnings
warnings.filterwarnings("ignore")


ep_db, desa_db = load_data()
# print(desa_db['TransplantID'].tolist())
# print(flatten2list(desa_db.Donor_HLA.values))
# ['TransplantID', 'DESA_Epitope', '#DESA', 'EpvsHLA_Pos', 'EpvsHLA_Donor',
#        'Failure', 'Survival[Y]', 'DESA->Donor_HLA', 'Donor_HLA',
#        'Donor_HLA_Class']
# print(desa_db[desa_db['TransplantID'] == 4080].EpvsHLA_Donor.values[0].keys())
# print(desa_db[desa_db.TransplantID == 4080].EpvsHLA_Donor.values)
# print(ep_db.info())

# def test_polymorphic_residues():


def test_get_hla_locus():
    data = {'A*27:01': 'A', 'DRB3*01:01': 'DRB',
            'DQA*01:01': 'DQA', 'DQB*01:01': 'DQB'}
    for hla in data.keys():
        locus = get_hla_locus(data[hla])
    assert locus == data[hla] 


def test_hla_to_filename():
    test_set = {
        'A*12:12': {'Locus': 'A', 'name':'A_12_12_V1.pdb'},
        'DRB1*17:19': {'Locus': 'DRB', 'name':'DRB1_17_19_V1.pdb'},
        'DQA*01:01': {'Locus': 'DQA', 'name':'DQA_01_01_V1.pdb'},
    }
    for hla in test_set.keys():
        locus, filename = hla_to_filename(hla)
        assert locus == test_set[hla].get('Locus')
        assert filename == test_set[hla].get('name')

def test_find_molecule_path():
    all_hla = set(flatten2list(desa_db.Donor_HLA.values))
    print(f"All the donor HLA's are: {len(all_hla)}")
    for hla in all_hla:
        locus, filename = hla_to_filename(hla)
        pdb_exist, pdb_path = find_molecule_path(locus, filename)
        assert pdb_exist == True


def test_hlavsdesa_donor():
    TxIDs = desa_db['TransplantID'].values.tolist()
    hlavsdesa = hlavsdesa_donor(ep_db, desa_db, TxIDs, rAb=True, mAb=True)
    for TxID in TxIDs:
        desavshla = desa_db[desa_db['TransplantID'] == TxID].EpvsHLA_Donor.values[0]
        assert set(hlavsdesa[TxID].keys()) == set(desavshla.values())
        assert set(desavshla.keys()) == set(flatten2list([hlavsdesa[TxID][key]['desa'] 
                                                            for key in hlavsdesa[TxID].keys()]))

@pytest.mark.skip(reason='No Reason')
def test_3dview_data_preparation():
    TxIDs = [4080, 5552, 833] 
    # TxIDs = desa_db['TransplantID'].values.tolist()
    # print(f'TxID:{TxIDs}')
    hlavsdesa = hlavsdesa_donor(ep_db, desa_db, TxIDs, rAb=True, mAb=True)
    _3d_models_data = _3dview_data_preparation(hlavsdesa, 'sphere')

test_3dview_data_preparation()
# test_find_molecule_path()
# pdb_exist, pdb_path = find_molecule_path(locus, filename)
# test_hla_locus_filename()
