""" The utility functions used by the app """
import os
import logging
import pandas as pd


###################################################################
# Auxiliary functions for analysis
###################################################################
def load_data():
    """ This function loads DESA and EpitopeDB data frames"""

    return load_epitope_db(), load_desa_db()

def load_desa_db():
    """ This function loads DESA data frame"""

    desa_path = os.path.expanduser('./data/desa_3d_view.pickle')
    desa_db = pd.read_pickle(desa_path)
    # desa_db = desa_db.drop('Status', axis=1)
    return desa_db

def load_epitope_db():
    """ This function loads EpitopeDB data frame"""

    ep_path = os.path.expanduser('./data/20201123_EpitopevsHLA.pickle')
    epitope_db = pd.read_pickle(ep_path)
    return epitope_db

def flatten2list(nested_object) -> list:
    """ This function flattens objects in a nested structure and returns list"""
    gather = []
    for item in nested_object:
        if isinstance(item, (list, set)):
            gather.extend(flatten2list(item))
        else:
            gather.append(item)
    return gather

def flatten2set(nested_object) -> set:
    """ This function flattens objects in a nested structure and returns a set"""

    return set(flatten2list(nested_object))
