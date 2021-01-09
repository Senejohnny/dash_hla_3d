# import sys
# print(sys.path)
from app.visualisation import VisualiseHLA

eps = ['105S', '113HN', '114H', '114Q', '116L', '131S', '144QL', '44RME', '62EE', '62QE', '63NI', '65QIA', '66IS', '66IY', '66NH', '70IAQ', '71TD', '74Y', '77D','99S', '9H']

vis = VisualiseHLA()
vis.from_epitopes(set(eps))

print(vis.vis_data_from_epitopes.keys())

# TxIDs = {1402}
# vis = HLA_Epitope_Vis()
# vis.from_transplant(TxIDs)
# print(vis.txvshlavsep)