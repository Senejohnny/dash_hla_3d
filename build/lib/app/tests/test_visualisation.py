""" test scripts for visualisation """
from app.visualisation import VisualiseHLA


def test_from_transplants():
    TxIDs = {1402}
    hlas = {'A*11:01', 'C*03:03'}
    vis = VisualiseHLA()
    vis.from_transplant(TxIDs)
    for TxID in TxIDs:
        assert TxID in vis.txvshlavsep.keys()
        for hla in hlas:
            assert hla in vis.txvshlavsep.get(TxID).keys()

# test_from_transplants()

def test_from_epitopes():
    epitopes = set(['105S', '113HN', '114H', '114Q', '116L', '131S', '144QL',
                    '44RME', '62EE', '62QE', '63NI', '65QIA', '66IS', '66IY',
                    '66NH', '70IAQ', '71TD', '74Y', '77D','99S', '9H'])
    vis = VisualiseHLA()
    vis.from_epitopes(epitopes)
    print(vis.hlavsep)

# print(vis.vis_data_from_epitopes.keys()),,
test_from_epitopes()
