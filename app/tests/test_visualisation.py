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
    epitopes_classI = set(['105S', '113HN', '114H', '114Q', '116L', '131S', '144QL',
                    '44RME', '62EE', '62QE', '63NI', '65QIA', '66IS', '66IY',
                    '66NH', '70IAQ', '71TD', '74Y', '77D','99S', '9H'])
    epitopes_classII = set(['140TV', '149H', '160AD', '180VMP', '180VTP', '181T', '25Q', '25R', '28D',
                            '30G[DR]', '31FH', '32H', '37N', '37S', '46VY', '47F', '4R', '57A', '57S',
                            '66IT', '6C', '74A', '75IL', '96EV', '96HK', '98ES', '98KS'])

    ep = set(['74EL', '73GQ', '104AK', '52LL', '56PA', '55PPA', '125G', '125SQ', '75S', '26L[DR]', 
              '67VT', '32YN', '130A', '74R', '26L[DQ]', '66EV', '116I', '74E', '38L', '40GR', '37YV[DQ]', 
              '70RE', '32Y', '86A', '87Y', '67I', '70R', '37L', '55RPD', '13GM', '87F'])

    vis = VisualiseHLA(ignore_hla={'B*13:01', 'DRB1*13:05'})
    vis.from_epitopes(ep)
    print({key:len(vis.hlavsep[key]) for key in vis.hlavsep.keys()})
    # vis.from_epitopes(epitopes_classII)
    # print({key:len(vis.hlavsep[key]) for key in vis.hlavsep.keys()})

# print(vis.vis_data_from_epitopes.keys()),,
test_from_epitopes()
