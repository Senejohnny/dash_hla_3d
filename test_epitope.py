""" test scripts for app.epitopes.py """
from app.epitope import Epitope 
from app.utilities import flatten2set


def test_min_hlavsep():
    epitopes = set(['105S', '113HN', '114H', '114Q', '116L', '131S', '144QL',
                    '44RME', '62EE', '62QE', '63NI', '65QIA', '66IS', '66IY',
                    '66NH', '70IAQ', '71TD', '74Y', '77D','99S', '9H'])
    epitope = Epitope()
    _min_hlavsep = epitope.min_hlavsep(epitopes)
    print(_min_hlavsep)
    # for hla in _min_hlavsep.keys():
    #     assert hla in set(['B*55:01', 'A*31:01', 'B*37:01', 'A*23:02', 'B*13:02'])
    assert flatten2set(_min_hlavsep.values()).intersection(epitopes) == epitopes

test_min_hlavsep()