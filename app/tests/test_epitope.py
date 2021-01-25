""" test scripts for app.epitopes.py """
import pytest
from app.epitope import Epitope
from app.common.utilities import flatten2set
import sys
print(sys.path)


@pytest.mark.parametrize('test_input, expected_output',
[
    ({'3C':'A*12:01', '25I':'A*12:01'}, {'A*12:01': set(['3C', '25I'])}),

])
def test_epvshla2hlavsep(test_input, expected_output):
    assert Epitope().epvshla2hlavsep(test_input) == expected_output



def test_min_hlavsep():
    epitopes = set(['105S', '113HN', '114H', '114Q', '116L', '131S', '144QL',
                    '44RME', '62EE', '62QE', '63NI', '65QIA', '66IS', '66IY',
                    '66NH', '70IAQ', '71TD', '74Y', '77D','99S', '9H'])
    epitope = Epitope()
    _min_hlavsep = epitope.min_hlavsep(epitopes)
    # print(epitope.path)
    print(_min_hlavsep)
    assert flatten2set(_min_hlavsep.values()).intersection(epitopes) == epitopes

test_min_hlavsep()