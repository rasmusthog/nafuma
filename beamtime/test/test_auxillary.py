import beamtime.auxillary as aux

def test_swap_values():


    dict = {'test1': 1, 'test2': 2}
    key1 = 'test1'
    key2 = 'test2'

    oldval1 = dict[key1]
    oldval2 = dict[key2]

    new_dict = aux.swap_values(dict=dict, key1=key1, key2=key2) 

    assert (dict[key1] == oldval2) and (dict[key2] == oldval1)


def test_ceil() -> None:

    assert aux.ceil(1.05, 0.5) == 1.5
    assert aux.ceil(1.05, 1) == 2.0
    assert aux.ceil(1.1, 0.2) == 1.2


def test_floor() -> None:

    assert aux.floor(2.02, 1) == 2.0
    assert aux.floor(2.02, 0.01) == 2.02
    assert aux.floor(2.013, 0.01) == 2.01