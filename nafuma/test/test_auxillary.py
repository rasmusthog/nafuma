import nafuma.auxillary as aux
import os

def test_swap_values():


    test_dict = {'test1': 1, 'test2': 2}
    key1 = 'test1'
    key2 = 'test2'

    oldval1 = test_dict[key1]
    oldval2 = test_dict[key2]

    new_dict = aux.swap_values(options=test_dict, key1=key1, key2=key2) 

    assert (test_dict[key1] == oldval2) and (test_dict[key2] == oldval1)


def test_ceil() -> None:

    assert aux.ceil(1.05, 0.5) == 1.5
    assert aux.ceil(1.05, 1) == 2.0
    assert aux.ceil(1.1, 0.2) == 1.2


def test_floor() -> None:

    assert aux.floor(2.02, 1) == 2.0
    assert aux.floor(2.02, 0.01) == 2.02
    assert aux.floor(2.013, 0.01) == 2.01



def test_options() -> None:


    options = {}

    default_options = {
        'test1': 1,
        'test2': 2, 
        'test3': 3,
        'test4': 4,
        'test5': 5,
    }


    options = aux.update_options(options=options, default_options=default_options)

    assert options['test1'] == default_options['test1']


def test_save_options() -> None:
    
    options = {'test1': 1, 'test2': 2}
    path = 'tmp.dat'

    aux.save_options(options, path)

    assert os.path.isfile(path)

    os.remove(path)


def test_load_options() -> None:
    
    options = {'test1': 1, 'test2': 2}
    path = 'tmp.dat'

    aux.save_options(options, path)

    loaded_options = aux.load_options(path)

    assert (loaded_options['test1'] == 1) and (loaded_options['test2'] == 2)

    os.remove(path)
