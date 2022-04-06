import beamtime.plotting as btp
from cycler import cycler
import itertools


def test_generate_colours() -> None:

    assert type(btp.generate_colours('black', kind='single')) == itertools.cycle