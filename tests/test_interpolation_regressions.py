import numpy as np
import pytest
from microturbogeo.geom import reinterpolate_curve


def test_shifted_parameter_and_duplicate_points():
    x, y, s = reinterpolate_curve([0, 1, 1, 2], [0, 0, 0, 0], [5, 6, 6, 7], N_new=4)
    assert len(s) == 7
    np.testing.assert_allclose(x, s - 5)
    np.testing.assert_allclose(y, 0)
    np.testing.assert_allclose(s[[0, -1]], [5, 7])
    assert np.all(np.diff(s) > 0)


@pytest.mark.parametrize('s', [[0, 0, 1], [0, 2, 1], [0, np.nan, 1]])
def test_invalid_parameters(s):
    with pytest.raises(ValueError):
        reinterpolate_curve([0, 1, 2], [0, 0, 0], s)


@pytest.mark.parametrize('kwargs', [{'base': 1}, {'N_new': 1}, {'N_new': 2.5}, {'N_new': True}])
def test_invalid_options(kwargs):
    with pytest.raises(ValueError):
        reinterpolate_curve([0, 1], [0, 0], [0, 1], **kwargs)
