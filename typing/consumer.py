# SPDX-FileCopyrightText: 2026 Albert M. W. Yau
# SPDX-License-Identifier: GPL-2.0-or-later
"""Representative consumer calls used to validate the packaged stubs."""

from typing import cast

import numpy as np
from numpy.typing import NDArray

from spharm_cpp import (
    Spharmt,
    cpu_info,
    gaussian_lats_wts,
    getgeodesicpts,
    getspecindx,
    legendre,
    regrid,
    specintrp,
)


def exercise_public_types() -> None:
    grid: NDArray[np.float64] = np.zeros((8, 16), dtype=np.float64)
    transform = Spharmt(16, 8, dtype=np.float64)
    spectrum = transform.grdtospec(grid, ntrunc=3)
    restored = transform.spectogrd(spectrum)
    vorticity, divergence = transform.getvrtdivspec(grid, grid, ntrunc=3)
    zonal, meridional = transform.getuv(vorticity, divergence)
    gradient_zonal, gradient_meridional = transform.getgrad(spectrum)
    smooth = np.ones(8, dtype=np.float64)
    smoothed = transform.specsmooth(grid, smooth)
    copied = regrid(transform, transform, grid, ntrunc=3, smooth=smooth)

    latitudes, weights = gaussian_lats_wts(8)
    zonal_indices, degrees = getspecindx(3)
    one_latitude = legendre(20.0, 3, dtype=np.float64)
    latitude_batch = legendre(np.array([-30.0, 0.0, 30.0], dtype=np.float64), 3)
    latitude_batch32 = cast(
        NDArray[np.float32],
        legendre(np.array([-30.0, 30.0], dtype=np.float32), 3, dtype=np.float32),
    )
    point_latitudes, point_longitudes = getgeodesicpts(2)
    interpolated = specintrp(0.0, spectrum, one_latitude)
    interpolated_many = specintrp(
        np.array([0.0, 90.0], dtype=np.float64),
        spectrum,
        np.stack([one_latitude, one_latitude]),
    )
    spectrum32: NDArray[np.complex64] = np.zeros(10, dtype=np.complex64)
    interpolated_many32 = specintrp(
        np.array([0.0, 90.0], dtype=np.float64), spectrum32, latitude_batch32
    )
    diagnostics = cpu_info()

    _ = (
        restored,
        zonal,
        meridional,
        gradient_zonal,
        gradient_meridional,
        smoothed,
        copied,
        latitudes,
        weights,
        zonal_indices,
        degrees,
        latitude_batch,
        point_latitudes,
        point_longitudes,
        interpolated,
        interpolated_many,
        interpolated_many32,
        diagnostics,
    )
