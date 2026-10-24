# SPDX-License-Identifier: GPL-2.0-or-later
"""Compare the float32 API with a separately installed PySpharm reference.

Run this in an environment containing both PySpharm 1.0.9 and an installed
pyspharm-cpp wheel. The packages use distinct import names, so neither shadows
the other. The regular-grid highest-degree analysis exclusion is intentional;
the reference comparison uses only the safely recoverable degrees.
"""

import numpy as np
import spharm
import spharm_cpp


def compare():
    for m in range(1, 20):
        lat_ref, lon_ref = spharm.getgeodesicpts(m)
        lat_new, lon_new = spharm_cpp.getgeodesicpts(m)
        assert lat_ref.shape == lat_new.shape == (10 * (m - 1) ** 2 + 2,)
        np.testing.assert_allclose(lat_ref, lat_new, rtol=0, atol=4e-5)
        # Both 0 and 360 degrees describe the same seam. The legacy routine's
        # float32 coordinate roundoff can select either endpoint.
        lon_error = (lon_ref.astype(np.float64) - lon_new.astype(np.float64) + 180) % 360 - 180
        assert np.max(np.abs(lon_error)) < 1e-4

    for nlat in (1, 2, 3, 8, 9, 16):
        for ref, new in zip(spharm.gaussian_lats_wts(nlat),
                            spharm_cpp.gaussian_lats_wts(nlat)):
            np.testing.assert_allclose(ref, new, rtol=2e-14, atol=2e-14)

    for truncation in (0, 1, 8, 32):
        m_ref, n_ref = spharm.getspecindx(truncation)
        m_new, n_new = spharm_cpp.getspecindx(truncation)
        np.testing.assert_array_equal(m_ref, m_new)
        np.testing.assert_array_equal(n_ref, n_new)
        for latitude in (-90.0, -89.999, -37.0, 0.0, 37.25, 89.999, 90.0):
            np.testing.assert_allclose(
                spharm.legendre(latitude, truncation),
                spharm_cpp.legendre(latitude, truncation, dtype=np.float32),
                rtol=2e-5,
                atol=2e-5,
            )

    for gridtype in ("regular", "gaussian"):
        old = spharm.Spharmt(18, 9, 2.0, gridtype=gridtype)
        new = spharm_cpp.Spharmt(18, 9, 2.0, gridtype=gridtype,
                                 dtype=np.float32)
        truncation = 6
        rng = np.random.default_rng(41)
        u = rng.normal(size=(9, 18)).astype(np.float32)
        v = rng.normal(size=(9, 18)).astype(np.float32)
        spec_old = old.grdtospec(u, truncation)
        spec_new = new.grdtospec(u, truncation)
        np.testing.assert_allclose(spec_old, spec_new, rtol=5e-5, atol=5e-5)
        np.testing.assert_allclose(old.spectogrd(spec_old),
                                   new.spectogrd(spec_old),
                                   rtol=5e-5, atol=5e-5)
        pnm_old = spharm.legendre(21.0, truncation)
        pnm_new = spharm_cpp.legendre(21.0, truncation, dtype=np.float32)
        np.testing.assert_allclose(
            spharm.specintrp(50.0, spec_old, pnm_old),
            spharm_cpp.specintrp(50.0, spec_old, pnm_new),
            rtol=5e-5,
            atol=5e-5,
        )

        vrt_old, div_old = old.getvrtdivspec(u, v, truncation)
        vrt_new, div_new = new.getvrtdivspec(u, v, truncation)
        np.testing.assert_allclose(vrt_old, vrt_new, rtol=8e-5, atol=8e-5)
        np.testing.assert_allclose(div_old, div_new, rtol=8e-5, atol=8e-5)
        for ref, actual in zip(old.getuv(vrt_old, div_old),
                               new.getuv(vrt_old, div_old)):
            np.testing.assert_allclose(ref, actual, rtol=8e-5, atol=8e-5)
        for ref, actual in zip(old.getgrad(spec_old), new.getgrad(spec_old)):
            np.testing.assert_allclose(ref, actual, rtol=8e-5, atol=8e-5)
        for ref, actual in zip(old.getpsichi(u, v, truncation),
                               new.getpsichi(u, v, truncation)):
            np.testing.assert_allclose(ref, actual, rtol=1e-4, atol=1e-4)

        if gridtype == "gaussian":
            smooth = np.ones(9, dtype=np.float32)
            # squeeze only the legacy trailing singleton axis documented in
            # the progress notes; the new API follows the documented rank-2.
            np.testing.assert_allclose(
                np.squeeze(old.specsmooth(u, smooth)),
                new.specsmooth(u, smooth),
                rtol=6e-5,
                atol=6e-5,
            )

    source_old = spharm.Spharmt(18, 9, gridtype="gaussian")
    target_old = spharm.Spharmt(16, 8, gridtype="regular")
    source_new = spharm_cpp.Spharmt(18, 9, gridtype="gaussian",
                                   dtype=np.float32)
    target_new = spharm_cpp.Spharmt(16, 8, gridtype="regular",
                                   dtype=np.float32)
    field = np.random.default_rng(42).normal(size=(9, 18)).astype(np.float32)
    smooth = np.ones(8, dtype=np.float32)
    np.testing.assert_allclose(
        np.squeeze(spharm.regrid(source_old, target_old, field, 5, smooth)),
        spharm_cpp.regrid(source_new, target_new, field, 5, smooth),
        rtol=8e-5,
        atol=8e-5,
    )

    print("PySpharm 1.0.9 float32 reference comparisons passed.")


if __name__ == "__main__":
    compare()
