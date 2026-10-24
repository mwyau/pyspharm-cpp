import numpy
import numpy.typing

# SPDX-License-Identifier: GPL-2.0-or-later

class Spharmt:
    """Spherical harmonic transforms on regular or Gaussian latitude grids."""

    def __init__(
        self,
        nlon: int,
        nlat: int,
        rsphere: float = 6371200.0,
        gridtype: str = "regular",
        legfunc: str = "computed",
        *,
        dtype: object | None = None,
    ) -> None:
        """
        Create a transform for a grid with nlat latitudes and nlon longitudes.

        dtype selects float32/complex64 or float64/complex128 arrays.
        """

    @property
    def nlon(self) -> int: ...
    @property
    def nlat(self) -> int: ...
    @property
    def rsphere(self) -> float: ...
    @property
    def gridtype(self) -> str: ...
    @property
    def legfunc(self) -> str: ...
    @property
    def dtype(self) -> str: ...
    def __repr__(self) -> str: ...
    def grdtospec(
        self,
        datagrid: numpy.typing.NDArray[numpy.float32]
        | numpy.typing.NDArray[numpy.float64],
        ntrunc: int | None = None,
    ) -> numpy.typing.NDArray[numpy.complex64] | numpy.typing.NDArray[numpy.complex128]:
        """
        Analyze a latitude-longitude grid into triangular spectral coefficients.
        """

    def spectogrd(
        self,
        dataspec: numpy.typing.NDArray[numpy.complex64]
        | numpy.typing.NDArray[numpy.complex128],
    ) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
        """Synthesize triangular spectral coefficients on this grid."""

    def getvrtdivspec(
        self,
        ugrid: numpy.typing.NDArray[numpy.float32]
        | numpy.typing.NDArray[numpy.float64],
        vgrid: numpy.typing.NDArray[numpy.float32]
        | numpy.typing.NDArray[numpy.float64],
        ntrunc: int | None = None,
    ) -> tuple[
        numpy.typing.NDArray[numpy.complex64] | numpy.typing.NDArray[numpy.complex128],
        numpy.typing.NDArray[numpy.complex64] | numpy.typing.NDArray[numpy.complex128],
    ]:
        """
        Analyze zonal and meridional wind components into vorticity and divergence coefficients. DUCC spin-1 analysis excludes the even-longitude Nyquist order.
        """

    def getuv(
        self,
        vrtspec: numpy.typing.NDArray[numpy.complex64]
        | numpy.typing.NDArray[numpy.complex128],
        divspec: numpy.typing.NDArray[numpy.complex64]
        | numpy.typing.NDArray[numpy.complex128],
    ) -> tuple[
        numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
        numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
    ]:
        """
        Synthesize zonal and meridional winds from vorticity and divergence coefficients.
        """

    def getgrad(
        self,
        chispec: numpy.typing.NDArray[numpy.complex64]
        | numpy.typing.NDArray[numpy.complex128],
    ) -> tuple[
        numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
        numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
    ]:
        """
        Synthesize the zonal and meridional gradient components of a scalar potential.
        """

    def getpsichi(
        self,
        ugrid: numpy.typing.NDArray[numpy.float32]
        | numpy.typing.NDArray[numpy.float64],
        vgrid: numpy.typing.NDArray[numpy.float32]
        | numpy.typing.NDArray[numpy.float64],
        ntrunc: int | None = None,
    ) -> tuple[
        numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
        numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
    ]:
        """
        Return streamfunction and velocity-potential maps for a wind field. The analysis band excludes the even-longitude Nyquist order.
        """

    def specsmooth(
        self,
        datagrid: numpy.typing.NDArray[numpy.float32]
        | numpy.typing.NDArray[numpy.float64],
        smooth: numpy.typing.NDArray[numpy.float32]
        | numpy.typing.NDArray[numpy.float64],
    ) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
        """Apply a degree-dependent spectral smoothing window to a grid field."""

def gaussian_lats_wts(
    nlat: int,
) -> tuple[numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64]]:
    """Return Gaussian latitudes in degrees and their quadrature weights."""

def getspecindx(
    ntrunc: int,
) -> tuple[numpy.typing.NDArray[numpy.intp], numpy.typing.NDArray[numpy.intp]]:
    """
    Return zonal-wavenumber and degree indices in triangular coefficient order.
    """

def legendre(
    lat: float, ntrunc: int, *, dtype: object | None = None
) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
    """Return normalized associated Legendre functions at one latitude."""

def legendre_many(
    latitudes: numpy.typing.NDArray[numpy.float32]
    | numpy.typing.NDArray[numpy.float64],
    ntrunc: int,
    *,
    dtype: object | None = None,
) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
    """
    Return normalized associated Legendre functions for a batch of latitudes.
    """

def getgeodesicpts(
    m: int, *, dtype: object | None = None
) -> tuple[
    numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
    numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
]:
    """
    Return latitude and longitude coordinates of an icosahedral grid. Defaults to float64; dtype=float32 selects the legacy-precision calculation.
    """

def specintrp(
    lon: float,
    dataspec: numpy.typing.NDArray[numpy.complex64]
    | numpy.typing.NDArray[numpy.complex128],
    legfuncs: numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
) -> float:
    """
    Evaluate spectral coefficients at one longitude using precomputed Legendre values.
    """

def specintrp_many(
    longitudes: numpy.typing.NDArray[numpy.float32]
    | numpy.typing.NDArray[numpy.float64],
    dataspec: numpy.typing.NDArray[numpy.complex64]
    | numpy.typing.NDArray[numpy.complex128],
    legfuncs: numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
) -> numpy.typing.NDArray[numpy.float64]:
    """
    Evaluate spectral coefficients at multiple longitudes using precomputed Legendre values.
    """

def regrid(
    grdin: Spharmt,
    grdout: Spharmt,
    datagrid: numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
    ntrunc: int | None = None,
    smooth: numpy.typing.NDArray[numpy.float32]
    | numpy.typing.NDArray[numpy.float64]
    | None = None,
) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
    """Analyze a field on one grid and synthesize it on another."""

def cpu_info() -> dict[str, object]:
    """
    Return information about the current CPU and DUCC dispatch configuration.

    The ``features`` entry lists SIMD/vector features recognized by DUCC that are available to the current process.

    For multiarch builds, the result also reports the compiled and available x86-64 profiles, the configured profile limit, and the profile selected for the current process.
    """
