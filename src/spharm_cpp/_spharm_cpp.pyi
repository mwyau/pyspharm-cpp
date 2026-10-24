import numpy
import numpy.typing

# SPDX-License-Identifier: GPL-2.0-or-later

class Spharmt:
    """
    Scalar and vector spherical harmonic transforms on latitude-longitude grids.

    Grid maps are ordered north to south and eastward in longitude.
    Input fields have shape (nlat, nlon) or (nlat, nlon, nfields).
    Spectral coefficients have shape (ncoeff,) or (ncoeff, nfields)
    and use triangular degree and order packing.

    See Numerical conventions for normalization, coordinate conventions,
    and the separate analysis and synthesis truncation limits.
    """

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
        Create a transform on a regular or Gaussian latitude-longitude grid.

        Parameters
        ----------
        nlon : int
            Number of equally spaced longitudes, at least 4.
        nlat : int
            Number of latitude rings, at least 3.
        rsphere : float, default 6371200.0
            Positive sphere radius. Derivatives and wind quantities use this
            radius in the chosen length units.
        gridtype : {'regular', 'gaussian'}, default 'regular'
            Clenshaw-Curtis latitudes including both poles, or Gauss-Legendre
            latitudes excluding the poles.
        legfunc : {'computed'}, default 'computed'
            Associated Legendre functions are computed during each transform.
        dtype : numpy dtype, optional
            float64 (default) or float32 for maps; spectral coefficients
            use the corresponding complex128 or complex64 precision.

        Notes
        -----
        The public coefficient ordering follows PySpharm. For an explicit
        truncation ntrunc, there are (ntrunc + 1) * (ntrunc + 2) / 2 coefficients.
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
        Analyze gridded scalar fields into spherical harmonic coefficients.

        Parameters
        ----------
        datagrid : ndarray
            Real-valued array with shape (nlat, nlon) or (nlat, nlon, nfields).
            Its dtype must match the transform.
        ntrunc : int, optional
            Maximum spherical harmonic degree. The default is the grid's
            maximum analysis truncation.

        Returns
        -------
        ndarray
            Complex coefficients with shape (ncoeff,) or (ncoeff, nfields).
            Order m is the outer packed index; degree l runs from m to ntrunc.

        Notes
        -----
        Clenshaw-Curtis analysis permits ntrunc <= min(nlat - 2, nlon // 2).
        Gaussian analysis permits ntrunc <= min(nlat - 1, nlon // 2).
        """

    def spectogrd(
        self,
        dataspec: numpy.typing.NDArray[numpy.complex64]
        | numpy.typing.NDArray[numpy.complex128],
    ) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
        """
        Synthesize scalar fields from triangular spherical harmonic coefficients.

        Parameters
        ----------
        dataspec : ndarray
            Complex coefficients of shape (ncoeff,) or (ncoeff, nfields).
            The first dimension determines the harmonic truncation.

        Returns
        -------
        ndarray
            Real grid fields with shape (nlat, nlon) or
            (nlat, nlon, nfields), in the transform's precision.
        """

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
        Analyze eastward and northward winds into relative vorticity and divergence.

        Parameters
        ----------
        ugrid, vgrid : ndarray
            Eastward and northward wind components of matching shape
            (nlat, nlon) or (nlat, nlon, nfields), in the transform's dtype.
        ntrunc : int, optional
            Maximum degree of the vector analysis.

        Returns
        -------
        (vrtspec, divspec) : tuple of ndarray
            Spectral coefficients of relative vorticity and divergence,
            each with shape (ncoeff,) or (ncoeff, nfields).

        Notes
        -----
        The sphere radius determines the derivative scale. Vector analysis
        uses ntrunc <= (nlon - 1) // 2, excluding the even-longitude Nyquist
        order, and the corresponding grid's latitude truncation limit.
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
        Synthesize eastward and northward winds from vorticity and divergence.

        Parameters
        ----------
        vrtspec, divspec : ndarray
            Complex relative-vorticity and divergence coefficients with
            identical shape (ncoeff,) or (ncoeff, nfields).

        Returns
        -------
        (u, v) : tuple of ndarray
            Eastward and northward wind maps, each with shape
            (nlat, nlon) or (nlat, nlon, nfields).

        Notes
        -----
        The synthesis uses the transform's sphere radius and returns real
        arrays of its selected precision.
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
        Compute the eastward and northward gradient of a scalar field.

        Parameters
        ----------
        chispec : ndarray
            Complex scalar spectral coefficients with shape (ncoeff,)
            or (ncoeff, nfields).

        Returns
        -------
        (eastward, northward) : tuple of ndarray
            Real components of the spherical gradient on the transform grid.
            Each component has the corresponding grid shape.

        Notes
        -----
        Each component is a spatial derivative divided by the sphere radius.
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
        Compute streamfunction and velocity potential from wind components.

        Parameters
        ----------
        ugrid, vgrid : ndarray
            Eastward and northward wind maps of the same shape and dtype.
        ntrunc : int, optional
            Maximum degree used to analyze the wind.

        Returns
        -------
        (psi, chi) : tuple of ndarray
            Real streamfunction and velocity-potential maps with the input
            grid shape and precision.

        Notes
        -----
        The degree-zero inverse Laplacian coefficient is set to zero.
        The vector analysis uses degrees below the even-longitude Nyquist
        order.
        """

    def specsmooth(
        self,
        datagrid: numpy.typing.NDArray[numpy.float32]
        | numpy.typing.NDArray[numpy.float64],
        smooth: numpy.typing.NDArray[numpy.float32]
        | numpy.typing.NDArray[numpy.float64],
    ) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
        """
        Apply an isotropic spectral filter to one or more scalar maps.

        Parameters
        ----------
        datagrid : ndarray
            Real map with shape (nlat, nlon) or (nlat, nlon, nfields).
        smooth : ndarray
            Real smoothing factors of length nlat; smooth[l] multiplies
            every spherical harmonic coefficient of degree l. Its dtype
            must match the transform.

        Returns
        -------
        ndarray
            Filtered map or maps with the input grid shape and dtype.
        """

def gaussian_lats_wts(
    nlat: int,
) -> tuple[numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64]]:
    """
    Return Gauss-Legendre latitudes and quadrature weights.

    Parameters
    ----------
    nlat : int
        Positive number of Gaussian latitude rings.

    Returns
    -------
    (latitudes, weights) : tuple of ndarray
        Float64 arrays of length nlat. Latitude is in degrees and
        decreases from north to south. Weights integrate over sine
        latitude and sum to 2 up to floating-point roundoff.
    """

def getspecindx(
    ntrunc: int,
) -> tuple[numpy.typing.NDArray[numpy.intp], numpy.typing.NDArray[numpy.intp]]:
    """
    Return longitudinal orders and degrees of packed spectral coefficients.

    Parameters
    ----------
    ntrunc : int
        Nonnegative maximum spherical harmonic degree.

    Returns
    -------
    (order, degree) : tuple of ndarray
        Integer arrays of length (ntrunc + 1) * (ntrunc + 2) / 2.
        Order m is the outer index and degree l runs from m to ntrunc.
    """

def legendre(
    lat: float, ntrunc: int, *, dtype: object | None = None
) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
    """
    Evaluate normalized associated Legendre functions at one latitude.

    Parameters
    ----------
    lat : float
        Latitude in degrees, between -90 and 90.
    ntrunc : int
        Nonnegative maximum spherical harmonic degree.
    dtype : numpy dtype, optional
        Output float32 or float64 (default).

    Returns
    -------
    ndarray
        Real values in triangular coefficient order, of length
        (ntrunc + 1) * (ntrunc + 2) / 2.

    Notes
    -----
    The normalization excludes the Condon-Shortley phase. In particular,
    the degree-zero value is 1 / sqrt(2).
    """

def legendre_many(
    latitudes: numpy.typing.NDArray[numpy.float32]
    | numpy.typing.NDArray[numpy.float64],
    ntrunc: int,
    *,
    dtype: object | None = None,
) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
    """
    Evaluate normalized associated Legendre functions at several latitudes.

    Parameters
    ----------
    latitudes : ndarray
        One-dimensional real latitude array in degrees, between -90 and 90.
    ntrunc : int
        Nonnegative maximum spherical harmonic degree.
    dtype : numpy dtype, optional
        Output float32 or float64 (default).

    Returns
    -------
    ndarray
        Real array of shape (npoints, ncoeff), using the same
        normalization and packing as legendre.
    """

def getgeodesicpts(
    m: int, *, dtype: object | None = None
) -> tuple[
    numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
    numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
]:
    """
    Return the vertices of a subdivided icosahedral geodesic grid.

    Parameters
    ----------
    m : int
        Number of points per icosahedral edge, including endpoints;
        must be at least 1.
    dtype : numpy dtype, optional
        float64 (default) or float32 output coordinates.

    Returns
    -------
    (latitudes, longitudes) : tuple of ndarray
        Coordinate arrays in degrees with 10 * (m - 1)**2 + 2 points.
        Coordinates are computed in float64 before optional conversion.
    """

def specintrp(
    lon: float,
    dataspec: numpy.typing.NDArray[numpy.complex64]
    | numpy.typing.NDArray[numpy.complex128],
    legfuncs: numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
) -> float:
    """
    Evaluate a scalar spherical harmonic expansion at one point.

    Parameters
    ----------
    lon : float
        Longitude in degrees.
    dataspec : ndarray
        One-dimensional complex coefficients in triangular packing.
    legfuncs : ndarray
        One-dimensional associated Legendre values at the point's
        latitude, as returned by legendre. Precision must match
        dataspec (float32/complex64 or float64/complex128).

    Returns
    -------
    float
        Real interpolated value at the specified longitude and
        latitude represented by legfuncs.
    """

def specintrp_many(
    longitudes: numpy.typing.NDArray[numpy.float32]
    | numpy.typing.NDArray[numpy.float64],
    dataspec: numpy.typing.NDArray[numpy.complex64]
    | numpy.typing.NDArray[numpy.complex128],
    legfuncs: numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64],
) -> numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]:
    """
    Evaluate one scalar spherical harmonic expansion at several points.

    Parameters
    ----------
    longitudes : ndarray
        One-dimensional longitude array in degrees, length npoints.
    dataspec : ndarray
        One-dimensional complex coefficients in triangular packing.
    legfuncs : ndarray
        Associated Legendre values of shape (npoints, ncoeff), e.g.,
        from legendre_many. Precision must match dataspec.

    Returns
    -------
    ndarray
        Interpolated values of shape (npoints,), with float32 or float64
        precision matching dataspec and legfuncs.
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
    """
    Regrid scalar fields by spherical harmonic analysis and synthesis.

    Parameters
    ----------
    grdin, grdout : Spharmt
        Source and target transforms. They may use different grid
        types, dimensions, or floating-point precision.
    datagrid : ndarray
        Source-grid fields of shape (grdin.nlat, grdin.nlon) or
        (grdin.nlat, grdin.nlon, nfields).
    ntrunc : int, optional
        Shared harmonic truncation. The default uses the source
        analysis limit and target synthesis limit.
    smooth : ndarray, optional
        Degree-dependent smoothing factors of length grdout.nlat,
        in the target transform's real precision.

    Returns
    -------
    ndarray
        Real fields on the target grid with shape
        (grdout.nlat, grdout.nlon) or (grdout.nlat, grdout.nlon, nfields).
        The output dtype follows grdout.
    """

def cpu_info() -> dict[str, object]:
    """
    Return detected processor features and spherical transform dispatch metadata.

    Returns
    -------
    dict
        architecture and features describe the detected CPU.
        multiarch indicates whether multiple binary profiles are compiled.
        Multiarch builds also report compiled_profiles, available_profiles,
        configured_limit, and active_profile.

    Notes
    -----
    The active profile describes the current pyspharm-cpp implementation.
    """
