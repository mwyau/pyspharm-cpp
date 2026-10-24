// SPDX-License-Identifier: GPL-2.0-or-later
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>

#include "internal.hpp"

namespace SPHARM_CPP_NAMESPACE
{

using namespace nanobind::literals;

std::size_t checked_truncation(const py::object& arg, std::size_t fallback,
                               std::size_t max_value, const char* name)
{
  if (arg.is_none())
    return fallback;
  auto value = py::cast<std::int64_t>(arg);
  if (value < 0 || static_cast<std::uint64_t>(value) > max_value)
    throw py::value_error(
        (std::string(name) + " must be between 0 and " + std::to_string(max_value))
            .c_str());
  return static_cast<std::size_t>(value);
}

std::size_t truncation_from_count(std::size_t count)
{
  // Solve count=(L+1)(L+2)/2 and reject non-triangular counts exactly.
  std::size_t lo = 0, hi = 1;
  while (hi < count && hi <= std::numeric_limits<std::size_t>::max() / 2) hi *= 2;
  while (lo < hi)
  {
    const auto mid = lo + (hi - lo) / 2;
    const auto n = mid + 1;
    auto a = n, b = n + 1;
    if (a % 2 == 0)
      a /= 2;
    else
      b /= 2;
    const auto triangular = a > std::numeric_limits<std::size_t>::max() / b
                                ? std::numeric_limits<std::size_t>::max()
                                : a * b;
    if (triangular < count)
      lo = mid + 1;
    else
      hi = mid;
  }
  const auto n = lo + 1;
  if (n * (n + 1) / 2 != count)
    throw py::value_error("spectral coefficient count must be triangular");
  return lo;
}

std::size_t coefficient_count(std::size_t lmax)
{
  if (lmax > std::numeric_limits<std::size_t>::max() - 2)
    throw py::value_error("ntrunc is too large");
  auto a = lmax + 1, b = lmax + 2;
  if (a % 2 == 0)
    a /= 2;
  else
    b /= 2;
  if (a > std::numeric_limits<std::size_t>::max() / b)
    throw py::value_error("ntrunc is too large");
  return a * b;
}

std::vector<double> real_vector(const CNpArr& array, const char* name)
{
  if (array.ndim() != 1)
    throw py::value_error((std::string(name) + " must be rank one").c_str());
  std::vector<double> result(array.shape(0));
  if (DUCC0_NAMESPACE::detail_pybind::isPyarr<float>(array))
  {
    auto values = DUCC0_NAMESPACE::detail_pybind::to_cfmav<float>(array);
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = values(i);
    return result;
  }
  if (DUCC0_NAMESPACE::detail_pybind::isPyarr<double>(array))
  {
    auto values = DUCC0_NAMESPACE::detail_pybind::to_cfmav<double>(array);
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = values(i);
    return result;
  }
  throw py::type_error(
      (std::string(name) + " must have dtype numpy.float32 or numpy.float64").c_str());
}

namespace
{

template <typename T>
std::vector<std::size_t> map_shape(const CNpArr& map, std::size_t nlat,
                                   std::size_t nlon)
{
  if (map.ndim() != 2 && map.ndim() != 3)
    throw py::value_error("grid data must have rank 2 or 3");
  if (map.shape(0) != nlat || map.shape(1) != nlon)
    throw py::value_error("grid data must have shape (nlat, nlon[, nfields])");
  if (map.ndim() == 3 && map.shape(2) == 0)
    throw py::value_error("the number of fields must be positive");
  return {map.ndim() == 2 ? 1 : map.shape(2)};
}

}  // namespace

class Spharmt
{
  std::size_t nlon_;
  std::size_t nlat_;
  double rsphere_;
  std::string gridtype_;
  std::string legfunc_;
  bool float32_;

  std::size_t default_analysis_lmax() const
  {
    return std::min(gridtype_ == "regular" ? nlat_ - 2 : nlat_ - 1, nlon_ / 2);
  }

  std::size_t default_vector_analysis_lmax() const
  {
    return std::min(gridtype_ == "regular" ? nlat_ - 2 : nlat_ - 1, (nlon_ - 1) / 2);
  }

  void check_analysis_lmax(std::size_t lmax) const
  {
    if (lmax > nlon_ / 2)
      throw py::value_error(
          "analysis requires ntrunc <= floor(nlon/2) for longitude resolution");
    if (gridtype_ == "regular" && lmax > nlat_ - 2)
      throw py::value_error(
          "regular-grid analysis requires ntrunc <= nlat-2; the highest degree is not "
          "uniquely recoverable");
  }

  void check_vector_analysis_lmax(std::size_t lmax) const
  {
    if (lmax > (nlon_ - 1) / 2)
      throw py::value_error(
          "spin-1 vector analysis requires ntrunc <= floor((nlon-1)/2); "
          "DUCC analysis does not support the longitude Nyquist order");
    check_analysis_lmax(lmax);
  }

 public:
  Spharmt(std::size_t nlon, std::size_t nlat, double rsphere, std::string gridtype,
          std::string legfunc, py::object dtype)
      : nlon_(nlon),
        nlat_(nlat),
        rsphere_(rsphere),
        gridtype_(std::move(gridtype)),
        legfunc_(std::move(legfunc)),
        float32_(false)
  {
    if (nlon_ < 4)
      throw py::value_error("nlon must be at least 4");
    if (nlat_ < 3)
      throw py::value_error("nlat must be at least 3");
    if (!std::isfinite(rsphere_) || !(rsphere_ > 0.0))
      throw py::value_error("rsphere must be finite and positive");
    if (gridtype_ != "regular" && gridtype_ != "gaussian")
      throw py::value_error("gridtype must be 'regular' or 'gaussian'");
    if (legfunc_ != "computed")
      throw py::value_error("only legfunc='computed' is supported");
    if (dtype.is_none())
      dtype = py::module_::import_("numpy").attr("float64");
    auto dtype_name = py::cast<std::string>(
        py::module_::import_("numpy").attr("dtype")(dtype).attr("name"));
    if (dtype_name == "float32")
      float32_ = true;
    else if (dtype_name != "float64")
      throw py::type_error("dtype must be numpy.float32 or numpy.float64");
  }

  std::size_t nlon() const
  {
    return nlon_;
  }
  std::size_t nlat() const
  {
    return nlat_;
  }
  double rsphere() const
  {
    return rsphere_;
  }
  const std::string& gridtype() const
  {
    return gridtype_;
  }
  const std::string& legfunc() const
  {
    return legfunc_;
  }
  const char* dtype_name() const
  {
    return float32_ ? "float32" : "float64";
  }
  bool is_float32() const
  {
    return float32_;
  }
  std::size_t analysis_lmax() const
  {
    return default_analysis_lmax();
  }

  std::string repr() const
  {
    return "Spharmt(" + std::to_string(nlon_) + ", " + std::to_string(nlat_) + ", " +
           std::to_string(rsphere_) + ", '" + gridtype_ + "', '" + legfunc_ +
           "', dtype=" + dtype_name() + ")";
  }

  NpArr grdtospec(const CNpArr& grid, const py::object& ntrunc) const
  {
    map_shape<float>(grid, nlat_, nlon_);
    const auto fallback = default_analysis_lmax();
    const auto lmax = checked_truncation(ntrunc, fallback, nlat_ - 1, "ntrunc");
    check_analysis_lmax(lmax);
    if (float32_)
    {
      validate_array_dtype<float>(grid, "numpy.float32");
      return grid_to_spectral_f32(grid, nlat_, nlon_, gridtype_, lmax);
    }
    validate_array_dtype<double>(grid, "numpy.float64");
    return grid_to_spectral_f64(grid, nlat_, nlon_, gridtype_, lmax);
  }

  NpArr spectogrd(const CNpArr& spec) const
  {
    if (spec.ndim() != 1 && spec.ndim() != 2)
      throw py::value_error("spectogrd requires a rank 1 or 2 coefficient array");
    const std::size_t nalm = spec.shape(0);
    const std::size_t lmax = truncation_from_count(nalm);
    if (lmax > nlat_ - 1)
      throw py::value_error("ntrunc too large for this grid");
    if (spec.ndim() == 2 && spec.shape(1) == 0)
      throw py::value_error("the number of fields must be positive");
    if (float32_)
    {
      validate_array_dtype<std::complex<float>>(spec, "numpy.complex64");
      return spectral_to_grid_f32(spec, nlat_, nlon_, gridtype_, lmax);
    }
    validate_array_dtype<std::complex<double>>(spec, "numpy.complex128");
    return spectral_to_grid_f64(spec, nlat_, nlon_, gridtype_, lmax);
  }

  py::tuple getvrtdivspec(const CNpArr& ugrid, const CNpArr& vgrid,
                          const py::object& ntrunc) const
  {
    if (ugrid.ndim() != 2 && ugrid.ndim() != 3)
      throw py::value_error("getvrtdivspec requires rank 2 or 3 wind arrays");
    if (ugrid.ndim() != vgrid.ndim())
      throw py::value_error("ugrid and vgrid must have the same shape");
    for (std::size_t i = 0; i < ugrid.ndim(); ++i)
      if (ugrid.shape(i) != vgrid.shape(i))
        throw py::value_error("ugrid and vgrid must have the same shape");
    map_shape<float>(ugrid, nlat_, nlon_);
    const auto lmax =
        checked_truncation(ntrunc, default_vector_analysis_lmax(), nlat_ - 1, "ntrunc");
    check_vector_analysis_lmax(lmax);
    if (float32_)
    {
      validate_array_dtype<float>(ugrid, "numpy.float32");
      validate_array_dtype<float>(vgrid, "numpy.float32");
      auto result =
          analyze_winds_f32(ugrid, vgrid, nlat_, nlon_, gridtype_, rsphere_, lmax);
      return py::make_tuple(result.first, result.second);
    }
    validate_array_dtype<double>(ugrid, "numpy.float64");
    validate_array_dtype<double>(vgrid, "numpy.float64");
    auto result =
        analyze_winds_f64(ugrid, vgrid, nlat_, nlon_, gridtype_, rsphere_, lmax);
    return py::make_tuple(result.first, result.second);
  }

  py::tuple getuv(const CNpArr& vrtspec, const CNpArr& divspec) const
  {
    if ((vrtspec.ndim() != 1 && vrtspec.ndim() != 2) ||
        vrtspec.ndim() != divspec.ndim())
      throw py::value_error(
          "getuv requires rank 1 or 2 coefficient arrays of matching shape");
    if (vrtspec.shape(0) != divspec.shape(0) ||
        (vrtspec.ndim() == 2 && vrtspec.shape(1) != divspec.shape(1)))
      throw py::value_error("vrtspec and divspec must have the same shape");
    const auto nalm = vrtspec.shape(0);
    const auto lmax = truncation_from_count(nalm);
    if (lmax > nlat_ - 1)
      throw py::value_error("ntrunc too large for this grid");
    if (vrtspec.ndim() == 2 && vrtspec.shape(1) == 0)
      throw py::value_error("the number of fields must be positive");
    if (float32_)
    {
      validate_array_dtype<std::complex<float>>(vrtspec, "numpy.complex64");
      validate_array_dtype<std::complex<float>>(divspec, "numpy.complex64");
      auto result =
          synthesize_winds_f32(vrtspec, divspec, nlat_, nlon_, gridtype_, rsphere_);
      return py::make_tuple(result.first, result.second);
    }
    validate_array_dtype<std::complex<double>>(vrtspec, "numpy.complex128");
    validate_array_dtype<std::complex<double>>(divspec, "numpy.complex128");
    auto result =
        synthesize_winds_f64(vrtspec, divspec, nlat_, nlon_, gridtype_, rsphere_);
    return py::make_tuple(result.first, result.second);
  }

  py::tuple getgrad(const CNpArr& chispec) const
  {
    if (chispec.ndim() != 1 && chispec.ndim() != 2)
      throw py::value_error("getgrad requires a rank 1 or 2 coefficient array");
    const auto lmax = truncation_from_count(chispec.shape(0));
    if (lmax > nlat_ - 1)
      throw py::value_error("ntrunc too large for this grid");
    if (chispec.ndim() == 2 && chispec.shape(1) == 0)
      throw py::value_error("the number of fields must be positive");
    if (float32_)
    {
      validate_array_dtype<std::complex<float>>(chispec, "numpy.complex64");
      auto result = spectral_gradient_f32(chispec, nlat_, nlon_, gridtype_, rsphere_);
      return py::make_tuple(result.first, result.second);
    }
    validate_array_dtype<std::complex<double>>(chispec, "numpy.complex128");
    auto result = spectral_gradient_f64(chispec, nlat_, nlon_, gridtype_, rsphere_);
    return py::make_tuple(result.first, result.second);
  }

  py::tuple getpsichi(const CNpArr& ugrid, const CNpArr& vgrid,
                      const py::object& ntrunc) const
  {
    const auto lmax =
        checked_truncation(ntrunc, default_vector_analysis_lmax(), nlat_ - 1, "ntrunc");
    check_vector_analysis_lmax(lmax);
    auto components = getvrtdivspec(ugrid, vgrid, py::int_(lmax));
    auto vrt = py::cast<CNpArr>(components[0]);
    auto div = py::cast<CNpArr>(components[1]);
    if (float32_)
    {
      NpArr psi = apply_inverse_laplacian_f32(vrt, rsphere_);
      NpArr chi = apply_inverse_laplacian_f32(div, rsphere_);
      return py::make_tuple(spectogrd(CNpArr(psi)), spectogrd(CNpArr(chi)));
    }
    NpArr psi = apply_inverse_laplacian_f64(vrt, rsphere_);
    NpArr chi = apply_inverse_laplacian_f64(div, rsphere_);
    return py::make_tuple(spectogrd(CNpArr(psi)), spectogrd(CNpArr(chi)));
  }

  NpArr specsmooth(const CNpArr& grid, const CNpArr& smooth) const
  {
    map_shape<float>(grid, nlat_, nlon_);
    if (smooth.ndim() != 1 || smooth.shape(0) != nlat_)
      throw py::value_error("smooth must be rank 1 with length nlat");
    CNpArr spec(grdtospec(grid, py::none()));
    if (float32_)
    {
      validate_array_dtype<float>(smooth, "numpy.float32");
      NpArr filtered = apply_smoothing_f32(spec, smooth);
      return spectogrd(CNpArr(filtered));
    }
    validate_array_dtype<double>(smooth, "numpy.float64");
    NpArr filtered = apply_smoothing_f64(spec, smooth);
    return spectogrd(CNpArr(filtered));
  }
};

void add_api(py::module_& m)
{
  py::class_<Spharmt>(
      m, "Spharmt",
      R"doc(Scalar and vector spherical harmonic transforms on latitude-longitude grids.

Grid maps are ordered north to south and eastward in longitude.
Input fields have shape (nlat, nlon) or (nlat, nlon, nfields).
Spectral coefficients have shape (ncoeff,) or (ncoeff, nfields)
and use triangular degree and order packing.

See Numerical conventions for normalization, coordinate conventions,
and the separate analysis and synthesis truncation limits.)doc",
      py::dynamic_attr())
      .def(py::init<std::size_t, std::size_t, double, std::string, std::string,
                    py::object>(),
           "nlon"_a, "nlat"_a, "rsphere"_a = 6.3712e6, "gridtype"_a = "regular",
           "legfunc"_a = "computed", py::kw_only(), "dtype"_a = py::none(),
           py::sig("def __init__(self, nlon: int, nlat: int, rsphere: float = "
                   "6371200.0, gridtype: str = 'regular', legfunc: str = 'computed', "
                   "*, dtype: object | None = None) -> None"),
           R"doc(Create a transform on a regular or Gaussian latitude-longitude grid.

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
truncation ntrunc, there are (ntrunc + 1) * (ntrunc + 2) / 2 coefficients.)doc")
      .def_prop_ro("nlon", &Spharmt::nlon)
      .def_prop_ro("nlat", &Spharmt::nlat)
      .def_prop_ro("rsphere", &Spharmt::rsphere)
      .def_prop_ro("gridtype", &Spharmt::gridtype)
      .def_prop_ro("legfunc", &Spharmt::legfunc)
      .def_prop_ro("dtype", &Spharmt::dtype_name)
      .def("__repr__", &Spharmt::repr)
      .def("grdtospec", &Spharmt::grdtospec, "datagrid"_a, "ntrunc"_a = py::none(),
           py::sig("def grdtospec(self, datagrid: numpy.typing.NDArray[numpy.float32] "
                   "| numpy.typing.NDArray[numpy.float64], ntrunc: int | None = None) "
                   "-> numpy.typing.NDArray[numpy.complex64] | "
                   "numpy.typing.NDArray[numpy.complex128]"),
           R"doc(Analyze gridded scalar fields into spherical harmonic coefficients.

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
Gaussian analysis permits ntrunc <= min(nlat - 1, nlon // 2).)doc")
      .def(
          "spectogrd", &Spharmt::spectogrd, "dataspec"_a,
          py::sig("def spectogrd(self, dataspec: numpy.typing.NDArray[numpy.complex64] "
                  "| numpy.typing.NDArray[numpy.complex128]) -> "
                  "numpy.typing.NDArray[numpy.float32] | "
                  "numpy.typing.NDArray[numpy.float64]"),
          R"doc(Synthesize scalar fields from triangular spherical harmonic coefficients.

Parameters
----------
dataspec : ndarray
    Complex coefficients of shape (ncoeff,) or (ncoeff, nfields).
    The first dimension determines the harmonic truncation.

Returns
-------
ndarray
    Real grid fields with shape (nlat, nlon) or
    (nlat, nlon, nfields), in the transform's precision.)doc")
      .def(
          "getvrtdivspec", &Spharmt::getvrtdivspec, "ugrid"_a, "vgrid"_a,
          "ntrunc"_a = py::none(),
          py::sig("def getvrtdivspec(self, ugrid: numpy.typing.NDArray[numpy.float32] "
                  "| numpy.typing.NDArray[numpy.float64], vgrid: "
                  "numpy.typing.NDArray[numpy.float32] | "
                  "numpy.typing.NDArray[numpy.float64], ntrunc: int | None = None) -> "
                  "tuple[numpy.typing.NDArray[numpy.complex64] | "
                  "numpy.typing.NDArray[numpy.complex128], "
                  "numpy.typing.NDArray[numpy.complex64] | "
                  "numpy.typing.NDArray[numpy.complex128]]"),
          R"doc(Analyze eastward and northward winds into relative vorticity and divergence.

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
order, and the corresponding grid's latitude truncation limit.)doc")
      .def("getuv", &Spharmt::getuv, "vrtspec"_a, "divspec"_a,
           py::sig("def getuv(self, vrtspec: numpy.typing.NDArray[numpy.complex64] | "
                   "numpy.typing.NDArray[numpy.complex128], divspec: "
                   "numpy.typing.NDArray[numpy.complex64] | "
                   "numpy.typing.NDArray[numpy.complex128]) -> "
                   "tuple[numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64], "
                   "numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64]]"),
           R"doc(Synthesize eastward and northward winds from vorticity and divergence.

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
arrays of its selected precision.)doc")
      .def("getgrad", &Spharmt::getgrad, "chispec"_a,
           py::sig("def getgrad(self, chispec: numpy.typing.NDArray[numpy.complex64] | "
                   "numpy.typing.NDArray[numpy.complex128]) -> "
                   "tuple[numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64], "
                   "numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64]]"),
           R"doc(Compute the eastward and northward gradient of a scalar field.

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
Each component is a spatial derivative divided by the sphere radius.)doc")
      .def("getpsichi", &Spharmt::getpsichi, "ugrid"_a, "vgrid"_a,
           "ntrunc"_a = py::none(),
           py::sig("def getpsichi(self, ugrid: numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64], vgrid: "
                   "numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64], ntrunc: int | None = None) -> "
                   "tuple[numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64], "
                   "numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64]]"),
           R"doc(Compute streamfunction and velocity potential from wind components.

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
order.)doc")
      .def("specsmooth", &Spharmt::specsmooth, "datagrid"_a, "smooth"_a,
           py::sig("def specsmooth(self, datagrid: numpy.typing.NDArray[numpy.float32] "
                   "| numpy.typing.NDArray[numpy.float64], smooth: "
                   "numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64]) -> "
                   "numpy.typing.NDArray[numpy.float32] | "
                   "numpy.typing.NDArray[numpy.float64]"),
           R"doc(Apply an isotropic spectral filter to one or more scalar maps.

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
    Filtered map or maps with the input grid shape and dtype.)doc");

  bind_numerics(m);

  m.def(
      "regrid",
      [](const Spharmt& grdin, const Spharmt& grdout, const CNpArr& datagrid,
         const py::object& ntrunc, const py::object& smooth)
      {
        const auto max_input = grdin.nlat() - 1;
        const auto max_output = grdout.nlat() - 1;
        const auto default_lmax = std::min(grdin.analysis_lmax(), max_output);
        const auto lmax = checked_truncation(ntrunc, default_lmax,
                                             std::min(max_input, max_output), "ntrunc");
        NpArr spec = grdin.grdtospec(datagrid, py::int_(lmax));
        if (smooth.is_none())
        {
          if (grdin.is_float32() != grdout.is_float32())
            spec = grdin.is_float32()
                       ? convert_spectral_precision_f64_from_f32(CNpArr(spec))
                       : convert_spectral_precision_f32_from_f64(CNpArr(spec));
          return grdout.spectogrd(CNpArr(spec));
        }
        CNpArr factors = py::cast<CNpArr>(smooth);
        if (factors.ndim() != 1 || factors.shape(0) != grdout.nlat())
          throw py::value_error("smooth must be rank 1 with length grdout.nlat");
        if (grdout.is_float32())
        {
          if (grdin.is_float32())
            spec = apply_smoothing_f32(CNpArr(spec), factors);
          else
          {
            spec = convert_spectral_precision_f32_from_f64(CNpArr(spec));
            spec = apply_smoothing_f32(CNpArr(spec), factors);
          }
        }
        else
        {
          if (grdin.is_float32())
            spec = convert_spectral_precision_f64_from_f32(CNpArr(spec));
          spec = apply_smoothing_f64(CNpArr(spec), factors);
        }
        return grdout.spectogrd(CNpArr(spec));
      },
      "grdin"_a, "grdout"_a, "datagrid"_a, "ntrunc"_a = py::none(),
      "smooth"_a = py::none(),
      py::sig(
          "def regrid(grdin: Spharmt, grdout: Spharmt, datagrid: "
          "numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64], "
          "ntrunc: int | None = None, smooth: numpy.typing.NDArray[numpy.float32] | "
          "numpy.typing.NDArray[numpy.float64] | None = None) -> "
          "numpy.typing.NDArray[numpy.float32] | numpy.typing.NDArray[numpy.float64]"),
      R"doc(Regrid scalar fields by spherical harmonic analysis and synthesis.

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
    The output dtype follows grdout.)doc");
}

}  // namespace SPHARM_CPP_NAMESPACE
