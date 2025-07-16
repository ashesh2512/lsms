#ifndef LSMS_REAL_H
#define LSMS_REAL_H

#include <cmath>

#include "TypeTraits.hpp"

#ifndef __CUDACC__
#include <mpi.h>
#include <hdf5.h>

#if CMAKE_SP == 1
  #define SP 1
#else
  #define SP 0
#endif

#if SP
  typedef float Real;
#else
  typedef double Real;
#endif

template <typename T>
inline Real toReal(T value) {
  return static_cast<Real>(value);
}

template <typename T>
inline T exp_t(T x) {
    return std::exp(x);  // for general type
}
template <>
inline float exp_t<float>(float x) {
    return std::expf(x);  // float-specific
}

template <typename T>
inline T pow_t(T x, T y) {
    return std::pow(x, y);  // for general type
}
template <typename T>
inline T pow_t(T x, int y) {
#if SP
  return std::powf(x, y);
#else
  return std::pow(x, y);
#endif
}

template <typename T>
inline T sqrt_t(T x) {
    return std::sqrt(x);  // for general type
}
template <>
inline float sqrt_t<float>(float x) {
    return std::sqrtf(x);  // float-specific
}


template<>
class TypeTraits<Real>
{
//  static hid_t hdf5_type;
public:
  inline static MPI_Datatype mpiType(void) {
#if SP
    return MPI_FLOAT;
#else
    return MPI_DOUBLE;
#endif
  }

  inline static hid_t hdf5Type(void) {
#if SP
    return H5T_NATIVE_FLOAT;
#else
    return H5T_NATIVE_DOUBLE;
#endif
  }
};

template<>
class TypeTraits<int>
{
//   static hid_t hdf5_type;
public:
  inline static MPI_Datatype mpiType(void) {return MPI_INT;}
  inline static hid_t hdf5Type(void) { return H5T_NATIVE_INT;}
};
// hid_t TypeTraits<int>::hdf5_type = H5Tcopy(H5T_NATIVE_INT);

template<>
class TypeTraits<char>
{
//   static hid_t hdf5_type;
public:
  inline static MPI_Datatype mpiType(void) {return MPI_CHAR;}
  inline static hid_t hdf5Type(void) { return H5T_NATIVE_CHAR;}
};
// hid_t TypeTraits<char>::hdf5_type = H5Tcopy(H5T_NATIVE_CHAR);

#endif

#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510 
#endif

