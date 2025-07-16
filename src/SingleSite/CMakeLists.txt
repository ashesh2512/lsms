
target_sources(
  lsmscore
  PUBLIC RadialSrelSolver.cpp
         RadialSrelSolver.hpp
         SingleSiteScattering.cpp
         readSingleAtomData_hdf5.cpp
         readSingleAtomData_bigcell.cpp
         F_readSingleAtomData_bigcell.f90
         F_readSingleAtomData_bigcell_sp.f90
         single_scatterer_nonrel.f
         single_scatterer_nonrel_sp.f
         semrel_c.f
         semrel_c_sp.f
         scalar_m.f
         scalar_m_sp.f
         single_scatterer_rel.f
         single_scatterer_rel_sp.f
         spzwafu.f
         spzwafu_sp.f
         csbf.f
         csbf_sp.f
         matops.f
         matops_sp.f
         gjinv.f
         gjinv_sp.f
         dirmag1-op.f
         dirmag1-op_sp.f
         dirmag2-op.f
         dirmag2-op_sp.f
         brmat.f
         brmat_sp.f
         writeSingleAtomData_hdf5.cpp
         writeSingleAtomData_bigcell.cpp
         F_writeSingleAtomData_bigcell.f90
         F_writeSingleAtomData_bigcell_sp.f90
         checkAntiFerromagneticStatus.cpp
         AtomData.hpp
         checkAntiFerromagneticStatus.hpp
         readSingleAtomData.hpp
         SingleSiteScattering.hpp
         writeSingleAtomData.hpp)

target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
