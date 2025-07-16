

target_sources(lsmscore
        PUBLIC
        atomic_dft.cpp
        atomic_dft.hpp
        calculateCoreStates.hpp
        calculateCoreStates.cpp
        coreSolver.cpp
        coreSolver.hpp
        corslv_c.f
        corslv_c_sp.f
        deepst_c.f
        deepst_c_sp.f
        getcor_c.f
        getcor_c_sp.f
        invals_c.f
        invals_c_sp.f
        inwhnk_c.f
        inwhnk_c_sp.f
        inws_c.f
        inws_c_sp.f
        outws_c.f
        outws_c_sp.f
        richnk_c.f
        richnk_c_sp.f
        semcst_c.f
        semcst_c_sp.f
        radialSolver.cpp
        radialSolver.hpp
        States.cpp
        States.hpp
        )

target_include_directories(lsmscore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
