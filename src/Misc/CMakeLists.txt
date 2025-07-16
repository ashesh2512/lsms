


target_sources(lsmscore
        PUBLIC

        SphericalHarmonics.cpp
        SphericalHarmonics.hpp

        associatedLegendreFunction.cpp
        associatedLegendreFunction.hpp
        bulirsch_stoer.f
        bulirschStoerIntegrator.cpp
        bulirschStoerIntegrator.hpp
        calculateGauntCoeficients.cpp
        calculateGaussLegendrePoints.hpp

        cgaunt_c.f90
        cgaunt_c_sp.f90

        clock_time.c
        ClebschGordan.hpp
        ClebschGordan.cpp
        cgc.h
        cinterp.f
        cinterp_sp.f
        clebsch.f
        clebsch_sp.f
        cmtruni.f
        cnewint.f
        Coeficients.cpp
        Coeficients.hpp
        congauss_c.f
        congauss_c_sp.f
        constraint.f
        constraint_sp.f
        dfv.f
        dfv_sp.f
        dfv_new.f
        dfv_new_sp.f
        fit.f
        fit_sp.f
        fitpot.f
        fitpot_sp.f
        fnpi.f
        fstop.f
        gaunt.f
        gaunt_sp.f
        Gaunt_gmp.hpp
        Gaunt.hpp
        getclm.f
        ifacts_c.f
        ifacts_c_sp.f

        Indices.hpp
        Indices.cpp

        spherical_harmonics.hpp
        spherical_harmonics.f90

        initwave.f
        initwave_sp.f
        integrateOneDim.hpp
        interp.f
        interp_sp.f
        matr.f
        matr_sp.f
        matrot1.f
        matrot1_sp.f
        mbeqa.f
        mbeqa_sp.f
        mod_midpoint.f
        newder.f
        newint.f
        newint_sp.f
        plglmax.f
        quadrature.cpp
        rationalFit.hpp
        readLastLine.cpp
        readLastLine.hpp
        ricbes.f
        ricbes_sp.f
        rotmat.f
        rotmat_sp.f
        rsimp.f
        rwave.f
        rwave_sp.f
        rzextr.f
        spin_trafo.f
        spin_trafo_sp.f
        stop_with_backtrace.hpp
        trltog.f
        trltog_sp.f
        u_sigma_u.f
        u_sigma_u_sp.f
        v_plus_minus.f
        v_plus_minus_sp.f
        wrtmtx.f
        wrtmtx_sp.f
        ylag.f
        zeroout.f
        zeroout_sp.f
        zsphbes.f
        zsphbesj.f
        zsphbesj_sp.f
        zsphbesjh.f
        zsphbesjh_sp.f
        io.hpp

        utils.hpp
        integrator.hpp
        poisson.hpp
        diff.hpp

        PotentialType.hpp
        num_digits.hpp

        ChargePlotter.cpp
        ChargePlotter.hpp

        )

if(USE_ESSL)
  target_sources(lsmscore
        PUBLIC
        essl_workaround.f
        )
  message(STATUS "Adding essl_workaround.f")
endif()

target_include_directories(
        lsmscore
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)


