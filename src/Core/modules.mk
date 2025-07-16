
local_dir := src/Core

local_src := $(addprefix $(local_dir)/,\
        calculateCoreStates.cpp \
        coreSolver.cpp \
        corslv_c.f \
				corslv_c_sp.f \
        deepst_c.f \
				deepst_c_sp.f \
        getcor_c.f \
				getcor_c_sp.f \
        invals_c.f \
				invals_c_sp.f \
        inwhnk_c.f \
				inwhnk_c_sp.f \
        inws_c.f \
				inws_c_sp.f \
        outws_c.f \
				outws_c_sp.f \
        richnk_c.f \
				richnk_c_sp.f \
        semcst_c.f \
				semcst_c_sp.f)

sources += $(local_src)