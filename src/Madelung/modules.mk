
local_dir := src/Madelung

local_src := $(addprefix $(local_dir)/,\
		 cal_madelung_matrix.f \
		 cal_madelung_matrix_sp.f \
         getkncut.f \
				 getkncut_sp.f \
         getrscut.f \
				 getrscut_sp.f \
         getstruc.f \
				 getstruc_sp.f \
         interf.f \
         interfsmr.f \
         lattice.f \
				 lattice_sp.f \
         madewd.f \
				 madewd_sp.f \
         madewdj.f \
         madsum.f \
				 madsum_sp.f \
         ord3v.f \
				 ord3v_sp.f \
         pqintg_c.f \
         lmfacts.f \
         bessj.f \
				 bessj_sp.f \
         Madelung.cpp)

sources += $(local_src)