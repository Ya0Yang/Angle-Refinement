cd src/accelerate
mex accel_cross.cpp CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3'
mex accel_get_commonline_vectors.cpp CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3'
mex accel_matrix_33_31.cpp CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3'
mex accel_matrix_33_31_tmult.cpp CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3'
cd splinterp
compile_mex_script
cd ../../../