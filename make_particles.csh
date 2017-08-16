#!/bin/csh

set root=/lustre/lysun/models/MOM6-test/MOM6-examples/
set root_src=${root}/build/intel/shared/

set template=${root}/src/mkmf/templates/dt2-intel.mk

rm -f path_names
${root}/src/mkmf/bin/list_paths ./
${root}/src/mkmf/bin/mkmf -t ${template} -o '-I../MOM6-examples/build/intel/shared/repro' -p particles -l '-L../MOM6-examples/build/intel/shared/repro -lfms' -c "-Duse_libMPI -Duse_netCDF -DSPMD" path_names

source ${root}/build/intel/env.dt2; make NETCDF=3 REPOR=1 particles -j
