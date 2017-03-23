PROGRAM particle_driver

  !use fms_mod, only: read_data
  use netcdf
  
  use time_manager_mod, only: time_type
  use particles_mod, only : particles_run, particles_save_restart
  !use particles_io, only : read_restart_particles
  use particles_framework, only : particle, particles, particles_gridded
  use particles_extra, only : particles_framework_ini, read_restart_particles 

  implicit none

  type(particles), target  :: drifters
  type(particles), pointer :: parts
  type(particles_gridded),pointer :: grd
  type(time_type) :: time

  CHARACTER(LEN=38) :: filename='19000101.ocean_hourly_1900_01_01_06.nc' ! Modify this later
  CHARACTER(LEN=14)  :: topname='ocean_topog.nc'
  CHARACTER(LEN=15) :: pafilename='drifters_inp.nc' 
  
  parts=>drifters

  !call mpp_domains_init

  ! create particles framework
  call particles_framework_ini(filename,topname,parts)
  grd=>parts%grd
  
  ! create drifters
  call read_restart_particles(pafilename,parts) ! Modify from: drifters_input.f90: drifters_input_new

  call particles_run(parts,time,grd%uo,grd%vo) ! Run the particles model

  call particles_save_restart(parts)

END PROGRAM particle_driver
