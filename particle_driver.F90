PROGRAM particle_driver

  !use fms_mod, only: read_data
  use netcdf
  use MOM_file_parser, only : get_param, log_param, param_file_type
  use MOM_get_input, only : get_MOM_input, directories
  use MOM_domains, only : MOM_domains_init
  use MOM_grid, only : ocean_grid_type, MOM_grid_init
  use MOM_io, only : slasher, MOM_io_init
  use MOM_hor_index, only : hor_index_type, hor_index_init

  use time_manager_mod, only: time_type
  use particles_mod, only : particles_run, particles_save_restart
  !use particles_io, only : read_restart_particles
  use particles_framework, only : particle, particles, particles_gridded
  use particles_extra, only : particles_framework_ini, read_restart_particles 

  implicit none

  type(ocean_grid_type), target :: Grid
  type(hor_index_type)   :: HI ! A hor_index_type for array extents
  type(particles), target  :: drifters
  type(particles), pointer :: parts
  type(particles_gridded),pointer :: grd
  type(time_type) :: time

  CHARACTER(LEN=38) :: filename='19000101.ocean_hourly_1900_01_01_06.nc' ! Modify this later
  CHARACTER(LEN=14)  :: topname='ocean_topog.nc'
  CHARACTER(LEN=15) :: pafilename='drifters_inp.nc' 
  character(len=40) :: mod = 'MOM'
  type(param_file_type) :: PF
  type(directories) :: dirs

  call get_MOM_input(PF,dirs)
  call MOM_domains_init(Grid%domain,PF)
  call MOM_io_init(PF)
  call hor_index_init(Grid%Domain, HI, PF, &
                      local_indexing=.true.)
!  call verticalGridInit( param_file, CS%GV )
! Set up the parameters of the physical domain (i.e. the grid), G
!  call set_grid_metrics(G, PF)

! Set up the bottom depth, G%bathyT either analytically or from file
! This also sets G%max_depth based on the input parameter MAXIMUM_DEPTH,
! or, if absent, is diagnosed as G%max_depth = max( G%D(:,:) )
!  call MOM_initialize_topography(G%bathyT, G%max_depth, G, PF)

!  parts=>drifters

  !call mpp_domains_init

  ! create particles framework
!  call particles_framework_ini(filename,topname,parts)
!  grd=>parts%grd
  
  ! create drifters
!  call read_restart_particles(pafilename,parts) ! Modify from: drifters_input.f90: drifters_input_new

!  call particles_run(parts,time,grd%uo,grd%vo) ! Run the particles model

!  call particles_save_restart(parts)

END PROGRAM particle_driver
