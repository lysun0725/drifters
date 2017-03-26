PROGRAM particle_driver

  use MOM_file_parser, only : get_param, log_param, param_file_type
  use MOM_get_input, only : get_MOM_input, directories
  use MOM_domains, only : MOM_domains_init
  use MOM_grid, only : ocean_grid_type, MOM_grid_init
  use MOM_io, only : slasher, MOM_io_init
  use MOM_hor_index, only : hor_index_type, hor_index_init
  use time_manager_mod, only: time_type
!  use particles_mod, only : particles_run, particles_save_restart
!  use particles_io, only : read_restart_particles
!  use particles_framework, only : particle, particles, particles_gridded
!  use particles_extra, only : particles_framework_ini, read_restart_particles

  implicit none

  type(ocean_grid_type), target :: Grid
  type(ocean_vertical_grid_type), target :: GV
  type(dyn_horgrid_type), pointer :: dG => NULL()
  type(MOM_control_struct), target, dimension(:) :: CSp !<pointer of ensemble list of MOM control structures
  type(MOM_control_struct), pointer :: CS=NULL() !<pointer to an element in CSp
  type(hor_index_type)   :: HI ! A hor_index_type for array extents
  type(particles), target  :: drifters
  type(particles), pointer :: node=>NULL() !<pointer to an element in drifters list
  type(particles_gridded),pointer :: grd
  type(time_type) :: time, start_time, end_time
  real :: time_step
  character(len=128) :: history_file
  character(len=32)  :: topo_file
  character(len=32)  :: drifter_file
  character(len=8)  :: mod = 'MOM'
  type(param_file_type) :: PF
  type(directories) :: dirs
  !<ens_info(1)=ensemble_size; ens_info(2)=tot_pes
  !<ens_info(3)=ocn_pes;ens_info(4)=atm_pes
  !<ens_info(5)=land_pes;ens_info(6)=ice_pes
  integer :: ens_size, nPEs_ocn, nPEs_atm, nPEs_lnd, nPES_ice, ens_info(6)
  integer, dimension(:), allocatable :: ocn_pelist, atm_pelist, lnd_pelist, ice_pelist
  integer :: calendar_type=-1
  character(len=16) :: month='jan',calendar='julian'
  logical :: concurrent = .false. !<use concurrent PE execution of model components
  integer :: date_init(6), years=0, months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: yr, mon, day, hr, min, sec

  namelist /particle_driver_nml/ date_init, calendar, months, days, hours, minutes, seconds

  call MOM_infra_init(); call io_infra_init()
  call ensemble_manager_init(); ens_info = get_ensemble_size()
  ens_size=ens_info(1); nPEs_ocn=ens_info(2); nPEs_atm=0; nPEs_lnd=0; nPEs_ice=0

  call ensemble_pelist_setup(concurrent,nPEs_atm,nPES_ocn,nPES_lnd,npes_ice,&
       atm_pelist,ocn_pelist,lnd_pelist,ice_pelist)
  allocate(ocn_pelist(nPEs_ocn))
  call set_current_pelist(ocn_pelist)
  allocate(CS(nPEs_ocn))

  call open_file(unit,'input.nml',form=ASCII_FIE,acrion=READONLY_FILE)
  read(unit,particle_driver_nml,iostat=io_status)
  call close_file(unit)

  ierr = check_nml_error(io_status,'particle_driver_nml')

  call get_MOM_input(PF,dirs)

  if (file_exists(trim(dirs%restart_input_dir)//'ocean.res')) then
    Mak
    call open_file(unit,trim(dirs%restart_input_dir)//'ocean.res',form=ASCII_FILE)
    read(unit,*) calendar_type; read(unit,*) date_init; read(unit,*) date
    call close_file(unit)
  endif

  if (calendar(1:6)=='julian') calendar_type=JULIAN
  if (calendar(1:11)=='no_calendar') calendar_type=NO_CALENDAR
  if (calendar(1:11)=='noleap') calendar_type=NOLEAP
  if (calendar_type<0) call MOM_error(FATAL,'paricles_driver: Invalid namelist value for calendar')

  call set_calendar_type(calendar_type)

  time_start = set_time(0,days=0)
  if (sum_date_init)>0) time_start=set_date(date_init(1),date_init(2), date_init(3), &
       date_init(4), date_init(5), date_init(6))

  if (sum(date)>0) then
    time_start_segment = time_start
    time=time_start_segment
  else
    time=time_start
  endif

  call MOM_domains_init(Grid%domain,PF)
  call MOM_io_init(PF)
  call hor_index_init(Grid%Domain, HI, PF, &
                      local_indexing=.true.)
  call verticalGridInit( PF, CS%GV )
  call create_dyn_horgrid(dG,HI)
  call clone_MOM_domain(G%Domain, dG%Domain)

  is   = dG%isc   ; ie   = dG%iec  ; js   = dG%jsc  ; je   = dG%jec ; nz = GV%ke
  isd  = dG%isd   ; ied  = dG%ied  ; jsd  = dG%jsd  ; jed  = dG%jed
  IsdB = dG%IsdB  ; IedB = dG%IedB ; JsdB = dG%JsdB ; JedB = dG%JedB

  do n=1,nPEs_ocn
    CS=>CSp(n)
  ! Allocate and initialize space for primary MOM variables.
    ALLOC_(CS%u(IsdB:IedB,jsd:jed,nz))   ; CS%u(:,:,:) = 0.0
    ALLOC_(CS%v(isd:ied,JsdB:JedB,nz))   ; CS%v(:,:,:) = 0.0
    ALLOC_(CS%h(isd:ied,jsd:jed,nz))     ; CS%h(:,:,:) = GV%Angstrom
  enddo

  call MOM_initialize_fixed(G,CS%OBC,PF,.false.,dirs%output_directory)
  call MOM_initialize_coord(GV, PF, .false., &
       dirs%output_directory, CS%tv, dG%max_depth)

  !   Shift from using the temporary dynamic grid type to using the final
  ! (potentially static) ocean-specific grid type.
  !   The next line would be needed if G%Domain had not already been init'd above:
  !     call clone_MOM_domain(dG%Domain, G%Domain)
  call MOM_grid_init(G, param_file, HI, bathymetry_at_vel=bathy_at_vel)
  call copy_dyngrid_to_MOM_grid(dG, G)
  call destroy_dyn_horgrid(dG)


  call MOM_initialize_state(CS%u, CS%v, CS%h, CS%tv, Time, G, GV, PF, &
                            dirs, CS%restart_CSp, CS%ALE_CSp, CS%tracer_Reg, &
                            CS%sponge_CSp, CS%ALE_sponge_CSp, CS%OBC, Time_in)

  call create_group_pass(CS%pass_uv_T_S_h, CS%u, CS%v, G%Domain)

! Set up the parameters of the physical domain (i.e. the grid), G
!  call set_grid_metrics(G, PF)

! Set up the bottom depth, G%bathyT either analytically or from file
! This also sets G%max_depth based on the input parameter MAXIMUM_DEPTH,
! or, if absent, is diagnosed as G%max_depth = max( G%D(:,:) )
!  call MOM_initialize_topography(G%bathyT, G%max_depth, G, PF)

  node=>drifters

  !call mpp_domains_init

  ! create particles framework
!  call particles_framework_ini(filename,topname,parts)
!  grd=>parts%grd

  ! create drifters
!  call read_restart_particles(pafilename,parts) ! Modify from: drifters_input.f90: drifters_input_new

!  call particles_run(parts,time,grd%uo,grd%vo) ! Run the particles model

!  call particles_save_restart(parts)

END PROGRAM particle_driver
