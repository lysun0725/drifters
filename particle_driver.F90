PROGRAM particle_driver

  use MOM_file_parser, only : get_param, log_param, param_file_type
  use MOM_get_input, only : Get_MOM_Input, directories
  use MOM_domains, only : MOM_domains_init, MOM_infra_init, clone_MOM_domain, create_group_pass
  use MOM_grid, only : ocean_grid_type, MOM_grid_init
  use MOM_io, only : slasher, MOM_io_init, ASCII_FILE, READONLY_FILE, vardesc, var_desc
  use MOM_io, only : check_nml_error, file_exists, open_file, close_file, io_infra_init
  use MOM_error_handler, only : FATAL, WARNING, MOM_error
  use MOM_hor_index, only : hor_index_type, hor_index_init
  use MOM_grid, only : ocean_grid_type
  use MOM_grid_initialize, only : set_grid_metrics
  use MOM_fixed_initialization, only : MOM_initialize_fixed, MOM_initialize_topography
  use MOM_coord_initialization, only : MOM_initialize_coord
  use MOM_transcribe_grid,       only : copy_dyngrid_to_MOM_grid, copy_MOM_grid_to_dyngrid
  use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
  use MOM_verticalGrid, only : verticalGrid_type, verticalGridInit, get_thickness_units
  use MOM, only : MOM_control_struct
  use MOM_restart, only : restore_state, restart_init, register_restart_field, restart_end
  use MOM_diag_mediator, only : diag_mediator_infrastructure_init
  use time_manager_mod, only: time_type, set_time, set_date, JULIAN, NOLEAP, NO_CALENDAR, set_calendar_type
  use ensemble_manager_mod, only : get_ensemble_size, ensemble_manager_init, ensemble_pelist_setup
  use mpp_mod, only : set_current_pelist => mpp_set_current_pelist
  use particles_mod, only : particles_init
  use particles_mod, only : particles_run, particles_save_restart
!  use particles_io, only : read_restart_particles
  use particles_framework, only : particle, particles !, particles_gridded
  use particles_framework, only : really_debug
!  use particles_extra, only : particles_framework_ini, read_restart_particles
  implicit none

  type(ocean_grid_type), target :: Grid
  type(verticalGrid_type), target :: GV
  type(dyn_horgrid_type), pointer :: dG => NULL()
  type(MOM_control_struct), pointer, dimension(:) :: CSp !<pointer of ensemble list of MOM control structures
  type(MOM_control_struct), pointer :: CS=> NULL() !<pointer to an element in CSp
  type(hor_index_type)   :: HI ! A hor_index_type for array extents
  type(particles), target  :: drifters
  type(particles), pointer :: node=>NULL() !<pointer to an element in drifters list
  type(ocean_grid_type),  pointer :: ocn_grd
  type(time_type) :: time, time_start, time_start_segment, time_end, time_in
  real :: time_step
  character(len=128) :: history_file
  character(len=32)  :: topo_file
  character(len=32)  :: drifter_file
  character(len=48)  :: thickness_units
  character(len=8)  :: mod = 'MOM'
  type(param_file_type) :: PF
  type(directories) :: dirs
  type(vardesc) :: vd
  !<ens_info(1)=ensemble_size; ens_info(2)=tot_pes
  !<ens_info(3)=ocn_pes;ens_info(4)=atm_pes
  !<ens_info(5)=land_pes;ens_info(6)=ice_pes
  integer :: ens_size, nPEs_ocn, nPEs_atm, nPEs_lnd, nPES_ice, ens_info(6)
  integer, dimension(:), allocatable :: ocn_pelist, atm_pelist, lnd_pelist, ice_pelist
  integer :: calendar_type=-1
  character(len=16) :: month='jan',calendar='julian'
  logical :: concurrent = .false. !<use concurrent PE execution of model components
  integer :: date_init(6), date(6), years=0, months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: yr, mon, day, hr, min, sec
  integer :: is, ie, js, je, isd, ied, jsd, jed, IscB, IecB, JscB, JecB, IsdB, IedB, JsdB, JedB
  integer :: ierr, io_status, n, nz, unit

  ! Modify these later
  real :: dt=0.1
  integer :: axes(2) !< Diagnostic axes

  namelist /particle_driver_nml/ date_init, calendar, months, days, hours, minutes, seconds

  call MOM_infra_init(); call io_infra_init()
  call ensemble_manager_init(); ens_info = get_ensemble_size()
  ens_size=ens_info(1); nPEs_ocn=ens_info(2); nPEs_atm=nPEs_ocn; nPEs_lnd=nPEs_ocn; nPEs_ice=nPEs_ocn
  allocate(ocn_pelist(nPEs_ocn))
  allocate(atm_pelist(nPEs_ocn))
  allocate(lnd_pelist(nPEs_ocn))
  allocate(ice_pelist(nPEs_ocn))

  call ensemble_pelist_setup(concurrent,nPEs_atm,nPES_ocn,nPES_lnd,npes_ice,&
       atm_pelist,ocn_pelist,lnd_pelist,ice_pelist)
  call set_current_pelist(ocn_pelist)
  allocate(CSp(nPEs_ocn))

  call open_file(unit,'input.nml',form=ASCII_FILE,action=READONLY_FILE)
  read(unit,particle_driver_nml,iostat=io_status)
  call close_file(unit)

  ierr = check_nml_error(io_status,'particle_driver_nml')

  call Get_MOM_Input(PF,dirs)

  if (file_exists(trim(dirs%restart_input_dir)//'ocean_solo.res')) then
    call open_file(unit,trim(dirs%restart_input_dir)//'ocean_solo.res',form=ASCII_FILE,action=READONLY_FILE)
    read(unit,*) calendar_type; read(unit,*) date_init; read(unit,*) date
    call close_file(unit)
  endif

  if (calendar(1:6)=='julian') calendar_type=JULIAN
  if (calendar(1:11)=='no_calendar') calendar_type=NO_CALENDAR
  if (calendar(1:11)=='noleap') calendar_type=NOLEAP
  if (calendar_type<0) call MOM_error(FATAL,'paricles_driver: Invalid namelist value for calendar')

  call set_calendar_type(calendar_type)

  time_start = set_time(0,days=0)
  if (sum(date_init)>0) time_start=set_date(date_init(1),date_init(2), date_init(3), &
       date_init(4), date_init(5), date_init(6))

  if (sum(date)>0) then
    time_start_segment = time_start
    time=time_start_segment
  else
    time=time_start
  endif

  call MOM_domains_init(Grid%domain,PF)
  call MOM_io_init(PF)

  call diag_mediator_infrastructure_init()

  call hor_index_init(Grid%Domain, HI, PF, &
       local_indexing=.true.)

  call create_dyn_horgrid(dG,HI)
  call clone_MOM_domain(Grid%Domain, dG%Domain)


  do n=1,nPEs_ocn
    CS=>CSp(n)
    call verticalGridInit( PF, CS%GV )
!    call MOM_timing_init(CS)
!    call tracer_registry_init(PF,CS%tracer_Reg)
    call MOM_initialize_fixed(dG,CS%OBC,PF,.false.,dirs%output_directory)
    call MOM_initialize_coord(CS%GV, PF, .false., &
       dirs%output_directory, CS%tv, dG%max_depth)
  enddo

  is   = dG%isc   ; ie   = dG%iec  ; js   = dG%jsc  ; je   = dG%jec ; nz = CS%GV%ke
  isd  = dG%isd   ; ied  = dG%ied  ; jsd  = dG%jsd  ; jed  = dG%jed
  IsdB = dG%IsdB  ; IedB = dG%IedB ; JsdB = dG%JsdB ; JedB = dG%JedB

  call MOM_grid_init(Grid, PF, HI)
  call copy_dyngrid_to_MOM_grid(dG, Grid)

  thickness_units = get_thickness_units(GV)

  do n=1,nPEs_ocn
    CS=>CSp(n)
    CS%GV=>GV
    ! Allocate and initialize space for primary MOM variables.
    allocate(CS%u(IsdB:IedB,jsd:jed,nz))   ; CS%u(:,:,:) = 0.0
    allocate(CS%v(isd:ied,JsdB:JedB,nz))   ; CS%v(:,:,:) = 0.0
    allocate(CS%h(isd:ied,jsd:jed,nz))     ; CS%h(:,:,:) = GV%Angstrom

    call restart_init(PF,CS%restart_CSp)
    vd = var_desc("h",thickness_units,"Layer Thickness")
    call register_restart_field(CS%h, vd, .true., CS%restart_CSp)
    vd = var_desc("u","meter second-1","Zonal velocity",'u','L')
    call register_restart_field(CS%u, vd, .true., CS%restart_CSp)
    vd = var_desc("v","meter second-1","Meridional velocity",'v','L')
    call register_restart_field(CS%v, vd, .true., CS%restart_CSp)


    call restore_state(dirs%input_filename, dirs%restart_input_dir, Time, Grid,CS%restart_CSp)

  !   Shift from using the temporary dynamic grid type to using the final
  ! (potentially static) ocean-specific grid type.
  !   The next line would be needed if G%Domain had not already been init'd above:
  !     call clone_MOM_domain(dG%Domain, G%Domain)



    call create_group_pass(CS%pass_uv_T_S_h, CS%u, CS%v, Grid%Domain)

  enddo

  call destroy_dyn_horgrid(dG)
  node=>drifters
  ocn_grd=>Grid
  
  do n=1,nPEs_ocn
    CS=>CSp(n)
    call particles_init( node, ocn_grd, Time, dt, CS)
  enddo

  do n=1,nPEs_ocn
    CS=>CSp(n)
    call particles_run(node,time,CS%u(:,:,1),CS%v(:,:,1)) ! Run the particles model
  enddo

  call particles_save_restart(node)

END PROGRAM particle_driver
