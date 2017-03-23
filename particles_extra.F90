module particles_extra

  use constants_mod, only: radius, pi, omega, HLF 

  use mpp_domains_mod, only : domain2D
  use mpp_mod, only : mpp_npes, mpp_pe
  use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING

  use particles_framework, only : particle, particles, particles_gridded
  use particles_framework, only : nclasses
  use particles_framework, only : find_cell, find_cell_by_search
  use particles_framework, only : use_slow_find
  use particles_framework, only : pos_within_cell, add_new_part_to_list,print_part

  public particles_framework_ini,read_restart_particles

  !#ifdef _FILE_VERSION
  !  character(len=128) :: version = _FILE_VERSION
  !#else
  !  character(len=128) :: version = 'unknown'
  !#endif

  contains

  subroutine particles_framework_ini(infile,topfile,parts, maskmap) ! This can be used to replace "read_restart_particles" in particles_io.F90
    use netcdf
    use mpp_parameter_mod, only: SCALAR_PAIR, CGRID_NE, BGRID_NE, CORNER, AGRID
    use mpp_domains_mod, only: mpp_update_domains, mpp_define_domains
    use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain, mpp_get_global_domain
    use mpp_domains_mod, only: CYCLIC_GLOBAL_DOMAIN, FOLD_NORTH_EDGE
    use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST
    use mpp_domains_mod, only: mpp_define_io_domain

    use mpp_mod, only: mpp_clock_begin, mpp_clock_end, mpp_clock_id, input_nml_file
    use mpp_mod, only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP

    use fms_mod, only : open_namelist_file, close_file
    use fms_mod, only: clock_flag_default

    use diag_manager_mod, only: register_diag_field, register_static_field, send_data
    use diag_manager_mod, only: diag_axis_init

    use particles_framework, only: really_debug
    use particles_framework, only: unitTests

    ! Arguments
    type(particles), pointer:: parts
    character, intent(in) :: infile ! 19000101.ocean_hourly_1900_01_01_01.nc 
    character, intent(in) ::topfile
    logical, intent(in), optional :: maskmap(:,:)

    ! Local
    type(particles_gridded), pointer :: grd
    integer :: ncid, varid, ierr, iunit, npes
    integer :: gni,gnj,gnk
    integer :: id_class, axes3d(3)
    integer :: is, ie, js, je, i, j, debug=0, dr_lev=8 ! LUYU: modify this later
    real :: minl
    real, allocatable, dimension(:) :: ocean_lon, ocean_lat, ocean_lev, ocean_dx, ocean_dy
    real, allocatable, dimension(:,:) :: ocean_mask, cos_rot, sin_rot
    real, allocatable, dimension(:,:,:) :: u, v
    
    ! Local for now, arugument later. Need to be assigned values.
    integer :: layout(2), io_layout(2), axes(2)
    integer :: dom_x_flags, dom_y_flags
    real :: dt=1


    ! Namelist parameters (and defaults)
    integer :: halo=4 ! Width of halo region
    integer :: traj_sample_hrs=24 ! Period between sampling of position for trajectory storage
    integer :: traj_write_hrs=480 ! Period between writing sampled trajectories to disk
    integer :: verbose_hrs=24 ! Period between verbose messages
    real :: rho_parts=850. ! Density of particles
    real :: spring_coef=1.e-4  ! Spring contant for particle interactions - Alon
    real :: radial_damping_coef=1.e-4     ! Coef for relative particle motion damping (radial component) -Alon
    real :: tangental_damping_coef=2.e-5     ! Coef for relative particle motion damping (tangental component) -Alon
    real :: LoW_ratio=1.5 ! Initial ratio L/W for newly calved particles
    real :: party_bit_erosion_fraction=0. ! Fraction of erosion melt flux to divert to party bits
    real :: sicn_shift=0. ! Shift of sea-ice concentration in erosion flux modulation (0<sicn_shift<1)
    logical :: use_operator_splitting=.true. ! Use first order operator splitting for thermodynamics
    logical :: add_weight_to_ocean=.true. ! Add weight of particles + bits to ocean
    logical :: passive_mode=.false. ! Add weight of particles + bits to ocean
    logical :: time_average_weight=.false. ! Time average the weight on the ocean
    real :: speed_limit=0. ! CFL speed limit for a part
    real :: tau_calving=0. ! Time scale for smoothing out calving field (years)
    real :: tip_parameter=0. ! parameter to override particle rollilng critica ratio (use zero to get parameter directly from ice and seawater densities
    real :: grounding_fraction=0. ! Fraction of water column depth at which grounding occurs
    logical :: Runge_not_Verlet=.True.  !True=Runge Kutta, False=Verlet.  - Added by Alon 
    logical :: use_updated_rolling_scheme=.false. ! Use the corrected Rolling Scheme rather than the erronios one  
    logical :: use_new_predictive_corrective =.False.  !Flag to use Bob's predictive corrective particle scheme- Added by Alon 
    logical :: interactive_particles_on=.false.  !Turn on/off interactions between particles  - Added by Alon 
    logical :: critical_interaction_damping_on=.true.  !Sets the damping on relative particle velocity to critical value - Added by Alon 
    logical :: do_unit_tests=.false. ! Conduct some unit tests
    logical :: input_freq_distribution=.false. ! Flag to show if input distribution is freq or mass dist (=1 if input is a freq dist, =0 to use an input mass dist)
    logical :: read_old_restarts=.true. ! If true, read restarts prior to grid_of_lists and particle_num innovations
    real, dimension(nclasses) :: initial_mass=(/8.8e7, 4.1e8, 3.3e9, 1.8e10, 3.8e10, 7.5e10, 1.2e11, 2.2e11, 3.9e11, 7.4e11/) ! Mass thresholds between particle classes (kg)
    real, dimension(nclasses) :: distribution=(/0.24, 0.12, 0.15, 0.18, 0.12, 0.07, 0.03, 0.03, 0.03, 0.02/) ! Fraction of calving to apply to this class (non-dim) , 
    real, dimension(nclasses) :: mass_scaling=(/2000, 200, 50, 20, 10, 5, 2, 1, 1, 1/) ! Ratio between effective and real particle mass (non-dim)
    real, dimension(nclasses) :: initial_thickness=(/40., 67., 133., 175., 250., 250., 250., 250., 250., 250./) ! Total thickness of newly calved parts (m)
    namelist /particles_nml/ verbose, budget, halo, traj_sample_hrs, initial_mass, traj_write_hrs, &
             distribution, mass_scaling, initial_thickness, verbose_hrs, spring_coef, radial_damping_coef, tangental_damping_coef, &
             rho_parts, LoW_ratio, debug, really_debug, use_operator_splitting, party_bit_erosion_fraction, use_updated_rolling_scheme, &
             parallel_reprod, use_slow_find, sicn_shift, add_weight_to_ocean, passive_mode, ignore_ij_restart, use_new_predictive_corrective, tip_parameter, &
             time_average_weight, generate_test_particles, speed_limit, fix_restart_dates, Runge_not_Verlet, interactive_particles_on, critical_interaction_damping_on, &
             make_calving_reproduce,restart_input_dir, orig_read               ,do_unit_tests,grounding_fraction, input_freq_distribution, force_all_pes_traj, &
             read_old_restarts,tau_calving


    ! Get the stderr and stdlog unit numbers
    stderrunit=stderr()
    stdlogunit=stdlog()
    !write(stdlogunit,*) "particles_extra: "//trim(version)

    ! Read namelist parameters
    #ifdef INTERNAL_FILE_NML
      read (input_nml_file,nml=particles_nml,iostat=ierr)
    #else
      iunit = open_namelist_file()
      read (iunit,particles_nml,iostat=ierr)
      call close_file (iunit)
    #endif

    ! Allocate overall structure
    allocate(parts)
    allocate(parts%grd)
    grd=>parts%grd ! For convenience to avoid parts%grd%X
    allocate(grd%domain)

    ! Read the restart files     
    call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )
    call check( NF90_INQ_DIMID(ncid,'xh',varid) )
    call check( NF90_INQUIRE_DIMENSION(ncid,varid,len=gni) )
    call check( NF90_INQ_DIMID(ncid,'yh',varid) )
    call check( NF90_INQUIRE_DIMENSION(ncid,varid,len=gnj) )
    call check( NF90_INQ_DIMID(ncid,'zl',varid) )
    call check( NF90_INQUIRE_DIMENSION(ncid,varid,len=gnk) )

    allocate(ocean_lon(gni))
    allocate(ocean_lat(gnj))
    allocate(ocean_lev(gnk))
    allocate(ocean_dx(gni-1))
    allocate(ocean_dy(gnj-1))

    call check( NF90_INQ_VARID(ncid,'xh',varid) )
    call check( NF90_GET_VAR(ncid,varid,ocean_lon) )
    call check( NF90_INQ_VARID(ncid,'yh',varid) )
    call check( NF90_GET_VAR(ncid,varid,ocean_lat) )
    call check( NF90_INQ_VARID(ncid,'zl',varid) )
    call check( NF90_GET_VAR(ncid,varid,ocean_lev) )

    allocate(u(gni,gnj,gnk))
    allocate(v(gni,gnj,gnk))    

    call check( NF90_INQ_VARID(ncid,'u',varid) )
    call check( NF90_GET_VAR(ncid,varid,u) )
    call check( NF90_INQ_VARID(ncid,'v',varid) )
    call check( NF90_GET_VAR(ncid,varid,v) )
    call check( NF90_CLOSE(ncid) )

    ! Read the topog file
    allocate(ocean_mask(gni,gnj))
    call check( NF90_OPEN(topfile,NF90_NOWRITE,ncid) )
    call check( NF90_INQ_VARID(ncid,'wet',varid) )
    call check( NF90_GET_VAR(ncid,varid,ocean_mask) )
    call check( NF90_CLOSE(ncid) )

    ! Decompose domain
    !npe=mpp_npes()
    !call mpp_define_layout( (/1,gni,1,gnj/), npes, layout )
    !if (mpp_pe()==0) print *, 'LAYOUT', layout

    call mpp_define_domains( (/1,gni,1,gnj/), layout, grd%domain, & 
                              maskmap=maskmap, &
                              xflags=dom_x_flags, xhalo=halo, &
                              yflags=dom_y_flags, yhalo=halo, name='dimond' )
  
    call mpp_define_io_domain(grd%domain, io_layout)

    call mpp_get_compute_domain( grd%domain, grd%isc, grd%iec, grd%jsc, grd%jec )
    call mpp_get_data_domain( grd%domain, grd%isd, grd%ied, grd%jsd, grd%jed )
    call mpp_get_global_domain( grd%domain, grd%isg, grd%ieg, grd%jsg, grd%jeg )

    call mpp_get_neighbor_pe(grd%domain, NORTH, grd%pe_N)
    call mpp_get_neighbor_pe(grd%domain, SOUTH, grd%pe_S)
    call mpp_get_neighbor_pe(grd%domain, EAST, grd%pe_E)
    call mpp_get_neighbor_pe(grd%domain, WEST, grd%pe_W)


    folded_north_on_pe = ((dom_y_flags == FOLD_NORTH_EDGE) .and. (grd%jec == gnj))  

    allocate( grd%lon(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lon(:,:)=999.
    allocate( grd%lat(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lat(:,:)=999.
    allocate( grd%lonc(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lon(:,:)=999.
    allocate( grd%latc(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%lat(:,:)=999.
    allocate( grd%dx(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dx(:,:)=0.
    allocate( grd%dy(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%dy(:,:)=0.
    allocate( grd%msk(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%msk(:,:)=0.
    allocate( grd%cos(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%cos(:,:)=1.
    allocate( grd%sin(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%sin(:,:)=0.
    allocate( grd%ocean_depth(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%ocean_depth(:,:)=0.
    allocate( grd%uo(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%uo(:,:)=0.
    allocate( grd%vo(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%vo(:,:)=0.
    allocate( grd%parity_x(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_x(:,:)=1.
    allocate( grd%parity_y(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%parity_y(:,:)=1.
    allocate( grd%particle_num(grd%isd:grd%ied, grd%jsd:grd%jed) ); grd%particle_num(:,:)=1  

    ! Copy data declared on ocean model computational domain
    is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec
    do j = js, je
      grd%lon(is:ie,j)=ocean_lon(is:ie)
    end do
    do i = is, ie
      grd%lat(i,js:je)=ocean_lat(js:je)  
    end do

    grd%ocean_depth(is:ie,js:je)=ocean_lev(dr_lev) ! LUYU: Initialize this variable later and define 15m for now.

    ! Copy data declared on ocean model data domain ! LUYU: how to initialize this?
    is=grd%isc-1; ie=grd%iec+1; js=grd%jsc-1; je=grd%jec+1
    do j = js, je
      grd%dx(is:ie,j)=ocean_dx(is:ie)
    end do 
    do i = is, ie
      grd%dy(i,js:je)=ocean_dy(js:je)
    end do
    grd%msk(is:ie,js:je)=ocean_mask(is:ie,js:je)   
    grd%cos(is:ie,js:je)=cos_rot(is:ie,js:je)
    grd%sin(is:ie,js:je)=sin_rot(is:ie,js:je)

    ! Copy data declared on ocean velocity
    is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec
    grd%uo(is:ie,js:je)=u(is:ie,js:je,dr_lev)
    grd%vo(is:ie,js:je)=v(is:ie,js:je,dr_lev)   

    call mpp_update_domains(grd%lon, grd%domain, position=CORNER)
    call mpp_update_domains(grd%lat, grd%domain, position=CORNER)
    call mpp_update_domains(grd%dy, grd%dx, grd%domain, gridtype=CGRID_NE, flags=SCALAR_PAIR)

    call mpp_update_domains(grd%msk, grd%domain)
    call mpp_update_domains(grd%cos, grd%domain, position=CORNER)
    call mpp_update_domains(grd%sin, grd%domain, position=CORNER)
    call mpp_update_domains(grd%ocean_depth, grd%domain)
    call mpp_update_domains(grd%parity_x, grd%parity_y, grd%domain, gridtype=AGRID)

    ! Sanitize lon and lat at the SW edges
    do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
        if (grd%lon(i,j).gt.900.) grd%lon(i,j)=grd%lon(i,j+1)
        if (grd%lat(i,j).gt.900.) grd%lat(i,j)=2.*grd%lat(i,j+1)-grd%lat(i,j+2)
    enddo; enddo

    if (.not. present(maskmap)) then ! Using a maskmap causes tickles this sanity check
      do j=grd%jsd,grd%jed; do i=grd%isd,grd%ied
        if (grd%lon(i,j).gt.900.) write(stderrunit,*) 'bad lon: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1,grd%lon(i,j)
        if (grd%lat(i,j).gt.900.) write(stderrunit,*) 'bad lat: ',mpp_pe(),i-grd%isc+1,j-grd%jsc+1,grd%lat(i,j)
      enddo; enddo
    endif

   !The fix to reproduce across PE layout change, from AJA
    j=grd%jsc; do i=grd%isc+1,grd%ied
    minl=grd%lon(i-1,j)-180.
    if (abs(grd%lon(i,j)-(modulo(grd%lon(i,j)-minl,360.)+minl))>180.) &
       grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
    enddo
    j=grd%jsc; do i=grd%isc-1,grd%isd,-1
    minl=grd%lon(i+1,j)-180.
    if (abs(grd%lon(i,j)-(modulo(grd%lon(i,j)-minl,360.)+minl))>180.) &
       grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
    enddo
    do j=grd%jsc+1,grd%jed; do i=grd%isd,grd%ied
      minl=grd%lon(i,j-1)-180.
      if (abs(grd%lon(i,j)-(modulo(grd%lon(i,j)-minl,360.)+minl))>180.) &
         grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
    enddo; enddo
    do j=grd%jsc-1,grd%jsd,-1; do i=grd%isd,grd%ied
      minl=grd%lon(i,j+1)-180.
      if (abs(grd%lon(i,j)-(modulo(grd%lon(i,j)-minl,360.)+minl))>180.) &
         grd%lon(i,j)=modulo(grd%lon(i,j)-minl,360.)+minl
    enddo; enddo

    ! lonc, latc used for searches
    do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
      grd%lonc(i,j)=0.25*( (grd%lon(i,j)+grd%lon(i-1,j-1)) &
                          +(grd%lon(i-1,j)+grd%lon(i,j-1)) )
      grd%latc(i,j)=0.25*( (grd%lat(i,j)+grd%lat(i-1,j-1)) &
                          +(grd%lat(i-1,j)+grd%lat(i,j-1)) )
    enddo; enddo

    if (debug) then
      write(stderrunit,'(a,i3,a,4i4,a,4f8.2)') 'diamonds, particles_init: (',mpp_pe(),') [ij][se]c=', &
           grd%isc,grd%iec,grd%jsc,grd%jec, &
           ' [lon|lat][min|max]=', minval(grd%lon),maxval(grd%lon),minval(grd%lat),maxval(grd%lat)
    endif
    
!Added by Alon  - If a freq distribution is input, we have to convert the freq distribution to a mass flux distribution)
if (input_freq_distribution) then
     Total_mass=0.
     do j=1,nclasses
          Total_mass=Total_mass+(distribution(j)*initial_mass(j))
     enddo
     do j=1,nclasses
           distribution(j)=(distribution(j)*initial_mass(j))/Total_mass
     enddo
endif 


 ! Parameters
  parts%dt=dt
  parts%traj_sample_hrs=traj_sample_hrs
  parts%traj_write_hrs=traj_write_hrs
  parts%verbose_hrs=verbose_hrs
  parts%grd%halo=halo
  parts%rho_parts=rho_parts
  parts%spring_coef=spring_coef
  parts%radial_damping_coef=radial_damping_coef
  parts%tangental_damping_coef=tangental_damping_coef
  parts%LoW_ratio=LoW_ratio
  parts%use_operator_splitting=use_operator_splitting
  parts%party_bit_erosion_fraction=party_bit_erosion_fraction
  parts%sicn_shift=sicn_shift
  parts%passive_mode=passive_mode
  parts%time_average_weight=time_average_weight
  parts%speed_limit=speed_limit
  parts%tip_parameter=tip_parameter
  parts%Runge_not_Verlet=Runge_not_Verlet   !Alon
  parts%use_updated_rolling_scheme=use_updated_rolling_scheme  !Alon
  parts%critical_interaction_damping_on=critical_interaction_damping_on   !Alon
  parts%interactive_particles_on=interactive_particles_on   !Alon
  parts%use_new_predictive_corrective=use_new_predictive_corrective  !Alon
  parts%grounding_fraction=grounding_fraction
  parts%add_weight_to_ocean=add_weight_to_ocean
  parts%read_old_restarts=read_old_restarts
  allocate( parts%initial_mass(nclasses) ); parts%initial_mass(:)=initial_mass(:)
  allocate( parts%distribution(nclasses) ); parts%distribution(:)=distribution(:)
  allocate( parts%mass_scaling(nclasses) ); parts%mass_scaling(:)=mass_scaling(:)
  allocate( parts%initial_thickness(nclasses) ); parts%initial_thickness(:)=initial_thickness(:)
  allocate( parts%initial_width(nclasses) )
  allocate( parts%initial_length(nclasses) )
  parts%initial_width(:)=sqrt(initial_mass(:)/(LoW_ratio*rho_parts*initial_thickness(:)))
  parts%initial_length(:)=LoW_ratio*parts%initial_width(:)

  if (parts%read_old_restarts) call error_mesg('diamonds, ice_parts_framework_init', 'Setting "read_old_restarts=.true." can lead to non-reproducing checksums in restarts!', WARNING)

  ! Diagnostics
  id_class = diag_axis_init('mass_class', initial_mass, 'kg','Z', 'particle mass')
  axes3d(1:2)=axes
  axes3d(3)=id_class


  ! Static fields
  id_class=register_static_field('particles', 'lon', axes, &
               'longitude (corners)', 'degrees_E')
  if (id_class>0) lerr=send_data(id_class, grd%lon(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('particles', 'lat', axes, &
               'latitude (corners)', 'degrees_N')
  if (id_class>0) lerr=send_data(id_class, grd%lat(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('particles', 'mask', axes, &
               'wet point mask', 'none')
  if (id_class>0) lerr=send_data(id_class, grd%msk(grd%isc:grd%iec,grd%jsc:grd%jec))
  id_class=register_static_field('particles', 'ocean_depth', axes, &
               'ocean depth', 'm')
  if (id_class>0) lerr=send_data(id_class, grd%ocean_depth(grd%isc:grd%iec,grd%jsc:grd%jec))

  !if (debug) then
  !  call grd_chksum2(grd, grd%lon, 'init lon')
  !  call grd_chksum2(grd, grd%lat, 'init lat')
  !  call grd_chksum2(grd, grd%lonc, 'init lonc')
  !  call grd_chksum2(grd, grd%latc, 'init latc')
  !  call grd_chksum2(grd, grd%msk, 'init msk')
  !  call grd_chksum2(grd, grd%cos, 'init cos')
  !  call grd_chksum2(grd, grd%sin, 'init sin')
  !  call grd_chksum2(grd, grd%ocean_depth, 'init ocean_depth')
  !endif

  if (do_unit_tests) then
   if (unitTests(parts)) call error_mesg('diamonds, particles_init', 'Unit tests failed!', FATAL)
  endif

 !write(stderrunit,*) 'diamonds: done'
  call mpp_clock_end(parts%clock_ini)
  call mpp_clock_end(parts%clock)

   

  end subroutine particles_framework_ini

! ##############################################################################

subroutine read_restart_particles(infile,parts) ! This can be used to replace "read_restart_particles" in particles_io.F90

  use netcdf
  use particles_framework, only: really_debug

  ! Arguments
  type(particles), pointer :: parts
  character (len=15), intent(in) :: infile

  ! Local
  real, allocatable :: positions(:,:)
  integer, allocatable :: ids(:)
  integer :: ier, ncid, ndim, nfield, nparts, varid
  logical :: lres
  type(particles_gridded), pointer :: grd
  type(particle) :: localpart ! Not a pointer but an actual local variable
  integer :: stderrunit, iNg, jNg, i, j, k


  ! Get the stderr unit number
  stderrunit=stderr()
    
  ! For convenience
  grd=>parts%grd
  iNg=(grd%ieg-grd%isg+1) ! Total number of points globally in i direction, used with read_old_restarts=.true.
  jNg=(grd%jeg-grd%jsg+1) ! Total number of points globally in j direction, used with read_old_restarts=.true.
  
  call check( NF90_OPEN(infile,NF90_NOWRITE,ncid) )

  ! Read the dimension of this data file
  call check( NF90_INQ_DIMID(ncid,'nd',varid) )
  call check( NF90_INQUIRE_DIMENSION(ncid,varid,len=ndim) )
  call check( NF90_INQ_DIMID(ncid,'np',varid) )
  call check( NF90_INQUIRE_DIMENSION(ncid,varid,len=nparts) )
    
  ! Read the variable values of this data file
  allocate(positions(ndim,nparts))
  allocate(ids(nparts))

  call check( NF90_INQ_VARID(ncid, 'ids', varid) )
  call check( NF90_GET_VAR(ncid, varid, ids) )
  call check( NF90_INQ_VARID(ncid, 'positions', varid) )
  call check( NF90_GET_VAR(ncid, varid, positions) )    

  call check( NF90_CLOSE(ncid) )

  do k=1, nparts
    localpart%lon=positions(1,k)
    localpart%lat=positions(2,k)

    if (use_slow_find) then
       lres=find_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
    else
       lres=find_cell_by_search(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
    endif
    if (really_debug) then
      write(stderrunit,'(a,i8,a,2f9.4,a,i8)') 'diamonds, read_restart_particles: part ',k,' is at ',localpart%lon,localpart%lat,&
         & ' on PE ',mpp_pe()
      write(stderrunit,*) 'diamonds, read_restart_particles: lres = ',lres
    endif
    !if (really_debug) lres=is_point_in_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne, explain=.true.)
    lres=pos_within_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne, localpart%xi, localpart%yj)
    !call add_new_berg_to_list(bergs%first, localberg, quick=.true.)
    call add_new_part_to_list(parts%first, localpart)
    if (really_debug) call print_part(stderrunit, parts%first, 'read_restart_parts, add_new_part_to_list')
  end do


end subroutine read_restart_particles

!###########################################################################################################################

  subroutine check(status)
  !======================================================================
  ! Check the error status of the netcdf command
  !========================================================================
    use netcdf
    implicit none
    integer, intent (in) :: status
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check


end module
