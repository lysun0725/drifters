module particles_framework

use constants_mod, only: radius, pi, omega, HLF

use mpp_domains_mod, only: domain2D
use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_sync_self, mpp_pe, mpp_root_pe, mpp_chksum
use mpp_mod, only: COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4
use mpp_mod, only: COMM_TAG_5, COMM_TAG_6, COMM_TAG_7, COMM_TAG_8
use mpp_mod, only: COMM_TAG_9, COMM_TAG_10
use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING
use time_manager_mod, only: time_type, get_date, get_time, set_date, operator(-)

implicit none ; private

integer, parameter :: buffer_width=20 ! size of buffer dimension for comms
integer, parameter :: buffer_width_traj=20  ! LUYU: modify this later. use?
integer, parameter :: nclasses=10 ! Number of particles classes

logical :: folded_north_on_pe = .false.
logical :: verbose=.false. ! Be verbose to stderr
logical :: debug=.false. ! Turn on debugging
logical :: really_debug=.false. ! Turn on debugging
logical :: parallel_reprod=.true. ! Reproduce across different PE decompositions
logical :: use_slow_find=.true. ! Use really slow (but robust) find_cell for reading restarts
logical :: generate_test_particles=.false. ! Create particles in absence of a restart file
character(len=10) :: restart_input_dir = 'INPUT/'
integer, parameter :: delta_buf=25 ! Size by which to increment buffers
real, parameter :: pi_180=pi/180. ! Converts degrees to radians
logical :: fix_restart_dates=.true. ! After a restart, check that bergs were created before the current model date
logical :: do_unit_tests=.false. ! Conduct some unit tests
logical :: force_all_pes_traj=.false. ! Force all pes write trajectory files regardless of io_layout

!Public params !Niki: write a subroutine to expose these
public nclasses,buffer_width,buffer_width_traj
public verbose, really_debug, debug, restart_input_dir
public use_slow_find,generate_test_particles
public orig_read, force_all_pes_traj
public pi_180

!Public types
public particles_gridded, xyt, particle, particles, buffer 

!Public subs
public particles_framework_init
public send_parts_to_other_pes
public pack_part_into_buffer2, unpack_part_from_buffer2
public pack_traj_into_buffer2, unpack_traj_from_buffer2
public increase_buffer, increase_ibuffer, increase_ibuffer_traj, increase_buffer_traj
public add_new_part_to_list, count_out_of_order, check_for_duplicates
public insert_part_into_list, create_particle, delete_particle_from_list, destroy_particle
public print_fld,print_part,print_parts, record_posn, push_posn, append_posn, check_position
public move_trajectory, move_all_trajectories
public find_cell,find_cell_by_search,count_parts,is_point_in_cell,pos_within_cell
public sum_mass,sum_heat,bilin,yearday,parts_chksum
public checksum_gridded
public grd_chksum2,grd_chksum3
public fix_restart_dates
public offset_part_dates
public unitTests

type :: particles_gridded
  type(domain2D), pointer :: domain ! MPP domain
  integer :: halo ! Nominal halo width
  integer :: isc, iec, jsc, jec ! Indices of computational domain
  integer :: isd, ied, jsd, jed ! Indices of data domain
  integer :: isg, ieg, jsg, jeg ! Indices of global domain
  integer :: my_pe, pe_N, pe_S, pe_E, pe_W ! MPI PE identifiers
  real, dimension(:,:), pointer :: lon=>null() ! Longitude of cell corners
  real, dimension(:,:), pointer :: lat=>null() ! Latitude of cell corners
  real, dimension(:,:), pointer :: lonc=>null() ! Longitude of cell centers
  real, dimension(:,:), pointer :: latc=>null() ! Latitude of cell centers
  real, dimension(:,:), pointer :: dx=>null() ! Length of cell edge (m)
  real, dimension(:,:), pointer :: dy=>null() ! Length of cell edge (m)
  real, dimension(:,:), pointer :: msk=>null() ! Ocean-land mask (1=ocean)
  real, dimension(:,:), pointer :: cos=>null() ! Cosine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: sin=>null() ! Sine from rotation matrix to lat-lon coords
  real, dimension(:,:), pointer :: ocean_depth=>NULL() ! Depth of ocean (m)
  real, dimension(:,:), pointer :: uo=>null() ! Ocean zonal flow (m/s)
  real, dimension(:,:), pointer :: vo=>null() ! Ocean meridional flow (m/s)
  real, dimension(:,:), pointer :: tmp=>null() ! Temporary work space
  real, dimension(:,:), pointer :: tmpc=>null() ! Temporary work space
  real, dimension(:,:), pointer :: parity_x=>null() ! X component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  real, dimension(:,:), pointer :: parity_y=>null() ! Y component of vector point from i,j to i+1,j+1 (for detecting tri-polar fold)
  integer, dimension(:,:), pointer :: particle_num=>null() ! Counts particles created for naming purposes
end type particles_gridded

type :: xyt
  real :: lon, lat, day
  real :: uvel, vvel
  real :: axn, ayn, bxn, byn, uvel_old, vvel_old, lat_old, lon_old  !Explicit and implicit accelerations !Alon 
  real :: uo, vo
  integer :: year, particle_num
  type(xyt), pointer :: next=>null()
end type xyt

type :: particle
  type(particle), pointer :: prev=>null(), next=>null()
  ! State variables (specific to the particle, needed for restarts)
  real :: lon, lat, uvel, vvel, mass, thickness, width, length
  real :: axn, ayn, bxn, byn, uvel_old, vvel_old, lon_old, lat_old !Explicit and implicit accelerations !Alon 
  real :: start_lon, start_lat, start_day, start_mass, mass_scaling
  real :: mass_of_bits, heat_density
  integer :: start_year
  integer :: particle_num
  integer :: ine, jne ! nearest index in NE direction (for convenience)
  real :: xi, yj ! Non-dimensional coords within current cell (0..1)
  ! Environment variables (as seen by the particle)
  real :: uo, vo, ui, vi, ua, va, ssh_x, ssh_y, sst, cn, hi
  type(xyt), pointer :: trajectory=>null()
end type particle

type :: buffer
  integer :: size=0
  real, dimension(:,:), pointer :: data
end type buffer

type :: particles !; private!
  type(particles_gridded), pointer :: grd
  type(particle), pointer :: first=>null()
  type(xyt), pointer :: trajectories=>null()
  real :: dt           ! Time-step between particle calls (should make adaptive?)
  integer :: current_year
  real :: current_yearday ! 1.00-365.99
  integer :: traj_sample_hrs, traj_write_hrs
  integer :: verbose_hrs
  integer :: clock, clock_mom, clock_the, clock_int, clock_cal, clock_com, clock_ini, clock_ior, clock_iow, clock_dia ! ids for fms timers
  integer :: clock_trw, clock_trp
  real :: rho_parts ! Density of particles [kg/m^3]
  real :: spring_coef  ! Spring contant for particle interactions - Alon
  real :: radial_damping_coef     ! Coef for relative particle motion damping (radial component) -Alon
  real :: tangental_damping_coef     ! Coef for relative particle motion damping (tangental component) -Alon
  real :: LoW_ratio ! Initial ratio L/W for newly calved particles
  real :: party_bit_erosion_fraction ! Fraction of erosion melt flux to divert to party bits
  real :: sicn_shift ! Shift of sea-ice concentration in erosion flux modulation (0<sicn_shift<1)
  real, dimension(:), pointer :: initial_mass, distribution, mass_scaling
  real, dimension(:), pointer :: initial_thickness, initial_width, initial_length
  logical :: restarted=.false. ! Indicate whether we read state from a restart or not
  logical :: use_operator_splitting=.true. ! Use first order operator splitting for thermodynamics
  logical :: add_weight_to_ocean=.true. ! Add weight of parts to ocean
  logical :: passive_mode=.false. ! Add weight of particles + bits to ocean
  logical :: time_average_weight=.false. ! Time average the weight on the ocean
  logical :: use_updated_rolling_scheme=.false. ! True to use the aspect ratio based rolling scheme rather than incorrect version of WM scheme   (set tip_parameter=1000. for correct WM scheme)
  logical :: Runge_not_Verlet=.True.  !True=Runge Kuttai, False=Verlet.  - Added by Alon 
  logical :: use_new_predictive_corrective =.False.  !Flag to use Bob's predictive corrective particle scheme- Added by Alon 
  logical :: interactive_particles_on=.false.  !Turn on/off interactions between particles  - Added by Alon 
  logical :: critical_interaction_damping_on=.true.  !Sets the damping on relative particle velocity to critical value - Added by Alon 
  logical :: read_old_restarts=.true. ! If true, read restarts prior to grid_of_lists and particle_num innovation
  real :: speed_limit=0. ! CFL speed limit for a part [m/s]
  real :: tip_parameter=0. ! parameter to override particle rollilng critica ratio (use zero to get parameter directly from ice and seawater densities) 
  real :: grounding_fraction=0. ! Fraction of water column depth at which grounding occurs
  type(buffer), pointer :: obuffer_n=>null(), ibuffer_n=>null()
  type(buffer), pointer :: obuffer_s=>null(), ibuffer_s=>null()
  type(buffer), pointer :: obuffer_e=>null(), ibuffer_e=>null()
  type(buffer), pointer :: obuffer_w=>null(), ibuffer_w=>null()
  type(buffer), pointer :: obuffer_io=>null(), ibuffer_io=>null()
  ! Budgets
  real :: stored_start=0., stored_end=0.

  real :: party_mass_start=0., party_mass_end=0.
  real :: returned_mass_on_ocean=0.
  real :: net_melt=0., part_melt=0., party_src=0., party_melt=0.
  integer :: nparts_calved=0, nparts_melted=0, nparts_start=0, nparts_end=0
  integer :: nspeeding_tickets=0
  integer, dimension(:), pointer :: nparts_calved_by_class=>null()
end type particles

! Needs to be module global so can be public to particles_mod.
! Remove when backward compatibility no longer needed
logical :: orig_read=.false.

#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains


! ##############################################################################

subroutine particles_framework_init(parts, &
             gni, gnj, layout, io_layout, axes, dom_x_flags, dom_y_flags, &
             dt, Time, ice_lon, ice_lat, ice_wet, ice_dx, ice_dy, ice_area, &
             cos_rot, sin_rot, ocean_depth, maskmap, fractional_area)
! Arguments
type(particles), pointer :: parts
integer, intent(in) :: gni, gnj, layout(2), io_layout(2), axes(2)
integer, intent(in) :: dom_x_flags, dom_y_flags
real, intent(in) :: dt
type (time_type), intent(in) :: Time ! current time
real, dimension(:,:), intent(in) :: ice_lon, ice_lat, ice_wet
real, dimension(:,:), intent(in) :: ice_dx, ice_dy, ice_area
real, dimension(:,:), intent(in) :: cos_rot, sin_rot
real, dimension(:,:), intent(in),optional :: ocean_depth
logical, intent(in), optional :: maskmap(:,:)
logical, intent(in), optional :: fractional_area

end subroutine particles_framework_init

! ##############################################################################

subroutine offset_part_dates(parts,Time)
! Arguments
type(particles), pointer :: parts
type(time_type), intent(in) :: Time
! Local variables
type(particle), pointer :: this
integer :: iyr, imon, iday, ihr, imin, isec, yr_offset
real :: latest_start_year, part_start_year
real :: current_time_val

  call get_date(Time, iyr, imon, iday, ihr, imin, isec)
  latest_start_year=iyr-99999

  this=>parts%first
  if (associated(this)) then
   latest_start_year=float(this%start_year)+this%start_day/367.
  endif
  do while (associated(this))
    part_start_year=float(this%start_year)+this%start_day/367.
    if (part_start_year>latest_start_year) latest_start_year=part_start_year
    this=>this%next
  enddo
  call mpp_max(latest_start_year)

  current_time_val=float(iyr)+yearday(imon, iday, ihr, imin, isec)/367.
  if (latest_start_year<=current_time_val) return ! No conflicts!

  yr_offset=int(latest_start_year+1.)-iyr
  if (mpp_pe().eq.mpp_root_pe()) write(*,'(a,i8,a)') &
    'diamonds: Parts found with creation dates after model date! Adjusting part dates by ',yr_offset,' years'
  call parts_chksum(parts, 'before adjusting start dates')
  this=>parts%first
  do while (associated(this))
    this%start_year=this%start_year-yr_offset
    this=>this%next
  enddo
  call parts_chksum(parts, 'after adjusting start dates')

end subroutine offset_part_dates

! #############################################################################

subroutine send_parts_to_other_pes(parts)
! Arguments
type(particles), pointer :: parts
! Local variables
type(particle), pointer :: kick_the_bucket, this
integer :: nparts_to_send_e, nparts_to_send_w
integer :: nparts_to_send_n, nparts_to_send_s
integer :: nparts_rcvd_from_e, nparts_rcvd_from_w
integer :: nparts_rcvd_from_n, nparts_rcvd_from_s
type(particles_gridded), pointer :: grd
integer :: i, nparts_start, nparts_end
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>parts%grd

  if (debug) then
    nparts_start=count_parts(parts)
  endif

  ! Find number of parts that headed east/west
  nparts_to_send_e=0
  nparts_to_send_w=0
  if (associated(parts%first)) then
    this=>parts%first
    do while (associated(this))
      if (this%ine.gt.parts%grd%iec) then
        kick_the_bucket=>this
        this=>this%next
        nparts_to_send_e=nparts_to_send_e+1
        call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_e, nparts_to_send_e)
        call move_trajectory(parts, kick_the_bucket)
        call delete_particle_from_list(parts%first,kick_the_bucket)
      elseif (this%ine.lt.parts%grd%isc) then
        kick_the_bucket=>this
        this=>this%next
        nparts_to_send_w=nparts_to_send_w+1
        call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_w, nparts_to_send_w)
        call move_trajectory(parts, kick_the_bucket)
        call delete_particle_from_list(parts%first,kick_the_bucket)
      else
        this=>this%next
      endif
    enddo
  endif

  ! Send parts east
  if (grd%pe_E.ne.NULL_PE) then
    call mpp_send(nparts_to_send_e, plen=1, to_pe=grd%pe_E, tag=COMM_TAG_1)
    if (nparts_to_send_e.gt.0) then
      call mpp_send(parts%obuffer_e%data, nparts_to_send_e*buffer_width, grd%pe_E, tag=COMM_TAG_2)
    endif
  endif

  ! Send parts west
  if (grd%pe_W.ne.NULL_PE) then
    call mpp_send(nparts_to_send_w, plen=1, to_pe=grd%pe_W, tag=COMM_TAG_3)
    if (nparts_to_send_w.gt.0) then
      call mpp_send(parts%obuffer_w%data, nparts_to_send_w*buffer_width, grd%pe_W, tag=COMM_TAG_4)
    endif
  endif

  ! Receive parts from west
  if (grd%pe_W.ne.NULL_PE) then
    nparts_rcvd_from_w=-999
    call mpp_recv(nparts_rcvd_from_w, glen=1, from_pe=grd%pe_W, tag=COMM_TAG_1)
    if (nparts_rcvd_from_w.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_w,' from',grd%pe_W,' (W) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_w.gt.0) then
      call increase_ibuffer(parts%ibuffer_w, nparts_rcvd_from_w)
      call mpp_recv(parts%ibuffer_w%data, nparts_rcvd_from_w*buffer_width, grd%pe_W, tag=COMM_TAG_2)
      do i=1, nparts_rcvd_from_w
        call unpack_part_from_buffer2(parts%first, parts%ibuffer_w, i, grd)
      enddo
    endif
  else
    nparts_rcvd_from_w=0
  endif

  ! Receive parts from east
  if (grd%pe_E.ne.NULL_PE) then
    nparts_rcvd_from_e=-999
    call mpp_recv(nparts_rcvd_from_e, glen=1, from_pe=grd%pe_E, tag=COMM_TAG_3)
    if (nparts_rcvd_from_e.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_e,' from',grd%pe_E,' (E) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_e.gt.0) then
      call increase_ibuffer(parts%ibuffer_e, nparts_rcvd_from_e)
      call mpp_recv(parts%ibuffer_e%data, nparts_rcvd_from_e*buffer_width, grd%pe_E, tag=COMM_TAG_4)
      do i=1, nparts_rcvd_from_e
        call unpack_part_from_buffer2(parts%first, parts%ibuffer_e, i, grd)
      enddo
    endif
  else
    nparts_rcvd_from_e=0
  endif

  ! Find number of parts that headed north/south
  ! (note: this block should technically go ahead of the E/W recv block above
  !  to handle arbitrary orientation of PEs. But for simplicity, it is
  !  here to accomodate diagonal transfer of parts between PEs -AJA)
  nparts_to_send_n=0
  nparts_to_send_s=0
  if (associated(parts%first)) then
    this=>parts%first
    do while (associated(this))
      if (this%jne.gt.parts%grd%jec) then
        kick_the_bucket=>this
        this=>this%next
        nparts_to_send_n=nparts_to_send_n+1
        call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_n, nparts_to_send_n)
        call move_trajectory(parts, kick_the_bucket)
        call delete_particle_from_list(parts%first,kick_the_bucket)
      elseif (this%jne.lt.parts%grd%jsc) then
        kick_the_bucket=>this
        this=>this%next
        nparts_to_send_s=nparts_to_send_s+1
        call pack_part_into_buffer2(kick_the_bucket, parts%obuffer_s, nparts_to_send_s)
        call move_trajectory(parts, kick_the_bucket)
        call delete_particle_from_list(parts%first,kick_the_bucket)
      else
        this=>this%next
      endif
    enddo
  endif

  ! Send parts north
  if (grd%pe_N.ne.NULL_PE) then
    if(folded_north_on_pe) then
       call mpp_send(nparts_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_9)
    else 
       call mpp_send(nparts_to_send_n, plen=1, to_pe=grd%pe_N, tag=COMM_TAG_5)
    endif
    if (nparts_to_send_n.gt.0) then
       if(folded_north_on_pe) then
          call mpp_send(parts%obuffer_n%data, nparts_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
       else
          call mpp_send(parts%obuffer_n%data, nparts_to_send_n*buffer_width, grd%pe_N, tag=COMM_TAG_6)
       endif
    endif
  endif

  ! Send parts south
  if (grd%pe_S.ne.NULL_PE) then
    call mpp_send(nparts_to_send_s, plen=1, to_pe=grd%pe_S, tag=COMM_TAG_7)
    if (nparts_to_send_s.gt.0) then
      call mpp_send(parts%obuffer_s%data, nparts_to_send_s*buffer_width, grd%pe_S, tag=COMM_TAG_8)
    endif
  endif

  ! Receive parts from south
  if (grd%pe_S.ne.NULL_PE) then
    nparts_rcvd_from_s=-999
    call mpp_recv(nparts_rcvd_from_s, glen=1, from_pe=grd%pe_S, tag=COMM_TAG_5)
    if (nparts_rcvd_from_s.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_s,' from',grd%pe_S,' (S) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_s.gt.0) then
      call increase_ibuffer(parts%ibuffer_s, nparts_rcvd_from_s)
      call mpp_recv(parts%ibuffer_s%data, nparts_rcvd_from_s*buffer_width, grd%pe_S, tag=COMM_TAG_6)
      do i=1, nparts_rcvd_from_s
        call unpack_part_from_buffer2(parts%first, parts%ibuffer_s, i, grd) 
      enddo
    endif
  else
    nparts_rcvd_from_s=0
  endif

  ! Receive parts from north
  if (grd%pe_N.ne.NULL_PE) then
    nparts_rcvd_from_n=-999
    if(folded_north_on_pe) then
       call mpp_recv(nparts_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_9)
    else
       call mpp_recv(nparts_rcvd_from_n, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_7)
    endif
    if (nparts_rcvd_from_n.lt.0) then
      write(stderrunit,*) 'pe=',mpp_pe(),' received a bad number',nparts_rcvd_from_n,' from',grd%pe_N,' (N) !!!!!!!!!!!!!!!!!!!!!!'
    endif
    if (nparts_rcvd_from_n.gt.0) then
      call increase_ibuffer(parts%ibuffer_n, nparts_rcvd_from_n)
      if(folded_north_on_pe) then
         call mpp_recv(parts%ibuffer_n%data, nparts_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_10)
      else
         call mpp_recv(parts%ibuffer_n%data, nparts_rcvd_from_n*buffer_width, grd%pe_N, tag=COMM_TAG_8)
      endif
      do i=1, nparts_rcvd_from_n
        call unpack_part_from_buffer2(parts%first, parts%ibuffer_n, i, grd)
      enddo
    endif
  else
    nparts_rcvd_from_n=0
  endif

  if (debug) then
    nparts_end=count_parts(parts)
    i=nparts_rcvd_from_n+nparts_rcvd_from_s+nparts_rcvd_from_e+nparts_rcvd_from_w &
     -nparts_to_send_n-nparts_to_send_s-nparts_to_send_e-nparts_to_send_w
    if (nparts_end-(nparts_start+i).ne.0) then
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_end=',nparts_end,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_start=',nparts_start,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: delta=',i,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: error=',nparts_end-(nparts_start+i),' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_to_send_n=',nparts_to_send_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_to_send_s=',nparts_to_send_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_to_send_e=',nparts_to_send_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_to_send_w=',nparts_to_send_w,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_rcvd_from_n=',nparts_rcvd_from_n,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_rcvd_from_s=',nparts_rcvd_from_s,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_rcvd_from_e=',nparts_rcvd_from_e,' on PE',mpp_pe()
      write(stderrunit,'(a,i4,a,i4)') 'diamonds, send_parts_to_other_pes: nparts_rcvd_from_w=',nparts_rcvd_from_w,' on PE',mpp_pe()
      call error_mesg('diamonds, send_parts_to_other_pes:', 'We lost some parts!', FATAL)
    endif
  endif

  if (debug) then
    i=0
    this=>parts%first
    do while (associated(this))
      call check_position(grd, this, 'exchange (bot)')
      if (this%ine.lt.parts%grd%isc .or. &
          this%ine.gt.parts%grd%iec .or. &
          this%jne.lt.parts%grd%jsc .or. &
          this%jne.gt.parts%grd%jec) i=i+1
      this=>this%next
    enddo ! while
    call mpp_sum(i)
    if (i>0 .and. mpp_pe()==mpp_root_pe()) then
      write(stderrunit,'(a,i4)') 'diamonds, send_parts_to_other_pes: # of parts outside computational domain = ',i
      call error_mesg('diamonds, send_parts_to_other_pes:', 'there are parts still in halos!', FATAL)
    endif ! root_pe
  endif ! debug

  call mpp_sync_self()

end subroutine send_parts_to_other_pes

  subroutine pack_part_into_buffer2(part, buff, n)
  ! Arguments
  type(particle), pointer :: part
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  ! Local variables

    if (.not.associated(buff)) call increase_buffer(buff,delta_buf)
    if (n>buff%size) call increase_buffer(buff,delta_buf)

    buff%data(1,n)=part%lon
    buff%data(2,n)=part%lat
    buff%data(3,n)=part%uvel
    buff%data(4,n)=part%vvel
    buff%data(5,n)=part%xi
    buff%data(6,n)=part%yj
    buff%data(7,n)=part%start_lon
    buff%data(8,n)=part%start_lat
    buff%data(9,n)=float(part%start_year)
    buff%data(10,n)=part%start_day
    buff%data(11,n)=part%start_mass
    buff%data(12,n)=part%mass
    buff%data(13,n)=part%thickness
    buff%data(14,n)=part%width
    buff%data(15,n)=part%length
    buff%data(16,n)=part%mass_scaling
    buff%data(17,n)=part%mass_of_bits
    buff%data(18,n)=part%heat_density
    buff%data(19,n)=part%ine
    buff%data(20,n)=part%jne
    buff%data(21,n)=part%axn  !Alon
    buff%data(22,n)=part%ayn  !Alon
    buff%data(23,n)=part%bxn  !Alon
    buff%data(24,n)=part%byn  !Alon
    buff%data(25,n)=part%uvel_old  !Alon
    buff%data(26,n)=part%vvel_old  !Alon
    buff%data(27,n)=part%lon_old  !Alon
    buff%data(28,n)=part%lat_old  !Alon
    buff%data(29,n)=float(part%particle_num)

  end subroutine pack_part_into_buffer2

  subroutine increase_buffer(old,delta)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: delta
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size

    if (.not.associated(old)) then
      new_size=delta
    else
      new_size=old%size+delta
    endif
    allocate(new)
    allocate(new%data(buffer_width,new_size))
    new%size=new_size
    if (associated(old)) then
      new%data(:,1:old%size)=old%data(:,1:old%size)
      deallocate(old%data)
      deallocate(old)
    endif
    old=>new
   !write(stderr(),*) 'diamonds, increase_buffer',mpp_pe(),' increased to',new_size

  end subroutine increase_buffer

  subroutine unpack_part_from_buffer2(first, buff, n,grd, force_append)
  ! Arguments
  type(particle), pointer :: first
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  type(particles_gridded), pointer :: grd  
  logical, optional :: force_append
 ! Local variables
 !real :: lon, lat, uvel, vvel, xi, yj

 !real :: start_lon, start_lat, start_day, start_mass
 !integer :: ine, jne, start_year
  logical :: lres
  type(particle) :: localpart
  integer :: stderrunit
  logical :: force_app = .false.
  ! Get the stderr unit number
  stderrunit = stderr()

  if(present(force_append)) force_app = force_append
     
    localpart%lon=buff%data(1,n)
    localpart%lat=buff%data(2,n)
    localpart%uvel=buff%data(3,n)
    localpart%vvel=buff%data(4,n)
    localpart%xi=buff%data(5,n)
    localpart%yj=buff%data(6,n)
    localpart%start_lon=buff%data(7,n)
    localpart%start_lat=buff%data(8,n)
    localpart%start_year=nint(buff%data(9,n))
    localpart%start_day=buff%data(10,n)
    localpart%start_mass=buff%data(11,n)
    localpart%mass=buff%data(12,n)
    localpart%thickness=buff%data(13,n)
    localpart%width=buff%data(14,n)
    localpart%length=buff%data(15,n)
    localpart%mass_scaling=buff%data(16,n)
    localpart%mass_of_bits=buff%data(17,n)
    localpart%heat_density=buff%data(18,n)

    if(force_app) then !force append with origin ine,jne (for I/O)
       localpart%ine=buff%data(19,n) 
       localpart%jne=buff%data(20,n) 
       call add_new_part_to_list(first, localpart) 
    else
       
    localpart%axn=buff%data(21,n) !Alon
    localpart%ayn=buff%data(22,n) !Alon
    localpart%bxn=buff%data(23,n) !Alon
    localpart%byn=buff%data(24,n) !Alon
    localpart%uvel_old=buff%data(25,n) !Alon
    localpart%vvel_old=buff%data(26,n) !Alon
    localpart%lon_old=buff%data(27,n) !Alon
    localpart%lat_old=buff%data(28,n) !Alon
    localpart%particle_num=nint(buff%data(29,n))
   
     lres=find_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
    if (lres) then
      lres=pos_within_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne, localpart%xi, localpart%yj)
      call add_new_part_to_list(first, localpart)
    else
      lres=find_cell_wide(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
      if (lres) then
        lres=pos_within_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne, localpart%xi, localpart%yj)
        call add_new_part_to_list(first, localpart)
      else
        write(stderrunit,'("diamonds, unpack_part_from_buffer pe=(",i3,a,2i4,a,2f8.2)')&
             & mpp_pe(),') Failed to find i,j=',localpart%ine,localpart%jne,' for lon,lat=',localpart%lon,localpart%lat
        write(stderrunit,*) localpart%lon,localpart%lat
        write(stderrunit,*) localpart%uvel,localpart%vvel
        write(stderrunit,*) localpart%axn,localpart%ayn !Alon
        write(stderrunit,*) localpart%bxn,localpart%byn !Alon
        write(stderrunit,*) localpart%uvel_old,localpart%vvel_old !Alon
        write(stderrunit,*) localpart%lon_old,localpart%lat_old !Alon
        write(stderrunit,*) grd%isc,grd%iec,grd%jsc,grd%jec
        write(stderrunit,*) grd%isd,grd%ied,grd%jsd,grd%jed
        write(stderrunit,*) grd%lon(grd%isc-1,grd%jsc-1),grd%lon(grd%iec,grd%jsc)
        write(stderrunit,*) grd%lat(grd%isc-1,grd%jsc-1),grd%lat(grd%iec,grd%jec)
        write(stderrunit,*) grd%lon(grd%isd,grd%jsd),grd%lon(grd%ied,grd%jsd)
        write(stderrunit,*) grd%lat(grd%isd,grd%jsd),grd%lat(grd%ied,grd%jed)
        write(stderrunit,*) lres
        call error_mesg('diamonds, unpack_part_from_buffer', 'can not find a cell to place part in!', FATAL)
      endif
    endif
    endif

  end subroutine unpack_part_from_buffer2

  subroutine increase_ibuffer(old,delta)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: delta
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size, old_size

    if (.not.associated(old)) then
      new_size=delta+delta_buf
      old_size=0
    else
      old_size=old%size
      if (delta<old%size) then
        new_size=old%size+delta
      else
        new_size=delta+delta_buf
      endif
    endif

    if (old_size.ne.new_size) then
      allocate(new)
      allocate(new%data(buffer_width,new_size))
      new%size=new_size
      if (associated(old)) then
        new%data(:,1:old%size)=old%data(:,1:old%size)
        deallocate(old%data)
        deallocate(old)
      endif
      old=>new
     !write(stderr(),*) 'diamonds, increase_ibuffer',mpp_pe(),' increased to',new_size
    endif

  end subroutine increase_ibuffer

  subroutine increase_ibuffer_traj(old,delta)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: delta
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size, old_size

    if (.not.associated(old)) then
      new_size=delta+delta_buf
      old_size=0
    else
      old_size=old%size
      if (delta<old%size) then
        new_size=old%size+delta
      else
        new_size=delta+delta_buf
      endif
    endif

    if (old_size.ne.new_size) then
      allocate(new)
      allocate(new%data(buffer_width_traj,new_size))
      new%size=new_size
      if (associated(old)) then
        new%data(:,1:old%size)=old%data(:,1:old%size)
        deallocate(old%data)
        deallocate(old)
      endif
      old=>new
     !write(stderr(),*) 'diamonds, increase_ibuffer',mpp_pe(),' increased to',new_size
    endif

  end subroutine increase_ibuffer_traj

  subroutine increase_buffer_traj(old,delta)
  ! Arguments
  type(buffer), pointer :: old
  integer, intent(in) :: delta
  ! Local variables
  type(buffer), pointer :: new
  integer :: new_size

    if (.not.associated(old)) then
      new_size=delta
    else
      new_size=old%size+delta
    endif
    allocate(new)
    allocate(new%data(buffer_width_traj,new_size))
    new%size=new_size
    if (associated(old)) then
      new%data(:,1:old%size)=old%data(:,1:old%size)
      deallocate(old%data)
      deallocate(old)
    endif
    old=>new
   !write(stderr(),*) 'diamonds, increase_buffer',mpp_pe(),' increased to',new_size

  end subroutine increase_buffer_traj

  subroutine pack_traj_into_buffer2(traj, buff, n)
  ! Arguments
  type(xyt), pointer :: traj
  type(buffer), pointer :: buff
  integer, intent(in) :: n
  ! Local variables

    if (.not.associated(buff)) call increase_buffer_traj(buff,delta_buf)
    if (n>buff%size) call increase_buffer_traj(buff,delta_buf)

    buff%data(1,n)=traj%lon
    buff%data(2,n)=traj%lat
    buff%data(3,n)=float(traj%year)
    buff%data(4,n)=traj%day
    buff%data(5,n)=traj%uvel
    buff%data(6,n)=traj%vvel
    buff%data(13,n)=traj%uo
    buff%data(14,n)=traj%vo
    buff%data(24,n)=traj%axn !Alon
    buff%data(25,n)=traj%ayn !Alon
    buff%data(26,n)=traj%bxn !Alon
    buff%data(27,n)=traj%byn !Alon
    buff%data(28,n)=traj%uvel_old !Alon
    buff%data(29,n)=traj%vvel_old !Alon
    buff%data(30,n)=traj%lon_old !Alon
    buff%data(31,n)=traj%lat_old !Alon
    buff%data(32,n)=float(traj%particle_num)

  end subroutine pack_traj_into_buffer2

  subroutine unpack_traj_from_buffer2(first, buff, n)
  ! Arguments
  type(xyt), pointer :: first
  type(buffer), pointer :: buff
  integer, intent(in) :: n
 ! Local variables
  type(xyt) :: traj
  integer :: stderrunit
  ! Get the stderr unit number
  stderrunit = stderr()

    traj%lon=buff%data(1,n)
    traj%lat=buff%data(2,n)
    traj%year=nint(buff%data(3,n))
    traj%day=buff%data(4,n)
    traj%uvel=buff%data(5,n)
    traj%vvel=buff%data(6,n)
    traj%uo=buff%data(7,n)
    traj%vo=buff%data(8,n)
    traj%axn=buff%data(9,n) !Alon
    traj%ayn=buff%data(10,n) !Alon
    traj%bxn=buff%data(11,n) !Alon
    traj%byn=buff%data(12,n) !Alon
    traj%uvel_old=buff%data(13,n) !Alon
    traj%vvel_old=buff%data(14,n) !Alon
    traj%lon_old=buff%data(15,n) !Alon
    traj%lat_old=buff%data(16,n) !Alon
    traj%particle_num=nint(buff%data(17,n))

    call append_posn(first, traj)

  end subroutine unpack_traj_from_buffer2


! ##############################################################################

subroutine add_new_part_to_list(first, partvals, quick)
! Arguments
type(particle), pointer :: first
type(particle), intent(in) :: partvals
logical, intent(in), optional :: quick
! Local variables
type(particle), pointer :: new=>null()

  new=>null()
  call create_particle(new, partvals)

  if (present(quick)) then
    if(quick) call insert_part_into_list(first, new, quick=.true.)
  else
    call insert_part_into_list(first, new)
  endif

  !Clear new
  new=>null()

end subroutine add_new_part_to_list

! ##############################################################################

subroutine count_out_of_order(parts,label)
! Arguments
type(particles), pointer :: parts
character(len=*) :: label
! Local variables
type(particle), pointer :: this, next
integer :: i, icnt1, icnt2, icnt3

  icnt1=0; icnt3=0
  this=>parts%first
  next=>null()
  if (associated(this)) then
    if (associated(this%next)) next=>this%next
  endif
  do while (associated(next))
    if (.not. inorder(this,next)) icnt1=icnt1+1
    if (inorder(this,next).and.inorder(next,this)) icnt3=icnt3+1
    this=>next
    next=>next%next
  enddo
  call mpp_sum(icnt1)

  i=0;
  icnt2=0
  this=>parts%first
  do while (associated(this))
    i=1
    if (this%ine<parts%grd%isc .or. &
        this%ine>parts%grd%iec .or. &
        this%jne<parts%grd%jsc .or. &
        this%jne>parts%grd%jec) icnt2=icnt2+1
    this=>this%next
    if (i>1.and..not.associated(this%prev)) then
      call error_mesg('diamonds, count_out_of_order', 'Pointer %prev is unassociated. This should not happen!', FATAL)
    endif
  enddo
  call mpp_sum(icnt2)

  if ((debug.or.icnt1.ne.0).and.mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,3(x,a,i6),x,a)') 'diamonds, count_out_of_order:', &
      '# out of order=', icnt1,'# in halo=',icnt2,'# identicals=',icnt3,label
  endif

  call check_for_duplicates(parts,label)

end subroutine count_out_of_order

! ##############################################################################

subroutine check_for_duplicates(parts,label)
! Arguments
type(particles), pointer :: parts
character(len=*) :: label
! Local variables
type(particle), pointer :: this1, next1, this2, next2
integer :: icnt_id, icnt_same

  icnt_id=0
  icnt_same=0
  this1=>parts%first
  do while (associated(this1))
    this2=>this1%next
    do while (associated(this2))
      if (sameid(this1,this2)) icnt_id=icnt_id+1
      if (samepart(this1,this2)) icnt_same=icnt_same+1
      this2=>this2%next
    enddo
    this1=>this1%next
  enddo
  call mpp_sum(icnt_id)
  call mpp_sum(icnt_same)

  if ((debug.or.icnt_id>0.or.icnt_same>0).and.mpp_pe().eq.mpp_root_pe()) then
    write(*,'(a,2(x,a,i9),x,a)') 'diamonds, check_for_duplicates:', &
      '# with same id=', icnt_id,'# identical parts=',icnt_same,label
  endif

end subroutine check_for_duplicates

! ##############################################################################

subroutine insert_part_into_list(first, newpart, quick)
! Arguments
type(particle), pointer :: first, newpart
logical, intent(in), optional :: quick
! Local variables
type(particle), pointer :: this, prev
logical :: quickly = .false.

if(present(quick)) quickly = quick

  if (associated(first)) then
    if (.not. parallel_reprod .or. quickly) then
      newpart%next=>first
      first%prev=>newpart
      first=>newpart
    else
      if (inorder(newpart,first)) then
        ! Insert at front of list
        newpart%next=>first
        first%prev=>newpart
        first=>newpart
      else
        this=>first
        prev=>null()
        do while( associated(this) )
          if (inorder(newpart,this) ) then
            exit
          endif
          prev=>this
          this=>this%next
        enddo
        prev%next=>newpart
        newpart%prev=>prev
        if (associated(this)) this%prev=>newpart
        newpart%next=>this
      endif
    endif
  else
    ! list is empty so create it
    first=>newpart
  endif

end subroutine insert_part_into_list

! ##############################################################################

logical function inorder(part1, part2)  !MP Alon - Change to include particle_num
! Arguments
type(particle), pointer :: part1, part2
! Local variables
  if (part1%start_year<part2%start_year) then ! want newer first
    inorder=.true.
    return
  else if (part1%start_year>part2%start_year) then
    inorder=.false.
    return
  endif
  if (part1%start_day<part2%start_day) then ! want newer first
    inorder=.true.
    return
  else if (part1%start_day>part2%start_day) then
    inorder=.false.
    return
  endif
  if (part1%start_mass<part2%start_mass) then ! want lightest first
    inorder=.true.
    return
  else if (part1%start_mass>part2%start_mass) then
    inorder=.false.
    return
  endif
  if (part1%start_lon<part2%start_lon) then ! want eastward first
    inorder=.true.
    return
  else if (part1%start_lon>part2%start_lon) then
    inorder=.false.
    return
  endif
  if (part1%start_lat<part2%start_lat) then ! want southern first
    inorder=.true.
    return
  else if (part1%start_lat>part2%start_lat) then
    inorder=.false.
    return
  endif
  inorder=.true. ! passing the above tests mean the parts 1 and 2 are identical?
end function inorder

! ##############################################################################

  real function time_hash(part)!  Alon: Think about removing this.
  ! Arguments
  type(particle), pointer :: part
    time_hash=part%start_day+366.*float(part%start_year)
  end function time_hash

! ##############################################################################

  real function pos_hash(part)
  ! Arguments
  type(particle), pointer :: part
    pos_hash=part%start_lon+360.*(part%start_lat+90.)
  end function pos_hash

! ##############################################################################

logical function sameid(part1, part2) !  Alon: MP updat this.
! Arguments
type(particle), pointer :: part1, part2
! Local variables
  sameid=.false.
  if (part1%start_year.ne.part2%start_year) return
  if (part1%start_day.ne.part2%start_day) return
  if (part1%start_mass.ne.part2%start_mass) return
  if (part1%start_lon.ne.part2%start_lon) return
  if (part1%start_lat.ne.part2%start_lat) return
  sameid=.true. ! passing the above tests means that parts 1 and 2 have the same id
end function sameid

! ##############################################################################

logical function samepart(part1, part2)
! Arguments
type(particle), pointer :: part1, part2
! Local variables
  samepart=.false.
  if (.not. sameid(part1, part2)) return
  if (part1%lon.ne.part2%lon) return
  if (part1%lat.ne.part2%lat) return
  if (part1%mass.ne.part2%mass) return
  if (part1%uvel.ne.part2%uvel) return
  if (part1%vvel.ne.part2%vvel) return
  if (part1%thickness.ne.part2%thickness) return
  if (part1%width.ne.part2%width) return
  if (part1%length.ne.part2%length) return
  if (part1%axn.ne.part2%axn) return  !Alon
  if (part1%ayn.ne.part2%ayn) return  !Alon
  if (part1%bxn.ne.part2%bxn) return  !Alon
  if (part1%byn.ne.part2%byn) return  !Alon
  if (part1%uvel_old.ne.part2%uvel_old) return  !Alon
  if (part1%vvel_old.ne.part2%vvel_old) return  !Alon
  if (part1%lon_old.ne.part2%lon_old) return  !Alon
  if (part1%lat_old.ne.part2%lat_old) return  !Alon
  samepart=.true. ! passing the above tests mean that parts 1 and 2 are identical
end function samepart

! ##############################################################################

real function yearday(imon, iday, ihr, imin, isec)
! Arguments
integer, intent(in) :: imon, iday, ihr, imin, isec

  yearday=float(imon-1)*31.+float(iday-1)+(float(ihr)+(float(imin)+float(isec)/60.)/60.)/24.

end function yearday

! ##############################################################################

subroutine create_particle(part, partvals)
! Arguments
type(particle), pointer :: part
type(particle), intent(in) :: partvals
! Local variables
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  if (associated(part)) then
    write(stderrunit,*) 'diamonds, create_particle: part already associated!!!!',mpp_pe()
    call error_mesg('diamonds, create_particle', 'part already associated. This should not happen!', FATAL)
  endif
  allocate(part)
  part=partvals
  part%prev=>null()
  part%next=>null()

end subroutine create_particle


! ##############################################################################

subroutine delete_particle_from_list(first, part)
! Arguments
type(particle), pointer :: first, part
! Local variables

  ! Connect neighbors to each other
  if (associated(part%prev)) then
    part%prev%next=>part%next
  else
    first=>part%next
  endif
  if (associated(part%next)) part%next%prev=>part%prev

  ! Bye-bye part
  call destroy_particle(part)

end subroutine delete_particle_from_list

! ##############################################################################

subroutine destroy_particle(part)
! Arguments
type(particle), pointer :: part
! Local variables

  ! Bye-bye part
  deallocate(part)

end subroutine destroy_particle

! ##############################################################################

subroutine print_part(iochan, part, label)
! Arguments
integer, intent(in) :: iochan
type(particle), pointer :: part
character(len=*) :: label
! Local variables

  write(iochan,'("diamonds, print_part: ",a," pe=(",i3,") start lon,lat,yr,day,mass=",2f10.4,i5,i9,f7.2,es12.4)') &
    label, mpp_pe(), part%start_lon, part%start_lat, &
    part%start_year, part%particle_num, part%start_day, part%start_mass
  write(iochan,'("diamonds, print_part: ",a," pe=(",i3,a,2i5,7(a,2f10.4),a,2l2)') &
    label, mpp_pe(), ') i,j=',part%ine, part%jne, &
    ' xi,yj=', part%xi, part%yj, &
    ' lon,lat=', part%lon, part%lat, &
    ' u,v=', part%uvel, part%vvel, &
    ' axn,ayn=', part%axn, part%ayn, &
    ' bxn,byn=', part%bxn, part%byn, &
    ' uvel_old,vvel_old=', part%uvel_old, part%vvel_old, &
    ' lon_old,lat_old=', part%lon_old, part%lat_old, &
    ' p,n=', associated(part%prev), associated(part%next)
  write(iochan,'("diamonds, print_part: ",a," pe=(",i3,") ",6(a,2f10.4))') &
    label, mpp_pe(), 'uo,vo=', part%uo, part%vo, 'ua,va=', part%ua, part%va, 'ui,vi=', part%ui, part%vi
!Two lines above added by Alon
end subroutine print_part

! ##############################################################################

subroutine print_parts(iochan, parts, label)
! Arguments
integer, intent(in) :: iochan
type(particles), pointer :: parts
character(len=*) :: label
! Local variables
integer :: nparts, nnparts
type(particle), pointer :: this

  this=>parts%first
  do while(associated(this))
    call print_part(iochan, this, label)
    this=>this%next
  enddo
  nparts=count_parts(parts)
  nnparts=nparts
  call mpp_sum(nnparts)
  if (nparts.gt.0) write(iochan,'("diamonds, ",a," there are",i5," parts out of",i6," on PE ",i4)') label, nparts, nnparts, mpp_pe()

end subroutine print_parts

! ##############################################################################

integer function count_parts(parts)
! Arguments
type(particles), pointer :: parts
! Local variables
type(particle), pointer :: this

  count_parts=0
  this=>parts%first
  do while(associated(this))
    count_parts=count_parts+1
    this=>this%next
  enddo

end function count_parts

! ##############################################################################

subroutine record_posn(parts)
! Arguments
type(particles), pointer :: parts
! Local variables
type(xyt) :: posn
type(particle), pointer :: this

  this=>parts%first
  do while (associated(this))
    posn%lon=this%lon
    posn%lat=this%lat
    posn%year=parts%current_year
    posn%day=parts%current_yearday
    posn%uvel=this%uvel
    posn%vvel=this%vvel
    posn%uo=this%uo
    posn%vo=this%vo
    posn%axn=this%axn
    posn%ayn=this%ayn
    posn%bxn=this%bxn
    posn%byn=this%byn
    posn%uvel_old=this%uvel_old
    posn%vvel_old=this%vvel_old
    posn%lon_old=this%lon_old
    posn%lat_old=this%lat_old
    posn%particle_num=this%particle_num

    call push_posn(this%trajectory, posn)

    this=>this%next
  enddo

end subroutine record_posn

! ##############################################################################

subroutine push_posn(trajectory, posn_vals)
! Arguments
type(xyt), pointer :: trajectory
type(xyt) :: posn_vals
! Local variables
type(xyt), pointer :: new_posn

  allocate(new_posn)
  new_posn=posn_vals
  new_posn%next=>trajectory
  trajectory=>new_posn

end subroutine push_posn

subroutine append_posn(trajectory, posn_vals)
! This routine appends a new position leaf to the end of the given trajectory 
! Arguments
type(xyt), pointer :: trajectory
type(xyt) :: posn_vals
! Local variables
type(xyt), pointer :: new_posn,next,last

  allocate(new_posn)
  new_posn=posn_vals
  new_posn%next=>null()
  if(.NOT. associated(trajectory)) then
     trajectory=>new_posn
  else
     ! Find end of the trajectory and point it to the  new leaf
     next=>trajectory
     do while (associated(next))
        last=>next
        next=>next%next
     enddo
     last%next=>new_posn
  endif
end subroutine append_posn

! ##############################################################################

subroutine move_trajectory(parts, part)
! Arguments
type(particles), pointer :: parts
type(particle), pointer :: part
! Local variables
type(xyt), pointer :: next, last
type(xyt) :: vals

  ! If the trajectory is empty, ignore it
  if (.not.associated(part%trajectory)) return

  ! Push identifying info into first posn (note reverse order due to stack)
  vals%lon=part%start_lon
  vals%lat=part%start_lat
  vals%year=part%start_year
  vals%particle_num=part%particle_num
  vals%day=part%start_day
  call push_posn(part%trajectory, vals)
  vals%lon=0.
  vals%lat=99.
  vals%year=0
  vals%day=0.
  call push_posn(part%trajectory, vals)

  ! Find end of part trajectory and point it to start of existing trajectories
  next=>part%trajectory
  do while (associated(next))
    last=>next
    next=>next%next
  enddo
  last%next=>parts%trajectories

  parts%trajectories=>part%trajectory
  part%trajectory=>null()

end subroutine move_trajectory

! ##############################################################################

subroutine move_all_trajectories(parts, delete_parts)
! Arguments
type(particles),    pointer    :: parts
logical, optional, intent(in) :: delete_parts
! Local variables
type(particle), pointer :: this, next
logical :: delete_parts_after_moving_traj

  delete_parts_after_moving_traj = .false.
  if (present(delete_parts)) delete_parts_after_moving_traj = delete_parts
  this=>parts%first
  do while (associated(this))
    next=>this%next
    call move_trajectory(parts, this)
 !  if (delete_parts_after_moving_traj) call destroy_particle(this)
    this=>next
  enddo

end subroutine move_all_trajectories

! ##############################################################################

logical function find_cell_by_search(grd, x, y, i, j)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: x, y
integer, intent(inout) :: i, j
! Local variables
integer :: is,ie,js,je,di,dj,io,jo,icnt
real :: d0,d1,d2,d3,d4,d5,d6,d7,d8,dmin
logical :: explain=.false.

911 continue

  find_cell_by_search=.false.
  is=grd%isc; ie=grd%iec; js=grd%jsc; je=grd%jec

  ! Start at nearest corner
  d1=dcost(x,y,grd%lonc(is+1,js+1),grd%latc(is+1,js+1))
  d2=dcost(x,y,grd%lonc(ie-1,js+1),grd%latc(ie-1,js+1))
  d3=dcost(x,y,grd%lonc(ie-1,je-1),grd%latc(ie-1,je-1))
  d4=dcost(x,y,grd%lonc(is+1,je-1),grd%latc(is+1,je-1))
  dmin=min(d1,d2,d3,d4)
  if (d1==dmin) then; i=is+1; j=js+1
  elseif (d2==dmin) then; i=ie-1; j=js+1
  elseif (d3==dmin) then; i=ie-1; j=je-1
  elseif (d4==dmin) then; i=is+1; j=je-1
  else
    call error_mesg('diamonds, find_cell_by_search:', 'This should never EVER happen! (1)', FATAL)
  endif

  if (explain) then
    write(0,'(i3,a,2i4,f9.4)') mpp_pe(),'Initial corner i-is,j-js,cost=',i-is,j-js,dmin
    write(0,'(i3,a,3f9.4)') mpp_pe(),'cost ',d4,d3
    write(0,'(i3,a,3f9.4)') mpp_pe(),'cost ',d1,d2
  endif

  if (is_point_in_cell(grd, x, y, i, j)) then
    find_cell_by_search=.true.
    return
  endif
    
  do icnt=1, 1*(ie-is+je-js)
    io=i; jo=j

    d0=dcost(x,y,grd%lonc(io,jo),grd%latc(io,jo))
    d1=dcost(x,y,grd%lonc(io,jo+1),grd%latc(io,jo+1))
    d2=dcost(x,y,grd%lonc(io-1,jo+1),grd%latc(io-1,jo+1))
    d3=dcost(x,y,grd%lonc(io-1,jo),grd%latc(io-1,jo))
    d4=dcost(x,y,grd%lonc(io-1,jo-1),grd%latc(io-1,jo-1))
    d5=dcost(x,y,grd%lonc(io,jo-1),grd%latc(io,jo-1))
    d6=dcost(x,y,grd%lonc(io+1,jo-1),grd%latc(io+1,jo-1))
    d7=dcost(x,y,grd%lonc(io+1,jo),grd%latc(io+1,jo))
    d8=dcost(x,y,grd%lonc(io+1,jo+1),grd%latc(io+1,jo+1))

  ! dmin=min(d0,d1,d3,d5,d7)
    dmin=min(d0,d1,d2,d3,d4,d5,d6,d7,d8)
    if (d0==dmin) then; di=0; dj=0
    elseif (d2==dmin) then; di=-1; dj=1
    elseif (d4==dmin) then; di=-1; dj=-1
    elseif (d6==dmin) then; di=1; dj=-1
    elseif (d8==dmin) then; di=1; dj=1
    elseif (d1==dmin) then; di=0; dj=1
    elseif (d3==dmin) then; di=-1; dj=0
    elseif (d5==dmin) then; di=0; dj=-1
    elseif (d7==dmin) then; di=1; dj=0
    else
      call error_mesg('diamonds, find_cell_by_search:', 'This should never EVER happen!', FATAL)
    endif

    i=min(ie, max(is, io+di))
    j=min(je, max(js, jo+dj))

    if (explain) then
      write(0,'(i3,a,2i4,f9.5,a,2i4,a,2i4)') mpp_pe(),'Current position i,j,cost=',i,j,dmin,' di,dj=',di,dj,' old io,jo=',io,jo
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d2,d1,d8
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d3,d0,d7
      write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d4,d5,d6
    endif

    if (is_point_in_cell(grd, x, y, i, j)) then
      find_cell_by_search=.true.
      return
    endif
    
    if ((i==io.and.j==jo) &
        .and. .not.find_better_min(grd, x, y, 3, i, j) &
       ) then
      ! Stagnated
      find_cell_by_search=find_cell_loc(grd, x, y, is, ie, js, je, 1, i, j)
      if (.not. find_cell_by_search) find_cell_by_search=find_cell_loc(grd, x, y, is, ie, js, je, 3, i, j)
      i=min(ie, max(is, i))
      j=min(je, max(js, j))
      if (is_point_in_cell(grd, x, y, i, j)) then
        find_cell_by_search=.true.
      else
  !     find_cell_by_search=find_cell(grd, x, y, io, jo)
  !     if (find_cell_by_search) then
  !       if (explain) then
  !         call print_fld(grd, grd%lat, 'Lat')
  !         call print_fld(grd, grd%lon, 'Lat')
  !         do j=grd%jsd, grd%jed; do i=grd%isd, grd%ied
  !           grd%tmp(i,j) = dcost(x,y,grd%lonc(i,j),grd%latc(i,j))
  !         enddo; enddo
  !         call print_fld(grd, grd%tmp, 'Cost')
  !         stop 'Avoid recursing'
  !       endif
  !       write(0,'(i3,a,2i5,a,2i3,a,2f8.3)') mpp_pe(),'diamonds, find_cell_by_search: false negative io,jo=',io,jo,' di,dj=',di,dj,' targ=',x,y
  !       explain=.true.; goto 911
  !     endif
      endif
      return
    endif

  enddo

  find_cell_by_search=find_cell(grd, x, y, i, j)
  if (find_cell_by_search) then
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d2,d1,d8
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d3,d0,d7
    write(0,'(i3,a,3f9.5)') mpp_pe(),'cost ',d4,d5,d6
    write(0,'(i3,a,2f9.5)') mpp_pe(),'x,y ',x,y
    write(0,'(i3,a,4i5)') mpp_pe(),'io,jo ',io,jo,di,dj
    write(0,'(i3,a,2i5,a,2i3)') mpp_pe(),'diamonds, find_cell_by_search: false negative 2 i,j=',i-is,j-js,' di,dj=',di,dj
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'diamonds, find_cell_by_search: false negative 2 io,jo=',io,jo
    write(0,'(i3,a,2i5,a,2f8.3)') mpp_pe(),'diamonds, find_cell_by_search: false negative 2 i,j=',i,j,' targ=',x,y
    return
  endif
  find_cell_by_search=.false.

  contains

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

  real function dcost(x1, y1, x2, y2)
  ! Arguments
  real, intent(in) :: x1, x2, y1, y2
  ! Local variables
  real :: x1m

    x1m=modulo(x1-(x2-180.),360.)+(x2-180.)
  ! dcost=(x2-x1)**2+(y2-y1)**2
    dcost=(x2-x1m)**2+(y2-y1)**2
  end function dcost

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

  logical function find_better_min(grd, x, y, w, oi, oj)
  ! Arguments
  type(particles_gridded), intent(in) :: grd
  real, intent(in) :: x, y
  integer, intent(in) :: w
  integer, intent(inout) :: oi, oj
  ! Local variables
  integer :: i,j,xs,xe,ys,ye
  real :: dmin, dcst

  xs=max(grd%isc, oi-w)
  xe=min(grd%iec, oi+w)
  ys=max(grd%jsc, oj-w)
  ye=min(grd%jec, oj+w)

  find_better_min=.false.
  dmin=dcost(x,y,grd%lonc(oi,oj),grd%latc(oi,oj))
  do j=ys,ye; do i=xs,xe
      dcst=dcost(x,y,grd%lonc(i,j),grd%latc(i,j))
      if (dcst<dmin) then
        find_better_min=.true.
        dmin=dcst
        oi=i;oj=j
      endif
  enddo; enddo

  end function find_better_min

! # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

  logical function find_cell_loc(grd, x, y, is, ie, js, je, w, oi, oj)
  ! Arguments
  type(particles_gridded), intent(in) :: grd
  real, intent(in) :: x, y
  integer, intent(in) :: is, ie, js, je, w
  integer, intent(inout) :: oi, oj
  ! Local variables
  integer :: i,j,xs,xe,ys,ye

    xs=max(is, oi-w)
    xe=min(ie, oi+w)
    ys=max(js, oj-w)
    ye=min(je, oj+w)

    find_cell_loc=.false.
    do j=ys,ye; do i=xs,xe
        if (is_point_in_cell(grd, x, y, i, j)) then
          oi=i; oj=j; find_cell_loc=.true.
          return
        endif
    enddo; enddo

  end function find_cell_loc

end function find_cell_by_search

! ##############################################################################

logical function find_cell(grd, x, y, oi, oj)
! Arguments
type(particles_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(out) :: oi, oj
! Local variables
integer :: i,j

  find_cell=.false.; oi=-999; oj=-999

  do j=grd%jsc,grd%jec; do i=grd%isc,grd%iec
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell=.true.
        return
      endif
  enddo; enddo

end function find_cell

! ##############################################################################

logical function find_cell_wide(grd, x, y, oi, oj)
! Arguments
type(particles_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(out) :: oi, oj
! Local variables
integer :: i,j

  find_cell_wide=.false.; oi=-999; oj=-999

  do j=grd%jsd+1,grd%jed; do i=grd%isd+1,grd%ied
      if (is_point_in_cell(grd, x, y, i, j)) then
        oi=i; oj=j; find_cell_wide=.true.
        return
      endif
  enddo; enddo

end function find_cell_wide

! ##############################################################################

logical function is_point_in_cell(grd, x, y, i, j, explain)
! Arguments
type(particles_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
logical, intent(in), optional :: explain
! Local variables
real :: xlo, xhi, ylo, yhi
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  ! Safety check index bounds
  if (i-1.lt.grd%isd.or.i.gt.grd%ied.or.j-1.lt.grd%jsd.or.j.gt.grd%jed) then
    write(stderrunit,'(a,i3,(a,3i4))') &
                     'diamonds, is_point_in_cell: pe=(',mpp_pe(),') i,s,e=', &
                     i,grd%isd,grd%ied,' j,s,e=', j,grd%jsd,grd%jed
    call error_mesg('diamonds, is_point_in_cell', 'test is off the PE!', FATAL)
  endif

  is_point_in_cell=.false.

  ! Test crude bounds
  xlo=min( modulo(grd%lon(i-1,j-1)-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i  ,j-1)-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i-1,j  )-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i  ,j  )-(x-180.),360.)+(x-180.) )
  xhi=max( modulo(grd%lon(i-1,j-1)-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i  ,j-1)-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i-1,j  )-(x-180.),360.)+(x-180.), &
           modulo(grd%lon(i  ,j  )-(x-180.),360.)+(x-180.) )
  if (x.lt.xlo .or. x.gt.xhi) return
  ylo=min( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  yhi=max( grd%lat(i-1,j-1), grd%lat(i,j-1), grd%lat(i-1,j), grd%lat(i,j) )
  if (y.lt.ylo .or. y.gt.yhi) return
  
  if (grd%lat(i,j).gt.89.999) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, explain=explain) 
  elseif (grd%lat(i-1,j).gt.89.999) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i  ,j  ),grd%lat(i-1,j  ), &
                                        grd%lon(i-1,j-1),grd%lat(i-1,j  ), &
                                        x, y, explain=explain) 
  elseif (grd%lat(i-1,j-1).gt.89.999) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j  ),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, explain=explain) 
  elseif (grd%lat(i,j-1).gt.89.999) then
    is_point_in_cell=sum_sign_dot_prod5(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                        grd%lon(i-1,j-1),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j-1), &
                                        grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                        grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                        x, y, explain=explain) 
  else
  is_point_in_cell=sum_sign_dot_prod4(grd%lon(i-1,j-1),grd%lat(i-1,j-1), &
                                      grd%lon(i  ,j-1),grd%lat(i  ,j-1), &
                                      grd%lon(i  ,j  ),grd%lat(i  ,j  ), &
                                      grd%lon(i-1,j  ),grd%lat(i-1,j  ), &
                                      x, y, explain=explain) 
  endif

end function is_point_in_cell

! ##############################################################################

logical function sum_sign_dot_prod4(x0, y0, x1, y1, x2, y2, x3, y3, x, y, explain)
! Arguments
real, intent(in) :: x0, y0, x1, y1, x2, y2, x3, y3, x, y
logical, intent(in), optional :: explain
! Local variables
real :: p0,p1,p2,p3,xx
real :: l0,l1,l2,l3
real :: xx0,xx1,xx2,xx3
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  sum_sign_dot_prod4=.false.
  xx=modulo(x-(x0-180.),360.)+(x0-180.) ! Reference x to within 180 of x0
  xx0=modulo(x0-(x0-180.),360.)+(x0-180.) ! Reference x0 to within 180 of xx
  xx1=modulo(x1-(x0-180.),360.)+(x0-180.) ! Reference x1 to within 180 of xx
  xx2=modulo(x2-(x0-180.),360.)+(x0-180.) ! Reference x2 to within 180 of xx
  xx3=modulo(x3-(x0-180.),360.)+(x0-180.) ! Reference x3 to within 180 of xx

  l0=(xx-xx0)*(y1-y0)-(y-y0)*(xx1-xx0)
  l1=(xx-xx1)*(y2-y1)-(y-y1)*(xx2-xx1)
  l2=(xx-xx2)*(y3-y2)-(y-y2)*(xx3-xx2)
  l3=(xx-xx3)*(y0-y3)-(y-y3)*(xx0-xx3)

  p0=sign(1., l0); if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3); if (l3.eq.0.) p3=0.

  if ( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) .eq. abs((p0+p2)+(p1+p3)) ) then
    sum_sign_dot_prod4=.true.
  endif


  if (present(explain)) then
   if(explain) then
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: x=',mpp_pe(),':', &
                           x0,x1,x2,x3, x
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: X=',mpp_pe(),':', &
                           xx0,xx1,xx2,xx3, xx
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: y=',mpp_pe(),':', &
                           y0,y1,y2,y3, y
   write(stderrunit,'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod4: l=',mpp_pe(),':', &
                           l0,l1,l2,l3
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod4: p=',mpp_pe(),':', &
                           p0,p1,p2,p3, abs( (abs(p0)+abs(p2))+(abs(p1)+abs(p3)) - abs((p0+p2)+(p1+p3)) )
   endif
  endif

end function sum_sign_dot_prod4

! ##############################################################################

logical function sum_sign_dot_prod5(x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x, y, explain)
! Arguments
real, intent(in) :: x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x, y
logical, intent(in), optional :: explain
! Local variables
real :: p0,p1,p2,p3,p4,xx
real :: l0,l1,l2,l3,l4
real :: xx0,xx1,xx2,xx3,xx4
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  sum_sign_dot_prod5=.false.
  xx=modulo(x-(x0-180.),360.)+(x0-180.) ! Reference x to within 180 of x0
  xx0=modulo(x0-(x0-180.),360.)+(x0-180.) ! Reference x0 to within 180 of xx
  xx1=modulo(x1-(x0-180.),360.)+(x0-180.) ! Reference x1 to within 180 of xx
  xx2=modulo(x2-(x0-180.),360.)+(x0-180.) ! Reference x2 to within 180 of xx
  xx3=modulo(x3-(x0-180.),360.)+(x0-180.) ! Reference x3 to within 180 of xx
  xx4=modulo(x4-(x0-180.),360.)+(x0-180.) ! Reference x4 to within 180 of xx

  l0=(xx-xx0)*(y1-y0)-(y-y0)*(xx1-xx0)
  l1=(xx-xx1)*(y2-y1)-(y-y1)*(xx2-xx1)
  l2=(xx-xx2)*(y3-y2)-(y-y2)*(xx3-xx2)
  l3=(xx-xx3)*(y4-y3)-(y-y3)*(xx4-xx3)
  l4=(xx-xx4)*(y0-y4)-(y-y4)*(xx0-xx4)

  p0=sign(1., l0); if (l0.eq.0.) p0=0.
  p1=sign(1., l1); if (l1.eq.0.) p1=0.
  p2=sign(1., l2); if (l2.eq.0.) p2=0.
  p3=sign(1., l3); if (l3.eq.0.) p3=0.
  p4=sign(1., l4); if (l4.eq.0.) p4=0.

  if ( ((abs(p0)+abs(p2))+(abs(p1)+abs(p3)))+abs(p4) - abs(((p0+p2)+(p1+p3))+p4) .lt. 0.5 ) then
    sum_sign_dot_prod5=.true.
  endif

  if (present(explain)) then
   if(explain) then
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: x=',mpp_pe(),':', &
                           x0,x1,x2,x3,x4, x
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: X=',mpp_pe(),':', &
                           xx0,xx1,xx2,xx3,xx4, xx
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: y=',mpp_pe(),':', &
                           y0,y1,y2,y3,y4, y
   write(stderrunit,'(a,i3,a,1p10e12.4)') 'sum_sign_dot_prod5: l=',mpp_pe(),':', &
                           l0,l1,l2,l3,l4
   write(stderrunit,'(a,i3,a,10f12.4)') 'sum_sign_dot_prod5: p=',mpp_pe(),':', &
                           p0,p1,p2,p3,p4
   endif
  endif

end function sum_sign_dot_prod5

! ##############################################################################

logical function pos_within_cell(grd, x, y, i, j, xi, yj, explain)
! Arguments
type(particles_gridded), intent(in) :: grd
real, intent(in) :: x, y
integer, intent(in) :: i, j
real, intent(out) :: xi, yj
logical, intent(in), optional :: explain
! Local variables
real :: x1,y1,x2,y2,x3,y3,x4,y4,xx,yy,fac
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  pos_within_cell=.false.; xi=-999.; yj=-999.
  if (i-1<grd%isd) return
  if (j-1<grd%jsd) return
  if (i>grd%ied) return
  if (j>grd%jed) return

  x1=grd%lon(i-1,j-1)
  y1=grd%lat(i-1,j-1)
  x2=grd%lon(i  ,j-1)
  y2=grd%lat(i  ,j-1)
  x3=grd%lon(i  ,j  )
  y3=grd%lat(i  ,j  )
  x4=grd%lon(i-1,j  )
  y4=grd%lat(i-1,j  )

  if (present(explain)) then
    if(explain) then
    write(stderrunit,'(a,4f12.6)') 'pos_within_cell: x1..x4 ',x1,x2,x3,x4
    write(stderrunit,'(a,2f12.6)') 'pos_within_cell: x',x
    write(stderrunit,'(a,4f12.6)') 'pos_within_cell: y1..y4 ',y1,y2,y3,y4
    write(stderrunit,'(a,2f12.6)') 'pos_within_cell: y',y
    endif
  endif

  if (max(y1,y2,y3,y4)<89.999) then
    call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xi, yj, explain=explain)
  else
    if (debug) write(stderrunit,*) 'diamonds, pos_within_cell: working in tangential plane!'
    xx=(90.-y)*cos(x*pi_180)
    yy=(90.-y)*sin(x*pi_180)
    x1=(90.-y1)*cos(grd%lon(i-1,j-1)*pi_180)
    y1=(90.-y1)*sin(grd%lon(i-1,j-1)*pi_180)
    x2=(90.-y2)*cos(grd%lon(i  ,j-1)*pi_180)
    y2=(90.-y2)*sin(grd%lon(i  ,j-1)*pi_180)
    x3=(90.-y3)*cos(grd%lon(i  ,j  )*pi_180)
    y3=(90.-y3)*sin(grd%lon(i  ,j  )*pi_180)
    x4=(90.-y4)*cos(grd%lon(i-1,j  )*pi_180)
    y4=(90.-y4)*sin(grd%lon(i-1,j  )*pi_180)
    if (present(explain)) then
      if(explain) then
      write(stderrunit,'(a,4f12.6)') 'pos_within_cell: x1..x4 ',x1,x2,x3,x4
      write(stderrunit,'(a,2f12.6)') 'pos_within_cell: x',xx
      write(stderrunit,'(a,4f12.6)') 'pos_within_cell: y1..y4 ',y1,y2,y3,y4
      write(stderrunit,'(a,2f12.6)') 'pos_within_cell: y',yy
      endif
    endif
    call calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, xx, yy, xi, yj, explain=explain)
    if (is_point_in_cell(grd, x, y, i, j)) then
      if (abs(xi-0.5)>0.5.or.abs(yj-0.5)>0.5) then
        ! Scale internal coordinates to be consistent with is_point_in_cell()
        ! Note: this is intended to fix the inconsistency between the tangent plane
        ! and lat-lon calculations
        fac=2.*max( abs(xi-0.5), abs(yj-0.5) ); fac=max(1., fac)
        xi=0.5+(xi-0.5)/fac
        yj=0.5+(yj-0.5)/fac
        if (debug) call error_mesg('diamonds, pos_within_cell', 'in cell so scaling internal coordinates!', WARNING)
      endif
    else
      if (abs(xi-0.5)<=0.5.and.abs(yj-0.5)<=0.5) then
        if (debug) call error_mesg('diamonds, pos_within_cell', 'out of cell but coordinates <=0.5!', WARNING)
      endif
    endif
  endif

  if (present(explain)) then
	if(explain) write(stderrunit,'(a,2f12.6)') 'pos_within_cell: xi,yj=',xi,yj
  endif

 !if (.not. is_point_in_cell(grd, x, y, i, j) ) then
 !   write(stderrunit,'(a,i3,a,8f8.2,a)') 'diamonds, pos_within_cell: (',mpp_pe(),') ', &
 !                   x1, y1, x2, y2, x3, y3, x4, y4, ' NOT IN CELL!'
 !endif

  if (xi.ge.0. .and. xi.le.1. .and. yj.ge.0. .and. yj.le.1.) then
    pos_within_cell=is_point_in_cell(grd, x, y, i, j, explain=explain)
    if (.not. pos_within_cell .and. verbose) then
      if (debug) call error_mesg('diamonds, pos_within_cell', 'pos_within_cell is in cell BUT is_point_in_cell disagrees!', WARNING)
    endif
   !pos_within_cell=.true. ! commenting this out makes pos_within_cell agree with is_point_in_cell
  endif

  contains

  subroutine calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, xi, yj, explain)
  ! Arguments
  real,  intent(in) :: x1, x2, x3, x4, y1, y2, y3, y4, x, y
  real, intent(out) :: xi, yj
  logical, intent(in), optional :: explain
  ! Local variables
  real :: alpha, beta, gamma, delta, epsilon, kappa, a, b, c, d, dx, dy, yy1, yy2
  logical :: expl=.false.
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  expl=.false.
  if (present(explain)) then
     if(explain) expl=.true.
  endif
  alpha=x2-x1
  delta=y2-y1
  beta=x4-x1
  epsilon=y4-y1
  gamma=(x3-x1)-(alpha+beta)
  kappa=(y3-y1)-(delta+epsilon)
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs alpha,beta,gamma',alpha,beta,gamma,delta,epsilon,kappa
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs delta,epsilon,kappa',alpha,beta,gamma,delta,epsilon,kappa

  a=(kappa*beta-gamma*epsilon)
  dx=modulo(x-(x1-180.),360.)+(x1-180.)-x1
  dy=y-y1
  b=(delta*beta-alpha*epsilon)-(kappa*dx-gamma*dy)
  c=(alpha*dy-delta*dx)
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs dx,dy=',dx,dy
  if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs A,B,C=',a,b,c
  if (abs(a)>1.e-12) then
    d=0.25*(b**2)-a*c
    if (expl) write(stderrunit,'(a,1p6e12.4)') 'calc_xiyj: coeffs D=',d
    if (d.ge.0.) then
      if (expl) write(stderrunit,'(a,1p3e12.4)') 'Roots for b/2a, sqrt(d) = ',-0.5*b/a,sqrt(d)/a
      yy1=-(0.5*b+sqrt(d))/a
      yy2=-(0.5*b-sqrt(d))/a
      if (abs(yy1-0.5).lt.abs(yy2-0.5)) then; yj=yy1; else; yj=yy2; endif
      if (expl) write(stderrunit,'(a,1p3e12.4)') 'Roots for y = ',yy1,yy2,yj
    else
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: x1..x4 ',mpp_pe(),x1,x2,x3,x4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: y1..y4 ',mpp_pe(),y1,y2,y3,y4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderrunit,'(a,i3)') 'calc_xiyj: b<0 in quadratic root solver!!!!',mpp_pe()
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs a,b,c,d,dx,dy',mpp_pe(),a,b,c,d,dx,dy
      call error_mesg('diamonds, calc_xiyj', 'We have complex roots. The grid must be very distorted!', FATAL)
    endif
  else
    if (b.ne.0.) then
      yj=-c/b
    else
      yj=0.
    endif
  endif

  a=(alpha+gamma*yj)
  b=(delta+kappa*yj)
  if (a.ne.0.) then
    xi=(dx-beta*yj)/a
  elseif (b.ne.0.) then
    xi=(dy-epsilon*yj)/b
  else
    c=(epsilon*alpha-beta*delta)+(epsilon*gamma-beta*kappa)*yj
    if (c.ne.0.) then
      xi=(epsilon*dx-beta*dy)/c
    else
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: x1..x4 ',mpp_pe(),x1,x2,x3,x4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: x2..x4 - x1',mpp_pe(),x2-x1,x3-x1,x4-x1
      write(stderrunit,'(a,i3,4f8.2)') 'calc_xiyj: y1..y4 ',mpp_pe(),y1,y2,y3,y4
      write(stderrunit,'(a,i3,3f8.2)') 'calc_xiyj: y2..y4 - x1',mpp_pe(),y2-y1,y3-y1,y4-y1
      write(stderrunit,'(a,i3,1p6e12.4)') 'calc_xiyj: coeffs alpha..kappa',mpp_pe(),alpha,beta,gamma,delta,epsilon,kappa
      write(stderrunit,'(a,i3,1p2e12.4)') 'calc_xiyj: coeffs a,b',mpp_pe(),a,b
      call error_mesg('diamonds, calc_xiyj', 'Can not invert either linear equaton for xi! This should not happen!', FATAL)
    endif
  endif
  if (expl) write(stderrunit,'(a,2e12.4)') 'calc_xiyj: xi,yj=',xi,yj

  end subroutine calc_xiyj

end function pos_within_cell

! ##############################################################################

subroutine check_position(grd, part, label)
! Arguments
type(particles_gridded), pointer :: grd
type(particle), pointer :: part
character(len=*) :: label
! Local variables
real :: xi, yj
logical :: lret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  lret=pos_within_cell(grd, part%lon, part%lat, part%ine, part%jne, xi, yj)
  if (xi.ne.part%xi.or.yj.ne.part%yj) then
    write(stderrunit,'("diamonds: check_position (",i4,") b%x,x,-=",3(es12.4,x),a)') mpp_pe(),part%xi,xi,part%xi-xi,label
    write(stderrunit,'("diamonds: check_position (",i4,") b%y,y,-=",3(es12.4,x),a)') mpp_pe(),part%yj,yj,part%yj-yj,label
    call print_part(stderrunit, part, 'check_position')
    call error_mesg('diamonds, check_position','part has inconsistent xi,yj!',FATAL)
  endif

end subroutine check_position

! ##############################################################################

real function sum_mass(first,justbits,justparts)
! Arguments
type(particle), pointer :: first
logical, intent(in), optional :: justbits, justparts
! Local variables
type(particle), pointer :: this

  sum_mass=0.
  this=>first
  do while(associated(this))
    if (present(justparts)) then
      sum_mass=sum_mass+this%mass*this%mass_scaling
    elseif (present(justbits)) then
      sum_mass=sum_mass+this%mass_of_bits*this%mass_scaling
    else
      sum_mass=sum_mass+(this%mass+this%mass_of_bits)*this%mass_scaling
    endif
    this=>this%next
  enddo

end function sum_mass

! ##############################################################################

real function sum_heat(first,justbits,justparts)
! Arguments
type(particle), pointer :: first
logical, intent(in), optional :: justbits, justparts
! Local variables
type(particle), pointer :: this
real :: dm

  sum_heat=0.
  this=>first
  do while(associated(this))
    dm=0.
    if (present(justparts)) then
      dm=this%mass*this%mass_scaling
    elseif (present(justbits)) then
      dm=this%mass_of_bits*this%mass_scaling
    else
      dm=(this%mass+this%mass_of_bits)*this%mass_scaling
    endif
    sum_heat=sum_heat+dm*this%heat_density
    this=>this%next
  enddo

end function sum_heat


subroutine sanitize_field(arr,val)
! Arguments
real, dimension(:,:),intent(inout) :: arr
real, intent(in) :: val
! Local variables
integer :: i, j

  do j=lbound(arr,2), ubound(arr,2)
    do i=lbound(arr,1), ubound(arr,1)
      if (abs(arr(i,j)).ge.val) arr(i,j)=0.
    enddo
  enddo

end subroutine sanitize_field

! ##############################################################################




! ##############################################################################

subroutine checksum_gridded(grd, label)
! Arguments
type(particles_gridded), pointer :: grd
character(len=*) :: label
! Local variables

  if (mpp_pe().eq.mpp_root_pe()) write(*,'(2a)') 'diamonds: checksumming gridded data @ ',trim(label)

  ! external forcing
  call grd_chksum2(grd, grd%uo, 'uo')
  call grd_chksum2(grd, grd%vo, 'vo')
  
  ! static
  call grd_chksum2(grd, grd%lon, 'lon')
  call grd_chksum2(grd, grd%lat, 'lat')
  call grd_chksum2(grd, grd%lonc, 'lonc')
  call grd_chksum2(grd, grd%latc, 'latc')
  call grd_chksum2(grd, grd%dx, 'dx')
  call grd_chksum2(grd, grd%dy, 'dy')
  call grd_chksum2(grd, grd%msk, 'msk')
  call grd_chksum2(grd, grd%cos, 'cos')
  call grd_chksum2(grd, grd%sin, 'sin')
  call grd_chksum2(grd, grd%ocean_depth, 'depth')

end subroutine checksum_gridded

! ##############################################################################

subroutine grd_chksum3(grd, fld, txt)
! Arguments
type(particles_gridded), pointer :: grd
real, dimension(:,:,:), intent(in) :: fld
character(len=*), intent(in) :: txt
! Local variables
integer :: i, j, k, halo, icount, io, jo
real :: mean, rms, SD, minv, maxv
real, dimension(lbound(fld,1):ubound(fld,1), lbound(fld,2):ubound(fld,2), lbound(fld,3):ubound(fld,3)) :: tmp

  halo=grd%halo
  mean=0.
  rms=0.
  sd=0.
  icount=0
  i=lbound(fld,1)+halo
  j=lbound(fld,2)+halo
  k=lbound(fld,3)
  minv=fld(i,j,k)
  maxv=fld(i,j,k)
  tmp(:,:,:)=0.
  io=grd%isd-lbound(fld,1)
  jo=grd%jsd-lbound(fld,2)
  do k=lbound(fld,3), ubound(fld,3)
    do j=lbound(fld,2)+halo, ubound(fld,2)-halo
      do i=lbound(fld,1)+halo, ubound(fld,1)-halo
        icount=icount+1
        mean=mean+fld(i,j,k)
        rms=rms+fld(i,j,k)**2
        minv=min(minv,fld(i,j,k))
        maxv=max(maxv,fld(i,j,k))
        tmp(i,j,k)=fld(i,j,k)*float(i+io+2*(j+jo)+3*(k-1))
      enddo
    enddo
  enddo
  call mpp_sum(icount)
  call mpp_sum(mean)
  call mpp_sum(rms)
  call mpp_min(minv)
  call mpp_max(maxv)
  mean=mean/float(icount)
  rms=sqrt(rms/float(icount))
  do k=lbound(fld,3), ubound(fld,3)
    do j=lbound(fld,2)+halo, ubound(fld,2)-halo
      do i=lbound(fld,1)+halo, ubound(fld,1)-halo
        sd=sd+(fld(i,j,k)-mean)**2
      enddo
    enddo
  enddo
  call mpp_sum(sd)
  sd=sqrt(sd/float(icount))
  i=mpp_chksum( fld(lbound(fld,1)+halo:ubound(fld,1)-halo, &
                    lbound(fld,2)+halo:ubound(fld,2)-halo,:) )
  j=mpp_chksum( tmp(lbound(fld,1)+halo:ubound(fld,1)-halo, &
                    lbound(fld,2)+halo:ubound(fld,2)-halo,:) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum3: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  j=mpp_chksum( tmp(lbound(fld,1):ubound(fld,1), &
                    lbound(fld,2):ubound(fld,2),:) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum3* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9))') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd
#endif

end subroutine grd_chksum3

! ##############################################################################

subroutine grd_chksum2(grd, fld, txt)
! Arguments
type(particles_gridded), pointer :: grd
real, dimension(grd%isd:grd%ied,grd%jsd:grd%jed), intent(in) :: fld
character(len=*), intent(in) :: txt
! Local variables
integer :: i, j, icount
real :: mean, rms, SD, minv, maxv

  grd%tmp(:,:)=0.

  mean=0.
  rms=0.
  sd=0.
  icount=0
  minv=fld(grd%isc,grd%jsc)
  maxv=fld(grd%isc,grd%jsc)
  do j=grd%jsc, grd%jec
    do i=grd%isc, grd%iec
      icount=icount+1
      mean=mean+fld(i,j)
      rms=rms+fld(i,j)**2
      minv=min(minv,fld(i,j))
      maxv=max(maxv,fld(i,j))
      grd%tmp(i,j)=fld(i,j)*float(i+2*j)
    enddo
  enddo
  call mpp_sum(icount)
  call mpp_sum(mean)
  call mpp_sum(rms)
  call mpp_min(minv)
  call mpp_max(maxv)
  mean=mean/float(icount)
  rms=sqrt(rms/float(icount))
  do j=grd%jsc, grd%jec
    do i=grd%isc, grd%iec
      sd=sd+(fld(i,j)-mean)**2
    enddo
  enddo
  call mpp_sum(sd)
  sd=sqrt(sd/float(icount))
  i=mpp_chksum( fld(grd%isc:grd%iec,grd%jsc:grd%jec) )
  j=mpp_chksum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum2: ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i8)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#ifdef CHECKSUM_HALOS
  i=mpp_chksum( fld(grd%isd:grd%ied,grd%jsd:grd%jed) )
  j=mpp_chksum( grd%tmp(grd%isd:grd%ied,grd%jsd:grd%jed) )
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, grd_chksum2* ",a18,2(x,a,"=",i22),5(x,a,"=",es16.9),x,a,"=",i)') &
     txt, 'chksum', i, 'chksum2', j, 'min', minv, 'max', maxv, 'mean',  mean, 'rms', rms, 'sd', sd!, '#', icount
#endif

end subroutine grd_chksum2

! ##############################################################################

subroutine parts_chksum(parts, txt, ignore_halo_violation)
! Arguments
type(particles), pointer :: parts
character(len=*), intent(in) :: txt
logical, optional :: ignore_halo_violation
! Local variables
integer :: i, nparts, ichk1, ichk2, ichk3, ichk4, ichk5, ipart
real, allocatable :: fld(:,:), fld2(:,:)
integer, allocatable :: icnt(:,:)
type(particle), pointer :: this
type(particles_gridded), pointer :: grd
logical :: check_halo

! For convenience
  grd=>parts%grd

  nparts=count_parts(parts)
  call mpp_max(nparts)
  allocate( fld( nparts, 19 ) ) !Changed from 11 to 19 by Alon
  allocate( fld2( nparts, 19 ) ) !Changed from 11 to 19 by Alon
  allocate( icnt( grd%isd:grd%ied, grd%jsd:grd%jed ) )
  fld(:,:)=0.
  fld2(:,:)=0.
  icnt(:,:)=0
  grd%tmp(:,:)=0.

  this=>parts%first
  i=0; ichk5=0
  do while(associated(this))
    i=i+1
    ipart=part_chksum(this)
    fld(i,1) = this%lon
    fld(i,2) = this%lat
    fld(i,3) = this%uvel
    fld(i,4) = this%vvel
    fld(i,5) = this%mass
    fld(i,6) = this%thickness
    fld(i,7) = this%width
    fld(i,8) = this%length
    fld(i,9) = this%axn !added by Alon
    fld(i,10) = this%ayn !added by Alon
    fld(i,11) = this%bxn !added by Alon
    fld(i,12) = this%byn !added by Alon
    fld(i,13) = this%uvel_old !added by Alon
    fld(i,14) = this%vvel_old !added by Alon
    fld(i,15) = this%lon_old !added by Alon
    fld(i,16) = this%lat_old !added by Alon
    fld(i,17) = time_hash(this) !Changed from 9 to 17 by Alon
    fld(i,18) = pos_hash(this) !Changed from 10 to 18 by Alon
    fld(i,19) = float(ipart) !Changed from 11 to 19 by Alon
    icnt(this%ine,this%jne)=icnt(this%ine,this%jne)+1
    fld2(i,:) = fld(i,:)*float( icnt(this%ine,this%jne) ) !*float( i )
    grd%tmp(this%ine,this%jne)=grd%tmp(this%ine,this%jne)+time_hash(this)*pos_hash(this)+log(this%mass)
    ichk5=ichk5+ipart
    this=>this%next
  enddo

  ichk1=mpp_chksum( fld )
  ichk2=mpp_chksum( fld2 )
  ichk3=mpp_chksum( grd%tmp )
  ichk4=mpp_chksum( grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec) )
  call mpp_sum( ichk5 )
  nparts=count_parts(parts)

  if (nparts.ne.sum(icnt(:,:))) then
    write(*,'("diamonds, parts_chksum: ",2(a,i8))') &
      '# parts =', nparts, ' sum(icnt) =',sum(icnt(:,:))
    call error_mesg('diamonds, parts_chksum:', 'mismatch in part count!', FATAL)
  endif

  check_halo=.true.
  if (present(ignore_halo_violation)) then
    if (ignore_halo_violation) check_halo=.false.
  endif
  if (check_halo.and.nparts.ne.sum(icnt(grd%isc:grd%iec, grd%jsc:grd%jec))) then
    write(*,'("diamonds, parts_chksum: ",2(a,i8))') &
      '# parts =', nparts, ' sum(icnt(comp_dom)) =',sum(icnt(:,:))
    call error_mesg('diamonds, parts_chksum:', 'mismatch in part count on computational domain!', FATAL)
  endif

  call mpp_sum(nparts)
  if (mpp_pe().eq.mpp_root_pe()) &
    write(*,'("diamonds, parts_chksum: ",a18,6(x,a,"=",i22))') &
      txt, 'chksum', ichk1, 'chksum2', ichk2, 'chksum3', ichk3, 'chksum4', ichk4, 'chksum5', ichk5, '#', nparts

  grd%tmp(:,:)=real(icnt(:,:))
  call grd_chksum2(grd,grd%tmp,'# of parts/cell')

  deallocate( fld )
  deallocate( fld2 )
  deallocate( icnt )

  if (debug) call count_out_of_order(parts,txt)

end subroutine parts_chksum

! ##############################################################################

integer function part_chksum(part )
! Arguments
type(particle), pointer :: part
! Local variables
real :: rtmp(36) !Changed from 28 to 34 by Alon
integer :: itmp(36+4), i8=0, ichk1, ichk2, ichk3 !Changed from 28 to 34 by Alon
integer :: i

  rtmp(:)=0.
  rtmp(1)=part%lon
  rtmp(2)=part%lat
  rtmp(3)=part%uvel
  rtmp(4)=part%vvel
  rtmp(5)=part%mass
  rtmp(6)=part%thickness
  rtmp(7)=part%width
  rtmp(8)=part%length
  rtmp(9)=part%start_lon
  rtmp(10)=part%start_lat
  rtmp(11)=part%start_day
  rtmp(12)=part%start_mass
  rtmp(13)=part%mass_scaling
  rtmp(14)=part%mass_of_bits
  rtmp(15)=part%heat_density
  rtmp(16)=part%xi
  rtmp(17)=part%yj
  rtmp(19)=part%uo
  rtmp(20)=part%vo
  rtmp(21)=part%ui
  rtmp(22)=part%vi
  rtmp(23)=part%ua
  rtmp(24)=part%va
  rtmp(25)=part%ssh_x
  rtmp(26)=part%ssh_y
  rtmp(27)=part%cn
  rtmp(28)=part%hi
  rtmp(29)=part%axn !Added by Alon
  rtmp(30)=part%ayn !Added by Alon
  rtmp(31)=part%bxn !Added by Alon
  rtmp(32)=part%byn !Added by Alon
  rtmp(33)=part%uvel_old !Added by Alon
  rtmp(34)=part%vvel_old !Added by Alon
  rtmp(35)=part%lat_old !Added by Alon
  rtmp(36)=part%lon_old !Added by Alon

  itmp(1:36)=transfer(rtmp,i8) !Changed from 28 to 36 by Alon
  itmp(37)=part%start_year !Changed from 29 to 37 by Alon
  itmp(38)=part%ine !Changed from 30 to 38 by Alon
  itmp(39)=part%jne !Changed from 31 to 39 by Alon
  itmp(40)=part%particle_num !added  by Alon

  ichk1=0; ichk2=0; ichk3=0
  do i=1,37+3 !Changd from 28 to 37 by Alon
   ichk1=ichk1+itmp(i)
   ichk2=ichk2+itmp(i)*i
   ichk3=ichk3+itmp(i)*i*i
  enddo
  part_chksum=ichk1+ichk2+ichk3

end function part_chksum

! ##############################################################################

real function bilin(grd, fld, i, j, xi, yj)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed), xi, yj
integer, intent(in) :: i, j
! Local variables


    bilin=(fld(i,j  )*xi+fld(i-1,j  )*(1.-xi))*yj &
         +(fld(i,j-1)*xi+fld(i-1,j-1)*(1.-xi))*(1.-yj)

end function bilin

! ##############################################################################

subroutine print_fld(grd, fld, label)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(in) :: fld(grd%isd:grd%ied,grd%jsd:grd%jed)
character(len=*) :: label
! Local variables
integer :: i, j
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  write(stderrunit,'("pe=",i3,x,a8,32i10)') mpp_pe(),label,(i,i=grd%isd,grd%ied)
  do j=grd%jed,grd%jsd,-1
    write(stderrunit,'("pe=",i3,x,i8,32es10.2)') mpp_pe(),j,(fld(i,j),i=grd%isd,grd%ied)
  enddo

end subroutine print_fld

! ##############################################################################

logical function unitTests(parts)
  type(particles), pointer :: parts
  type(particles_gridded), pointer :: grd
  ! Local variables
  integer :: stderrunit,i,j

  ! This function returns True is a unit test fails
  unitTests=.false.
  ! For convenience
  grd=>parts%grd
  stderrunit=stderr()
  
  i=grd%isc; j=grd%jsc
  call localTest( bilin(grd, grd%lon, i, j, 0., 1.), grd%lon(i-1,j) )
  call localTest( bilin(grd, grd%lon, i, j, 1., 1.), grd%lon(i,j) )
  call localTest( bilin(grd, grd%lat, i, j, 1., 0.), grd%lat(i,j-1) )
  call localTest( bilin(grd, grd%lat, i, j, 1., 1.), grd%lat(i,j) )

  contains
  subroutine localTest(answer, rightAnswer)
  real, intent(in) :: answer, rightAnswer
  if (answer==rightAnswer) return
  unitTests=.true.
  write(stderrunit,*) 'a=',answer,'b=',rightAnswer
  end subroutine localTest
end function unitTests

! ##############################################################################

end module
