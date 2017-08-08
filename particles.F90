module particles_mod

use constants_mod, only: radius, pi, omega, HLF
use MOM_grid, only : ocean_grid_type

use fms_mod, only: field_exist, get_global_att_value
use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING
use fms_mod, only: write_version_number, read_data, write_data, file_exist
use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_sync_self, mpp_chksum
use mpp_mod, only: mpp_clock_begin, mpp_clock_end, mpp_clock_id
use mpp_mod, only: CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP

use mpp_mod, only: mpp_gather
use fms_mod, only: clock_flag_default
use fms_io_mod, only: get_instance_filename
use mpp_domains_mod, only: domain2D, mpp_update_domains, mpp_define_domains
use mpp_parameter_mod, only: CGRID_NE, BGRID_NE, CORNER, AGRID
use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_data_domain
use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST
use time_manager_mod, only: time_type, get_date, get_time, set_date, operator(-)
use diag_manager_mod, only: register_diag_field, register_static_field, send_data
use diag_manager_mod, only: diag_axis_init

use particles_framework, only: particles_framework_init
use particles_framework, only: particles_gridded, xyt, particle, particles, buffer
use particles_framework, only: verbose, really_debug,debug
use particles_framework, only: find_cell,find_cell_by_search,count_parts,is_point_in_cell,pos_within_cell
use particles_framework, only: nclasses
use particles_framework, only: bilin,yearday,count_parts,parts_chksum
use particles_framework, only: checksum_gridded,add_new_part_to_list
use particles_framework, only: send_parts_to_other_pes,move_trajectory,move_all_trajectories
use particles_framework, only: record_posn,check_position,print_part,print_parts,print_fld
use particles_framework, only: add_new_part_to_list,delete_particle_from_list,destroy_particle
use particles_framework, only: grd_chksum2,grd_chksum3
use particles_framework, only: offset_part_dates

use particles_io,        only: particles_io_init,write_restart,write_trajectory
!use particles_io,        only: read_restart_parts,read_restart_parts_orig,read_restart_calving

implicit none ; private

public particles_init !, particles_end, particles_run, particles_stock_pe, particles
public particles_end, particles_run, particles
public particles_save_restart


#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains

! ##############################################################################
subroutine particles_init(parts, Grid, Time, dt, axes)

 use particles_io, only: read_restart_parts
 
 type(particles), pointer :: parts
 type(ocean_grid_type), pointer :: Grid !< Grid type from parent model
 type(time_type), intent(in) :: Time !< Time type from parent model
 real, intent(in)            :: dt !< particle timestep in seconds
 integer, dimension(2), intent(in) :: axes !< diagnostic axis ids


!subroutine particles_init(parts, &
!             gni, gnj, layout, io_layout, axes, dom_x_flags, dom_y_flags, &
!             dt, Time, lon, lat, wet, dx, dy,area, &
!             cos_rot, sin_rot, ocean_depth)


! Arguments
! type(particles), pointer :: parts
! integer, intent(in) :: gni, gnj, layout(2), io_layout(2), axes(2)
! integer, intent(in) :: dom_x_flags, dom_y_flags
! real, intent(in) :: dt
! type (time_type), intent(in) :: Time ! current time
! real, dimension(:,:), intent(in) :: lon, lat, wet
! real, dimension(:,:), intent(in) :: dx, dy, area
! real, dimension(:,:), intent(in) :: cos_rot, sin_rot
! real, dimension(:,:), intent(in), optional :: ocean_depth

 
 integer :: stdlogunit, stderrunit
 integer :: gni, gnj ! Global extent of ocean grid

  ! Get the stderr and stdlog unit numbers
 stderrunit=stderr()
 stdlogunit=stdlog()
 write(stdlogunit,*) "particles: "//trim(version)

 gni = Grid%ieg - Grid%isg + 1
 gnj = Grid%jeg - Grid%jsg + 1

 call particles_framework_init(parts, &
             gni, gnj, Grid%Domain%layout, Grid%Domain%io_layout, axes, Grid%Domain%X_FLAGS, Grid%Domain%X_FLAGS, &
             dt, Time, Grid%geolonT, Grid%geolatT, Grid%mask2dT, Grid%dxT, Grid%dyT, Grid%areaT, &
             Grid%cos_rot, Grid%sin_rot, ocean_depth=Grid%bathyT)


! call mpp_clock_begin(parts%clock_ior)
! call particles_io_init(parts,io_layout)
 call read_restart_parts(parts,Time)
! call parts_chksum(parts, 'read_restart_particles')
! call mpp_clock_end(parts%clock_ior)

! if (really_debug) call print_parts(stderrunit,parts,'particles_init, initial status')


end subroutine particles_init


! ##############################################################################

subroutine interp_flds(grd, i, j, xi, yj, uo, vo)
! Arguments
 type(particles_gridded), pointer :: grd
 integer, intent(in) :: i, j
 real, intent(in) :: xi, yj
 real, intent(out) :: uo, vo
 ! Local variables
 real :: cos_rot, sin_rot
 real :: hxm, hxp


 cos_rot=bilin(grd, grd%cos, i, j, xi, yj) ! If true, uses the inverted bilin function
 sin_rot=bilin(grd, grd%sin, i, j, xi, yj)

 uo=bilin(grd, grd%uo, i, j, xi, yj)
 vo=bilin(grd, grd%vo, i, j, xi, yj)
  
 ! Rotate vectors from local grid to lat/lon coordinates
 call rotate(uo, vo, cos_rot, sin_rot)

 
contains


 subroutine rotate(u, v, cos_rot, sin_rot)
   ! Arguments
   real, intent(inout) :: u, v
   real, intent(in) :: cos_rot, sin_rot
   ! Local variables
   real :: u_old, v_old

   u_old=u
   v_old=v
   u=cos_rot*u_old+sin_rot*v_old
   v=cos_rot*v_old-sin_rot*u_old

 end subroutine rotate

end subroutine interp_flds

!######################################################################################

subroutine accel(parts, part, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, ax, ay, axn, ayn, bxn, byn, debug_flag) !Saving  acceleration for Verlet, Adding Verlet flag - Alon  MP1
!subroutine accel(parts, part, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, ax, ay, debug_flag) !old version commmented out by Alon
! Arguments
type(particles), pointer :: parts
type(particle), pointer :: part
integer, intent(in) :: i, j
real, intent(in) :: xi, yj, lat, uvel, vvel, uvel0, vvel0, dt
real, intent(inout) :: ax, ay
real, intent(inout) :: axn, ayn, bxn, byn ! Added implicit and explicit accelerations to output -Alon
logical, optional :: debug_flag

ax=0
ay=0
axn=0
ayn=0
bxn=0
byn=0

end subroutine accel
! ##############################################################################

subroutine particles_run(parts, time, uo, vo, stagger)
! Arguments
 type(particles), pointer :: parts
 type(time_type), intent(in) :: time
 real, dimension(:,:), intent(in) :: uo, vo
 integer,    optional, intent(in) :: stagger
! Local variables
 integer :: iyr, imon, iday, ihr, imin, isec, k
 type(particles_gridded), pointer :: grd
 logical :: lerr, sample_traj, write_traj, lverbose
 real :: tmpsum
 integer :: i, j, Iu, ju, iv, Jv, Iu_off, ju_off, iv_off, Jv_off
 real :: mask
 real, dimension(:,:), allocatable :: uC_tmp, vC_tmp
 integer :: vel_stagger
 integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()

  call mpp_clock_begin(parts%clock)
  call mpp_clock_begin(parts%clock_int)

  vel_stagger = CGRID_NE ; if (present(stagger)) vel_stagger = stagger

  ! For convenience
  grd=>parts%grd

  call get_date(time, iyr, imon, iday, ihr, imin, isec)
  parts%current_year=iyr
  parts%current_yearday=yearday(imon, iday, ihr, imin, isec)
  ! Turn on sampling of trajectories, verbosity, budgets
  sample_traj=.false.
  if (parts%traj_sample_hrs>0) then
     if (mod(24*iday+ihr,parts%traj_sample_hrs).eq.0) sample_traj=.true.
  end if
  write_traj=.false.
  if (parts%traj_write_hrs>0) then
     if (mod(24*iday+ihr,parts%traj_write_hrs).eq.0) write_traj=.true.
  end if
  lverbose=.false.
  if (parts%verbose_hrs>0) then
     if (mod(24*iday+ihr,parts%verbose_hrs).eq.0) lverbose=verbose
  end if
  if (mpp_pe()==mpp_root_pe().and.lverbose) write(*,'(a,3i5,a,3i5,a,i5,f8.3)') &
       'parts: y,m,d=',iyr, imon, iday,' h,m,s=', ihr, imin, isec, &
       ' yr,yrdy=', parts%current_year, parts%current_yearday


  if (vel_stagger == BGRID_NE) then
    ! Copy ocean velocities. They are already on B-grid u-points.
    grd%uo(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1) = uo(:,:)
    grd%vo(grd%isc-1:grd%iec+1,grd%jsc-1:grd%jec+1) = vo(:,:)
    call mpp_update_domains(grd%uo, grd%vo, grd%domain, gridtype=BGRID_NE)
  elseif (vel_stagger == CGRID_NE) then
    ! The u- and v- points will have different offsets with symmetric memory.
    Iu_off = (size(uo,1) - (grd%iec - grd%isc))/2 - grd%isc + 1
    ju_off = (size(uo,2) - (grd%jec - grd%jsc))/2 - grd%jsc + 1
    iv_off = (size(vo,1) - (grd%iec - grd%isc))/2 - grd%isc + 1
    Jv_off = (size(vo,2) - (grd%jec - grd%jsc))/2 - grd%jsc + 1
    do I=grd%isc-1,grd%iec ; do J=grd%jsc-1,grd%jec
      ! Interpolate ocean velocities from C-grid velocity points.
      Iu = i + Iu_off ; ju = j + ju_off ; iv = i + iv_off ; Jv = j + Jv_off
      ! This masking is needed for now to prevent particles from running up on to land.
      mask = min(grd%msk(i,j), grd%msk(i+1,j), grd%msk(i,j+1), grd%msk(i+1,j+1))
      grd%uo(I,J) = mask * 0.5*(uo(Iu,ju)+uo(Iu,ju+1))
      grd%vo(I,J) = mask * 0.5*(vo(iv,Jv)+vo(iv+1,Jv))
    enddo ; enddo
  else
    call error_mesg('particle_run', 'Unrecognized value of stagger!', FATAL)
  endif

  call mpp_update_domains(grd%uo, grd%vo, grd%domain, gridtype=vel_stagger)

  if (debug) call parts_chksum(parts, 'run parts (top)')
  if (debug) call checksum_gridded(parts%grd, 'top of s/r run')

  call mpp_clock_end(parts%clock_int)

  ! Update particle positing
  call mpp_clock_begin(parts%clock_mom)
  !if (associated(parts%first)) call evolve_particles(parts) !Modify later
  if (debug) call parts_chksum(parts, 'run parts (evolved)',ignore_halo_violation=.true.)
  if (debug) call checksum_gridded(parts%grd, 's/r run after evolve')
  call mpp_clock_end(parts%clock_mom)

  ! Send parts to other PEs
  call mpp_clock_begin(parts%clock_com)
  call send_parts_to_other_pes(parts)
  if (debug) call parts_chksum(parts, 'run parts (exchanged)')
  if (debug) call checksum_gridded(parts%grd, 's/r run after exchange')
  call mpp_clock_end(parts%clock_com)

  ! For each part, record
  call mpp_clock_begin(parts%clock_dia)
  !if (sample_traj.and.associated(parts%first)) call record_posn(parts) ! modify later
  if (write_traj) then
    call move_all_trajectories(parts)
    call write_trajectory(parts%trajectories)
  endif

  ! Gridded diagnostics
  !if (grd%id_uo>0) &
  !  lerr=send_data(grd%id_uo, grd%uo(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  !if (grd%id_vo>0) &
  !  lerr=send_data(grd%id_vo, grd%vo(grd%isc:grd%iec,grd%jsc:grd%jec), Time)

  ! Dump particles to screen
  if (really_debug) call print_parts(stderrunit,parts,'particles_run, status')
  call mpp_clock_end(parts%clock_dia)

  call mpp_clock_end(parts%clock_int)


  if (debug) call parts_chksum(parts, 'run parts (bot)')
  if (debug) call checksum_gridded(parts%grd, 'end of s/r run')
  call mpp_clock_end(parts%clock_dia)

  call mpp_clock_end(parts%clock)

 
end subroutine particles_run


! ##############################################################################
!> Evolves particles forward by updating velocity and position with a time-stepping scheme
subroutine evolve_particles(parts)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  ! Local variables
  type(particles_gridded), pointer :: grd
  type(particle), pointer :: part
  real :: uveln, vveln, lonn, latn
  real :: axn, ayn, bxn, byn          ! Added by Alon - explicit and implicit accelerations from the previous step
  real :: xi, yj
  integer :: i, j
  integer :: grdi, grdj
  integer :: stderrunit
  logical :: bounced, interactive_particles_on, Runge_not_Verlet

  

end subroutine evolve_particles




! ##############################################################################

subroutine particles_save_restart(parts)
! Arguments
type(particles), pointer :: parts
! Local variables

  if (.not.associated(parts)) return

  call mpp_clock_begin(parts%clock_iow)
  call parts_chksum(parts, 'write_restart parts')
  call write_restart(parts)
  call mpp_clock_end(parts%clock_iow)

end subroutine particles_save_restart

! ##############################################################################

subroutine particles_end(parts)
! Arguments
type(particles), pointer :: parts
! Local variables
type(particle), pointer :: this, next

  if (.not.associated(parts)) return

  ! particles_save_restart() is called directly by SIS1 and SIS2 so
  ! we do not need to call it a second time. If particles were controlled
  ! by the coupler then the particles would need to take responsibility for
  ! the restarts at the end of the run.
  !call particles_save_restart(parts)

  call mpp_clock_begin(parts%clock_ini)
  ! Delete parts and structures
  call move_all_trajectories(parts, delete_parts=.true.)

  call write_trajectory(parts%trajectories)

  deallocate(parts%grd%lon)
  deallocate(parts%grd%lat)
  deallocate(parts%grd%lonc)
  deallocate(parts%grd%latc)
  deallocate(parts%grd%dx)
  deallocate(parts%grd%dy)
  !deallocate(parts%grd%area)
  deallocate(parts%grd%msk)
  deallocate(parts%grd%cos)
  deallocate(parts%grd%sin)
  deallocate(parts%grd%ocean_depth)
  deallocate(parts%grd%uo)
  deallocate(parts%grd%vo)
  deallocate(parts%grd%domain)
  deallocate(parts%grd)
  call dealloc_buffer(parts%obuffer_n)
  call dealloc_buffer(parts%obuffer_s)
  call dealloc_buffer(parts%obuffer_e)
  call dealloc_buffer(parts%obuffer_w)
  call dealloc_buffer(parts%ibuffer_n)
  call dealloc_buffer(parts%ibuffer_s)
  call dealloc_buffer(parts%ibuffer_e)
  call dealloc_buffer(parts%ibuffer_w)
  call dealloc_buffer(parts%ibuffer_io)
  call dealloc_buffer(parts%ibuffer_io)
  call mpp_clock_end(parts%clock_ini)
  deallocate(parts)

  if (mpp_pe()==mpp_root_pe()) write(*,'(a,i8)') 'diamonds: particles_end complete',mpp_pe()

  contains

  subroutine dealloc_buffer(buff)
  ! Arguments
  type(buffer), pointer :: buff
  ! Local variables
    if (associated(buff)) then
      if (associated(buff%data)) deallocate(buff%data)
      deallocate(buff)
    endif
  end subroutine dealloc_buffer

end subroutine particles_end



end module particles_mod
