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
             gni, gnj, Grid%Domain%layout, Grid%Domain%io_layout, axes, Grid%Domain%dom_x_flags, Grid%Domain%dom_y_flags, &
             dt, Time, Grid%geolonT, Grid%geolatT, Grid%mask2dT, Grid%dxT, Grid%dyT, Grid%areaT, &
             Grid%cos_rot, Grid%sin_rot, ocean_depth=Grid%bathyT)


! call mpp_clock_begin(parts%clock_ior)
! call particles_io_init(parts,io_layout)
 call read_restart_parts(parts,Grid,Time)
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
!subroutine accel(bergs, berg, i, j, xi, yj, lat, uvel, vvel, uvel0, vvel0, dt, ax, ay, debug_flag) !old version commmented out by Alon
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
  if (associated(parts%first)) call evolve_particles(parts)
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
  if (sample_traj.and.associated(parts%first)) call record_posn(parts)
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

subroutine evolve_particles(parts)
  use particles_framework, only: pi_180
   
  ! Arguments
  type(particles), pointer :: parts
  ! Local variables
  type(particles_gridded), pointer :: grd
  real :: uvel1, vvel1, lon1, lat1, u1, v1, dxdl1, ax1, ay1, axn1, ayn1 
  real :: uvel2, vvel2, lon2, lat2, u2, v2, dxdl2, ax2, ay2, axn2, ayn2
  real :: uvel3, vvel3, lon3, lat3, u3, v3, dxdl3, ax3, ay3, axn3, ayn3
  real :: uvel4, vvel4, lon4, lat4, u4, v4, dxdl4, ax4, ay4, axn4, ayn4
  real :: uveln, vveln, lonn, latn, un, vn, dxdln
  real :: x1, xdot1, xddot1, y1, ydot1, yddot1, xddot1n, yddot1n 
  real :: x2, xdot2, xddot2, y2, ydot2, yddot2, xddot2n, yddot2n
  real :: x3, xdot3, xddot3, y3, ydot3, yddot3, xddot3n, yddot3n
  real :: x4, xdot4, xddot4, y4, ydot4, yddot4, xddot4n, yddot4n
  real :: xn, xdotn, yn, ydotn, xddotn, yddotn
  real :: bxddot, byddot                                               ! Added by Alon
  real :: axn, ayn, bxn, byn                                           ! Added by Alon - explicit and implicit accelations from the previous step
  real :: r180_pi, dt, dt_2, dt_6, dydl, Rearth
  integer :: i, j
  integer :: i1,j1,i2,j2,i3,j3,i4,j4
  real :: xi, yj
  logical :: bounced, on_tangential_plane, error_flag
  logical :: Runge_not_Verlet  ! Runge_not_Verlet=1 for Runge Kutta, =0 for Verlet method. Added by Alon
  type(particle), pointer :: part
  integer :: stderrunit


  ! 4th order Runge-Kutta to solve:
  !    d/dt X = V,  d/dt V = A
  ! with I.C.'s:
  !    X=X1 and V=V1
  !
  !  A1 = A(X1)
  !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
  !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
  !  X4 = X1+  dt*V3 ; V4 = V1+  dt*A3; A4=A(X4)
  !
  !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
  !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6
  
  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>parts%grd

  ! Common constants
  r180_pi=1./pi_180
  dt=parts%dt
  dt_2=0.5*dt
  dt_6=dt/6.
  Rearth=6360.e3

  !Choosing time stepping scheme - Alon
  !Runge_not_Verlet=.False.    !Loading manually: true=Runge Kutta, False=Verlet   , Alon
  Runge_not_Verlet=parts%Runge_not_Verlet  ! Loading directly from namelist/default , Alon

  part=>parts%first
  do while (associated(part)) ! loop over all parts

  if (.not. is_point_in_cell(parts%grd, part%lon, part%lat, part%ine, part%jne) ) then
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
    enddo
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
    enddo
    call print_part(stderrunit, part, 'evolve_particle, part is not in proper starting cell')
    write(stderrunit,'(a,i3,2(i4,3f8.2))') 'evolve_particle: pe,lon/lat(i,j)=', mpp_pe(), &
             part%ine,part%lon,grd%lon(part%ine-1,part%jne-1),grd%lon(part%ine,part%jne), &
             part%jne,part%lat,grd%lat(part%ine-1,part%jne-1),grd%lat(part%ine,part%jne)
    if (debug) call error_mesg('diamonds, evolve_particle','part is in wrong starting cell!',FATAL)
  endif

  if (debug) call check_position(grd, part, 'evolve_particle (top)')

  i=part%ine
  j=part%jne
  xi=part%xi
  yj=part%yj
  bounced=.false.
  on_tangential_plane=.false.
  if (part%lat>89.) on_tangential_plane=.true. ! This could be a problem for some grid MJH
  i1=i;j1=j

  if (Runge_not_Verlet) then !Start of the Runge-Kutta Loop -Added by Alon, MP2

 !Loading past acceleartions - Alon
  axn=part%axn; ayn=part%ayn !Alon
  axn1=axn; axn2=axn; axn3=axn; axn4=axn
  ayn1=ayn; ayn2=ayn; ayn3=ayn; ayn4=ayn



 ! A1 = A(X1)
  lon1=part%lon; lat1=part%lat
  if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1)
  dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  dydl=r180_pi/Rearth
  uvel1=part%uvel; vvel1=part%vvel
  if (on_tangential_plane) call rotvec_to_tang(lon1,uvel1,vvel1,xdot1,ydot1)
  u1=uvel1*dxdl1; v1=vvel1*dydl
  call accel(parts, part, i, j, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn1, ayn1, bxn, byn) !axn,ayn, bxn, byn  - Added by Alon
  !call accel(parts, part, i, j, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt, ax1, ay1, axn1, ayn1, bxn, byn) !Note change to dt. Markpoint_1
  if (on_tangential_plane) call rotvec_to_tang(lon1,ax1,ay1,xddot1,yddot1)
  if (on_tangential_plane) call rotvec_to_tang(lon1,axn1,ayn1,xddot1n,yddot1n) !Alon
  
  !  X2 = X1+dt/2*V1 ; V2 = V1+dt/2*A1; A2=A(X2)
  !if (debug) write(stderr(),*) 'diamonds, evolve: x2=...'
  if (on_tangential_plane) then
    x2=x1+dt_2*xdot1; y2=y1+dt_2*ydot1
    xdot2=xdot1+dt_2*xddot1; ydot2=ydot1+dt_2*yddot1
    call rotpos_from_tang(x2,y2,lon2,lat2)
    call rotvec_from_tang(lon2,xdot2,ydot2,uvel2,vvel2)
  else
    lon2=lon1+dt_2*u1; lat2=lat1+dt_2*v1
    uvel2=uvel1+dt_2*ax1; vvel2=vvel1+dt_2*ay1
  endif
  i=i1;j=j1;xi=part%xi;yj=part%yj
  call adjust_index_and_ground(grd, lon2, lat2, uvel2, vvel2, i, j, xi, yj, bounced, error_flag)
  i2=i; j2=j

  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon2,lat2,x2,y2)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(parts%grd, lon2, lat2, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   !call print_fld(grd, grd%msk, 'msk')
   !call print_fld(grd, grd%ssh, 'ssh')
   !call print_fld(grd, grd%sst, 'sst')
   !call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: i1,i2,i=',i1,i2,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: j1,j2,j=',j1,j2,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lon1,lon2=',lon1,lon2,part%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lat1,lat2=',lat1,lat2,part%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: u1,u2,u0=',uvel1,uvel2,part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: v1,v2,v0=',vvel1,vvel2,part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ax1,ax2=',dt*ax1,dt*ax2
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ay1,ay2=',dt*ay1,dt*ay2
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u1,u2,u0=',dt*uvel1,dt*uvel2,dt*part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v1,v2,v0=',dt*vvel1,dt*vvel2,dt*part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u1,u2 (deg)=',dt*u1,dt*u2
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v1,v2 (deg)=',dt*v1,dt*v2
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(parts, part, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn1, ayn1, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn,- Added by Alon
   call print_part(stderrunit, part, 'evolve_particle, out of position at 2')
   write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos2 i,j,lon,lat,xi,yj=',i,j,lon2,lat2,xi,yj
   write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos2 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
   bounced=is_point_in_cell(parts%grd, lon2, lat2, i, j, explain=.true.)
   call error_mesg('diamonds, evolve_particle','part is out of posn at 2!',FATAL)
  endif
  dxdl2=r180_pi/(Rearth*cos(lat2*pi_180))
  u2=uvel2*dxdl2; v2=vvel2*dydl
  call accel(parts, part, i, j, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, axn2, ayn2, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
  !call accel(parts, part, i, j, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt, ax2, ay2, axn2, ayn2, bxn, byn) !Note change to dt. Markpoint_1
  if (on_tangential_plane) call rotvec_to_tang(lon2,ax2,ay2,xddot2,yddot2)
  if (on_tangential_plane) call rotvec_to_tang(lon2,axn2,ayn2,xddot2n,yddot2n) !Alon
  
  !  X3 = X1+dt/2*V2 ; V3 = V1+dt/2*A2; A3=A(X3)
 !if (debug) write(stderr(),*) 'diamonds, evolve: x3=...'
  if (on_tangential_plane) then
    x3=x1+dt_2*xdot2; y3=y1+dt_2*ydot2
    xdot3=xdot1+dt_2*xddot2; ydot3=ydot1+dt_2*yddot2
    call rotpos_from_tang(x3,y3,lon3,lat3)
    call rotvec_from_tang(lon3,xdot3,ydot3,uvel3,vvel3)
  else
    lon3=lon1+dt_2*u2; lat3=lat1+dt_2*v2
    uvel3=uvel1+dt_2*ax2; vvel3=vvel1+dt_2*ay2
  endif
  i=i1;j=j1;xi=part%xi;yj=part%yj
  call adjust_index_and_ground(grd, lon3, lat3, uvel3, vvel3, i, j, xi, yj, bounced, error_flag)
  i3=i; j3=j
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon3,lat3,x3,y3)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(parts%grd, lon3, lat3, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   !call print_fld(grd, grd%msk, 'msk')
   !call print_fld(grd, grd%ssh, 'ssh')
   !call print_fld(grd, grd%sst, 'sst')
   !call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: i1,i2,i3,i=',i1,i2,i3,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: j1,j2,j3,j=',j1,j2,j3,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lon1,lon2,lon3=',lon1,lon2,lon3,part%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lat1,lat2,lat3=',lat1,lat2,lat3,part%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: u1,u2,u3,u0=',uvel1,uvel2,uvel3,part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: v1,v2,v3,v0=',vvel1,vvel2,vvel3,part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ax1,ax2,ax3=',dt*ax1,dt*ax2,dt*ax3
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ay1,ay2,ay3=',dt*ay1,dt*ay2,dt*ay3
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u1,u2,u3,u0=',dt*uvel1,dt*uvel2,dt*uvel3,dt*part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v1,v2,v3,v0=',dt*vvel1,dt*vvel2,dt*vvel3,dt*part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u1,u2,u3 (deg)=',dt*u1,dt*u2,dt*u3
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v1,v2,v3 (deg)=',dt*v1,dt*v2,dt*v3
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(parts, part, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn1, ayn1, bxn, byn, debug_flag=.true.) !axn, ayn - Added by Alon
   write(stderrunit,*) 'Acceleration terms for position 2'
   error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
   call accel(parts, part, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, axn2, ayn2, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
    call print_part(stderrunit, part, 'evolve_particle, out of position at 3')
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos3 i,j,lon,lat,xi,yj=',i,j,lon3,lat3,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos3 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    bounced=is_point_in_cell(parts%grd, lon2, lat2, i, j, explain=.true.)
    call error_mesg('diamonds, evolve_particle','part is out of posn at 3!',FATAL)
  endif
  dxdl3=r180_pi/(Rearth*cos(lat3*pi_180))
  u3=uvel3*dxdl3; v3=vvel3*dydl
  call accel(parts, part, i, j, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, ax3, ay3, axn3, ayn3, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
  if (on_tangential_plane) call rotvec_to_tang(lon3,ax3,ay3,xddot3,yddot3)
  if (on_tangential_plane) call rotvec_to_tang(lon3,axn3,ayn3,xddot3n,yddot3n) !Alon
  
  !  X4 = X1+dt*V3 ; V4 = V1+dt*A3; A4=A(X4)
 !if (debug) write(stderr(),*) 'diamonds, evolve: x4=...'
  if (on_tangential_plane) then
    x4=x1+dt*xdot3; y4=y1+dt*ydot3
    xdot4=xdot1+dt*xddot3; ydot4=ydot1+dt*yddot3
    call rotpos_from_tang(x4,y4,lon4,lat4)
    call rotvec_from_tang(lon4,xdot4,ydot4,uvel4,vvel4)
  else
    lon4=lon1+dt*u3; lat4=lat1+dt*v3
    uvel4=uvel1+dt*ax3; vvel4=vvel1+dt*ay3
  endif
  i=i1;j=j1;xi=part%xi;yj=part%yj
  call adjust_index_and_ground(grd, lon4, lat4, uvel4, vvel4, i, j, xi, yj, bounced, error_flag)
  i4=i; j4=j
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon4,lat4,x4,y4)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(parts%grd, lon4, lat4, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   !call print_fld(grd, grd%msk, 'msk')
   !call print_fld(grd, grd%ssh, 'ssh')
   !call print_fld(grd, grd%sst, 'sst')
   !call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: i1,i2,i3,i4,i=',i1,i2,i3,i4,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: j1,j2,j3,j4,j=',j1,j2,j3,j4,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lon1,lon2,lon3,lon4=',lon1,lon2,lon3,lon4,part%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lat1,lat2,lat3,lat4=',lat1,lat2,lat3,lat4,part%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: u1,u2,u3,u4,u0=',uvel1,uvel2,uvel3,uvel4,part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: v1,v2,v3,v4,v0=',vvel1,vvel2,vvel3,vvel4,part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ax1,ax2,ax3,ax4=',dt*ax1,dt*ax2,dt*ax3,dt*ax4
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ay1,ay2,ay3,ay4=',dt*ay1,dt*ay2,dt*ay3,dt*ay4
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u1,u2,u3,u4,u0=',dt*uvel1,dt*uvel2,dt*uvel3,dt*uvel4,dt*part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v1,v2,v3,v4,v0=',dt*vvel1,dt*vvel2,dt*vvel3,dt*vvel4,dt*part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u1,u2,u3,u4 (deg)=',dt*u1,dt*u2,dt*u3,dt*u4
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v1,v2,v3,v4 (deg)=',dt*v1,dt*v2,dt*v3,dt*v4
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(parts, part, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn1, ayn1, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
   write(stderrunit,*) 'Acceleration terms for position 2'
   error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
   call accel(parts, part, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, axn2, ayn2, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
   write(stderrunit,*) 'Acceleration terms for position 3'
   error_flag=pos_within_cell(grd, lon3, lat3, i3, j3, xi, yj)
   call accel(parts, part, i3, j3, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, ax3, ay3, axn3, ayn3, bxn, byn, debug_flag=.true.) !axn, ayn, bxn, byn - Added by Alon
    call print_part(stderrunit, part, 'evolve_particle, out of position at 4')
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'pos4 i,j,lon,lat,xi,yj=',i,j,lon4,lat4,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'pos4 box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    bounced=is_point_in_cell(parts%grd, lon2, lat2, i, j, explain=.true.)
    call error_mesg('diamonds, evolve_particle','part is out of posn at 4!',FATAL)
  endif
  dxdl4=r180_pi/(Rearth*cos(lat4*pi_180))
  u4=uvel4*dxdl4; v4=vvel4*dydl
  call accel(parts, part, i, j, xi, yj, lat4, uvel4, vvel4, uvel1, vvel1, dt, ax4, ay4, axn4, ayn4, bxn, byn) !axn, ayn, bxn, byn - Added by Alon
  if (on_tangential_plane) call rotvec_to_tang(lon4,ax4,ay4,xddot4,yddot4)
  if (on_tangential_plane) call rotvec_to_tang(lon4,axn4,ayn4,xddot4n,yddot4n)
  
  !  Xn = X1+dt*(V1+2*V2+2*V3+V4)/6
  !  Vn = V1+dt*(A1+2*A2+2*A3+A4)/6
  if (on_tangential_plane) then
    xn=x1+dt_6*( (xdot1+xdot4)+2.*(xdot2+xdot3) )
    yn=y1+dt_6*( (ydot1+ydot4)+2.*(ydot2+ydot3) )
    xdotn=xdot1+dt_6*( (xddot1+xddot4)+2.*(xddot2+xddot3) )
    ydotn=ydot1+dt_6*( (yddot1+yddot4)+2.*(yddot2+yddot3) )
    xddotn=( (xddot1n+xddot4n)+2.*(xddot2n+xddot3n) )/6.  !Alon  
    yddotn=( (yddot1n+yddot4n)+2.*(yddot2n+yddot3n) )/6.  !Alon  
    call rotpos_from_tang(xn,yn,lonn,latn)
    call rotvec_from_tang(lonn,xdotn,ydotn,uveln,vveln)
    call rotvec_from_tang(lonn,xddotn,yddotn,axn,ayn) !Alon
  else
    lonn=part%lon+dt_6*( (u1+u4)+2.*(u2+u3) )
    latn=part%lat+dt_6*( (v1+v4)+2.*(v2+v3) )
    uveln=part%uvel+dt_6*( (ax1+ax4)+2.*(ax2+ax3) )
    vveln=part%vvel+dt_6*( (ay1+ay4)+2.*(ay2+ay3) )
    axn=( (axn1+axn4)+2.*(axn2+axn3) )/6. !Alon
    ayn=( (ayn1+ayn4)+2.*(ayn2+ayn3) )/6. !Alon
    bxn=(((ax1+ax4)+2.*(ax2+ax3) )/6)  - (axn/2)
    byn=(((ay1+ay4)+2.*(ay2+ay3) )/6)  - (ayn/2)

  endif



  i=i1;j=j1;xi=part%xi;yj=part%yj
  call adjust_index_and_ground(grd, lonn, latn, uveln, vveln, i, j, xi, yj, bounced, error_flag)

  if (.not.error_flag) then
    if (.not. is_point_in_cell(parts%grd, lonn, latn, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   !call print_fld(grd, grd%msk, 'msk')
   !call print_fld(grd, grd%ssh, 'ssh')
   !call print_fld(grd, grd%sst, 'sst')
   !call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: i1,i2,i3,i4,i=',i1,i2,i3,i4,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: j1,j2,j3,j4,j=',j1,j2,j3,j4,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lon1,lon2,lon3,lon4,lonn=',lon1,lon2,lon3,lon4,lonn,part%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lat1,lat2,lat3,lat4,latn=',lat1,lat2,lat3,lat4,latn,part%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: u1,u2,u3,u4,un,u0=',uvel1,uvel2,uvel3,uvel4,uveln,part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: v1,v2,v3,v4,vn,v0=',vvel1,vvel2,vvel3,vvel4,vveln,part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ax1,ax2,ax3,ax4,axn=',&
        & dt*ax1,dt*ax2,dt*ax3,dt*ax4,dt_6*( (ax1+ax4)+2.*(ax2+ax3) )
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ay1,ay2,ay3,ay4,ayn=',&
        & dt*ay1,dt*ay2,dt*ay3,dt*ay4,dt_6*( (ay1+ay4)+2.*(ay2+ay3) )
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u1,u2,u3,u4,un,u0=',&
        & dt*uvel1,dt*uvel2,dt*uvel3,dt*uvel4,dt*uveln,dt*part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v1,v2,v3,v4,vn,v0=',&
        & dt*vvel1,dt*vvel2,dt*vvel3,dt*vvel4,dt*vveln,dt*part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u1,u2,u3,u4,u_rk (deg)=',&
        & dt*u1,dt*u2,dt*u3,dt*u4,dt_6*( (u1+u4)+2.*(u2+u3) )
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v1,v2,v3,v4,v_rk (deg)=',&
        & dt*v1,dt*v2,dt*v3,dt*v4,dt_6*( (v1+v4)+2.*(v2+v3) )
   write(stderrunit,*) 'diamonds, evolve_particle: on_tangential_plane=',on_tangential_plane
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(parts, part, i1, j1, xi, yj, lat1, uvel1, vvel1, uvel1, vvel1, dt_2, ax1, ay1, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
   write(stderrunit,*) 'Acceleration terms for position 2'
   error_flag=pos_within_cell(grd, lon2, lat2, i2, j2, xi, yj)
   call accel(parts, part, i2, j2, xi, yj, lat2, uvel2, vvel2, uvel1, vvel1, dt_2, ax2, ay2, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
   write(stderrunit,*) 'Acceleration terms for position 3'
   error_flag=pos_within_cell(grd, lon3, lat3, i3, j3, xi, yj)
   call accel(parts, part, i3, j3, xi, yj, lat3, uvel3, vvel3, uvel1, vvel1, dt, ax3, ay3, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
   write(stderrunit,*) 'Acceleration terms for position 4'
   error_flag=pos_within_cell(grd, lon4, lat4, i4, j4, xi, yj)
   call accel(parts, part, i4, j4, xi, yj, lat4, uvel4, vvel4, uvel1, vvel1, dt, ax4, ay4, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn - Added by Alon
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'posn i,j,lon,lat,xi,yj=',i,j,lonn,latn,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'posn box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    call print_part(stderrunit, part, 'evolve_particle, out of cell at end!')
    bounced=is_point_in_cell(parts%grd, lonn, latn, i, j, explain=.true.)
    if (debug) call error_mesg('diamonds, evolve_particle','part is out of posn at end!',FATAL)
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
    enddo
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
    enddo
  endif
 
  endif ! End of the Runge-Kutta Loop -added by Alon  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (.not.Runge_not_Verlet) then !Start of the Verlet time_stepping -Whole loop added by Alon

 ! In this scheme a_n and b_n are saved from the previous timestep, giving the explicit and implicit parts of the acceleration, and a_np1, b_np1 are for the next time step
 ! Note that ax1=a_np1/2 +b_np1, as calculated by the acceleration subrouting
 ! Positions and velocity is updated by
 ! X2 = X1+dt*V1+((dt^2)/2)*a_n +((dt^2)/2)*b_n = X1+dt*u_star +((dt^2)/2)*b_n 
 ! V2 = V1+dt/2*a_n +dt/2*a_np1 +dt*b_n = u_star + dt/2*a_np1 + dt*b_np1 = u_star +dt*ax

!print *, 'you are here!'

  lon1=part%lon; lat1=part%lat
  if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1)
  dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  dydl=r180_pi/Rearth
  uvel1=part%uvel; vvel1=part%vvel

!Loading past acceleartions - Alon
  axn=part%axn; ayn=part%ayn !Alon
  bxn=part%bxn; byn=part%byn !Alon

  print *, 'first', axn, bxn, lon1, lat1, uvel1, i, j ,xi, yj

! Velocities used to update the position
  uvel2=uvel1+(dt_2*axn)+(dt_2*bxn)                    !Alon
  vvel2=vvel1+(dt_2*ayn)+(dt_2*byn)                    !Alon

  if (on_tangential_plane) call rotvec_to_tang(lon1,uvel2,vvel2,xdot2,ydot2)
  u2=uvel2*dxdl1; v2=vvel2*dydl


!Solving for new position
  if (on_tangential_plane) then
    xn=x1+(dt*xdot2) ; yn=y1+(dt*ydot2)             !Alon
    call rotpos_from_tang(xn,yn,lonn,latn)
  else
    lonn=lon1+(dt*u2) ; latn=lat1+(dt*v2)  !Alon
  endif
  dxdln=r180_pi/(Rearth*cos(latn*pi_180))

! Turn the velocities into u_star, v_star.(uvel3 is v_star)
  uvel3=uvel1+(dt_2*axn)                  !Alon
  vvel3=vvel1+(dt_2*ayn)                  !Alon


!Adjusting mass...                      Alon decided to move this before calculating the new velocities (so that acceleration can be a fn(r_np1)
  i=i1;j=j1;xi=part%xi;yj=part%yj
  call adjust_index_and_ground(grd, lonn, latn, uvel3, vvel3, i, j, xi, yj, bounced, error_flag)  !Alon:"unclear which velocity to use here?"
!  call adjust_index_and_ground(grd, lonn, latn, uvel1, vvel1, i, j, xi, yj, bounced, error_flag)  !Alon:"unclear which velocity to use here?"

  if (bounced) then  !This is the case when the particle changes direction due to  topography
    axn=0.
    ayn=0.
    bxn=0.
    byn=0.
  endif


  i2=i; j2=j

  print *, 'second', axn, bxn, lon1, lat1, uvel1, i , j , xi, yj
!Calling the acceleration   (note that the velocity is converted to u_star inside the accel script)
  call accel(parts, part, i, j, xi, yj, latn, uvel1, vvel1, uvel1, vvel1, dt, ax1, ay1, axn, ayn, bxn, byn) !axn, ayn, bxn, byn - Added by Alon

  print *, 'third', axn, bxn, lon1, lat1, uvel1, i, j, xi, yj
!Solving for the new velocity
  if (on_tangential_plane) then
    call rotvec_to_tang(lonn,uvel3,vvel3,xdot3,ydot3)
    call rotvec_to_tang(lon1,ax1,ay1,xddot1,yddot1)
    xdotn=xdot3+(dt*xddot1); ydotn=ydot3+(dt*yddot1)                                    !Alon
    call rotvec_from_tang(lonn,xdotn,ydotn,uveln,vveln)
  else
    uvel4=uvel3+(dt*ax1); vvel4=vvel3+(dt*ay1)    !Alon , we call it uvel3, vvel3 until it is put into lat/long co-ordinates, where it becomes uveln, vveln
  endif
!  uveln=uvel4*dxdln; vveln=vvel4*dydl    !Converted to degrees.  (Perhaps this should not be here)
  uveln=uvel4
  vveln=vvel4 

  print *, 'forth', axn, bxn, lon1, lat1, uvel1, i, j, xi, yj, uveln
!Debugging
  if (.not.error_flag) then
    if (.not. is_point_in_cell(parts%grd, lonn, latn, i, j)) error_flag=.true.
  endif
  if (error_flag) then
   !call print_fld(grd, grd%msk, 'msk')
   !call print_fld(grd, grd%ssh, 'ssh')
   !call print_fld(grd, grd%sst, 'sst')
   !call print_fld(grd, grd%hi, 'hi')
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: isd,isc,iec,ied=',grd%isd,grd%isc,grd%iec,grd%ied
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: jsd,jsc,jec,jed=',grd%jsd,grd%jsc,grd%jec,grd%jed
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: i1,i2,i=',i1,i2,i
   write(stderrunit,'(a,6i5)') 'diamonds, evolve_particle: j1,j2,j=',j1,j2,j
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lon1,lonn=',lon1,lonn,part%lon
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: lat1,latn=',lat1,latn,part%lat
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: u3,un,u0=',uvel3,uveln,part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: v3,vn,v0=',vvel3,vveln,part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ax1=',&
        & dt*ax1
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* ay1=',&
        & dt*ay1
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u3,un,u0=',&
        & dt*uvel3,dt*uveln,dt*part%uvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v3,vn,v0=',&
        & dt*vvel3,dt*vveln,dt*part%vvel
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* u1,u_n (deg)=',&
        & dt*u1,dt*uveln
   write(stderrunit,'(a,6es9.3)') 'diamonds, evolve_particle: dt* v1,v_n (deg)=',&
        & dt*v1,dt*vveln
   write(stderrunit,*) 'diamonds, evolve_particle: on_tangential_plane=',on_tangential_plane
   write(stderrunit,*) 'Acceleration terms for position 1'
   error_flag=pos_within_cell(grd, lon1, lat1, i1, j1, xi, yj)
   call accel(parts, part, i2, j2, xi, yj, latn, uvel3, vvel3, uvel1, vvel1, dt_2, ax1, ay1, axn, ayn, bxn, byn, debug_flag=.true.)  !axn, ayn, bxn, byn - Added by Alon
   
    write(stderrunit,'(a,i3,a,2i4,4f8.3)') 'pe=',mpp_pe(),'posn i,j,lon,lat,xi,yj=',i,j,lonn,latn,xi,yj
    write(stderrunit,'(a,i3,a,4f8.3)') 'pe=',mpp_pe(),'posn box=',grd%lon(i-1,j-1),grd%lon(i,j),grd%lat(i-1,j-1),grd%lat(i,j)
    call print_part(stderrunit, part, 'evolve_particle, out of cell at end!')
    bounced=is_point_in_cell(parts%grd, lonn, latn, i, j, explain=.true.)
    if (debug) call error_mesg('diamonds, evolve_particle','part is out of posn at end!',FATAL)
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lon',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lon(i,j),i=grd%isd,grd%ied)
    enddo
    write(stderrunit,'(i4,a4,32i7)') mpp_pe(),'Lat',(i,i=grd%isd,grd%ied)
    do j=grd%jed,grd%jsd,-1
      write(stderrunit,'(2i4,32f7.1)') mpp_pe(),j,(grd%lat(i,j),i=grd%isd,grd%ied)
    enddo
  endif


  print *, 'fifth', axn, bxn, lon1, lat1, uvel1, i, j, xi, yj, uveln

  endif ! End of the Verlet Stepiing -added by Alon  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Saving all the particle variables.
  part%axn=axn !Alon
  part%ayn=ayn !Alon
  part%bxn=bxn !Alon
  part%byn=byn !Alon
  part%lon=lonn
  part%lat=latn
  part%uvel=uveln
  part%vvel=vveln
  part%ine=i
  part%jne=j
  part%xi=xi
  part%yj=yj
  !call interp_flds(grd, i, j, xi, yj, part%uo, part%vo, part%ui, part%vi, part%ua, part%va, part%ssh_x, part%ssh_y, part%sst)
  !if (debug) call print_part(stderr(), part, 'evolve_particle, final posn.')
  if (debug) call check_position(grd, part, 'evolve_particle (bot)')


  part=>part%next
  enddo ! loop over all parts

! When we are using interactive particles, we update the (old) particle positions and velocities in a second loop, all together (to make code order invarient)
  part=>parts%first

  do while (associated(part)) ! loop over all parts

      !Updating particle positions and velocities
      part%lon_old=part%lon
      part%lat_old=part%lat
      part%uvel_old=part%uvel
      part%vvel_old=part%vvel

      part=>part%next
  enddo ! loop over all parts



  contains

  subroutine rotpos_to_tang(lon, lat, x, y)
  ! Arguments
  real, intent(in) :: lon, lat
  real, intent(out) :: x, y
  ! Local variables
  real :: r,colat,clon,slon

    if (lat>90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat>90 already!',lat
      call error_mesg('diamonds, rotpos_to_tang','Something went very wrong!',FATAL)
    endif
    if (lat==90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat==90 already!',lat
      call error_mesg('diamonds, rotpos_to_tang','Something went wrong!',FATAL)
    endif

    colat=90.-lat
    r=Rearth*(colat*pi_180)
    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    x=r*clon
    y=r*slon

  end subroutine rotpos_to_tang

  subroutine rotpos_from_tang(x, y, lon, lat)
  ! Arguments
  real, intent(in) :: x, y
  real, intent(out) :: lon, lat
  ! Local variables
  real :: r

    r=sqrt(x**2+y**2)
    lat=90.-(r180_pi*r/Rearth)
    lon=r180_pi*acos(x/r)*sign(1.,y)

  end subroutine rotpos_from_tang

  subroutine rotvec_to_tang(lon, uvel, vvel, xdot, ydot)
  ! Arguments
  real, intent(in) :: lon, uvel, vvel
  real, intent(out) :: xdot, ydot
  ! Local variables
  real :: clon,slon

    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    xdot=-slon*uvel-clon*vvel
    ydot=clon*uvel-slon*vvel

  end subroutine rotvec_to_tang

  subroutine rotvec_from_tang(lon, xdot, ydot, uvel, vvel)
  ! Arguments
  real, intent(in) :: lon, xdot, ydot
  real, intent(out) :: uvel, vvel
  ! Local variables
  real :: clon,slon

    clon=cos(lon*pi_180)
    slon=sin(lon*pi_180)
    uvel=-slon*xdot+clon*ydot
    vvel=-clon*xdot-slon*ydot

  end subroutine rotvec_from_tang

! ##############################################################################

subroutine adjust_index_and_ground(grd, lon, lat, uvel, vvel, i, j, xi, yj, bounced, error)
! Arguments
type(particles_gridded), pointer :: grd
real, intent(inout) :: lon, lat, uvel, vvel, xi, yj
integer, intent(inout) :: i,j
logical, intent(out) :: bounced, error
! Local variables
logical lret, lpos
real, parameter :: posn_eps=0.05
integer :: icount, i0, j0, inm, jnm
real :: xi0, yj0, lon0, lat0

  bounced=.false.
  error=.false.
  lon0=lon; lat0=lat ! original position
  i0=i; j0=j ! original i,j
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  xi0=xi; yj0=yj ! original xi,yj
  if (debug) then
    !Sanity check lret, xi and yj
    lret=is_point_in_cell(grd, lon, lat, i, j)
    if (xi<0. .or. xi>1. .or. yj<0. .or. yj>1.) then
      if (lret) then
        write(stderrunit,*) 'diamonds, adjust: WARNING!!! lret=T but |xi,yj|>1',mpp_pe()
        write(stderrunit,*) 'diamonds, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'diamonds, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'diamonds, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'diamonds, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'diamonds, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'diamonds, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        error=.true.; return
     endif
    else
      if (.not.lret) then
        write(stderrunit,*) 'diamonds, adjust: WARNING!!! lret=F but |xi,yj|<1',mpp_pe()
        write(stderrunit,*) 'diamonds, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'diamonds, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'diamonds, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'diamonds, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'diamonds, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'diamonds, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        error=.true.; return
      endif
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  endif ! debug
  if (lret) return ! Part was already in cell

  ! Find inm, jnm (as if adjusting i,j) based on xi,yj
  ! ignoring the mand mask.
  ! NOTE:  This search appears to have *NO* active role
  ! in the algorithm other than to flag a warning.
  icount=0
  inm=i0; jnm=j0 ! original i,j
  do while (debug .and. .not.lret .and. icount<4)
    icount=icount+1
    if (xi.lt.0.) then
      if (inm>grd%isd) then
        inm=inm-1
      endif
    elseif (xi.gt.1.) then
      if (inm<grd%ied) then
        inm=inm+1
      endif
    endif
    if (yj.lt.0.) then
      if (jnm>grd%jsd) then
        jnm=jnm-1
      endif
    elseif (yj.gt.1.) then
      if (jnm<grd%jed) then
        jnm=jnm+1
      endif
    endif
    lret=pos_within_cell(grd, lon, lat, inm, jnm, xi, yj) ! Update xi and yj
  enddo
  if (abs(inm-i0)>1) then
    write(stderrunit,*) 'pe=',mpp_pe(),'diamonds, adjust: inm,i0,inm-i0=',inm,i0,inm-i0
   !stop 'Moved too far in i without mask!'
  endif
  if (abs(jnm-j0)>1) then
    write(stderrunit,*) 'pe=',mpp_pe(),'diamonds, adjust: jnm,i0,jnm-j0=',jnm,j0,inm-j0
   !stop 'Moved too far in j without mask!'
  endif

  ! Adjust i,j based on xi,yj while bouncing off of masked land cells
  icount=0
  lret=pos_within_cell(grd, lon, lat, i0, j0, xi, yj)
  do while ( .not.lret.and. icount<4 )
    icount=icount+1
    if (xi.lt.0.) then
      if (i>grd%isd) then
        if (grd%msk(i-1,j)>0.) then
          if (i>grd%isd+1) i=i-1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing part from west',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    elseif (xi.gt.1.) then
      if (i<grd%ied) then
        if (grd%msk(i+1,j)>0.) then
          if (i<grd%ied) i=i+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing part from east',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    endif
    if (yj.lt.0.) then
      if (j>grd%jsd) then
        if (grd%msk(i,j-1)>0.) then
          if (j>grd%jsd+1) j=j-1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing part from south',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    elseif (yj.gt.1.) then
      if (j<grd%jed) then
        if (grd%msk(i,j+1)>0.) then
          if (j<grd%jed) j=j+1
        else
         !write(stderr(),'(a,6f8.3,i)') 'diamonds, adjust: bouncing part from north',lon,lat,xi,yj,uvel,vvel,mpp_pe()
          bounced=.true.
        endif
      endif
    endif
    if (bounced) then
      if (xi>1.) xi=1.-posn_eps
      if (xi<0.) xi=posn_eps
      if (yj>1.) yj=1.-posn_eps
      if (yj<0.) yj=posn_eps
      lon=bilin(grd, grd%lon, i, j, xi, yj)
      lat=bilin(grd, grd%lat, i, j, xi, yj)
    endif
    if (debug) then
      if (grd%msk(i,j)==0.) stop 'diamonds, adjust: Part is in land! This should not happen...'
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj
  enddo

 !if (debug) then
 !  if (abs(i-i0)>2) then
 !    stop 'diamonds, adjust: Moved too far in i!'
 !  endif
 !  if (abs(j-j0)>2) then
 !    stop 'diamonds, adjust: Moved too far in j!'
 !  endif
 !endif

  if (.not.bounced.and.lret.and.grd%msk(i,j)>0.) return ! Landed in ocean without bouncing so all is well
  if (.not.bounced.and..not.lret) then ! This implies the part traveled many cells without getting far enough
                                       ! OR that it did not move at all (round-off problem)
    if (debug) then
      write(stderrunit,*) 'diamonds, adjust: lon0, lat0=',lon0,lat0
      write(stderrunit,*) 'diamonds, adjust: xi0, yj0=',xi0,yj0
      write(stderrunit,*) 'diamonds, adjust: i0,j0=',i0,j0 
      write(stderrunit,*) 'diamonds, adjust: lon, lat=',lon,lat
      write(stderrunit,*) 'diamonds, adjust: xi,yj=',xi,yj 
      write(stderrunit,*) 'diamonds, adjust: i,j=',i,j
      write(stderrunit,*) 'diamonds, adjust: inm,jnm=',inm,jnm
      write(stderrunit,*) 'diamonds, adjust: icount=',icount
      lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
      write(stderrunit,*) 'diamonds, adjust: lret=',lret
    endif
    if (abs(i-i0)+abs(j-j0)==0) then
       ! This is a special case due to round off where is_point_in_cell()
       ! returns false but xi and yj are between 0 and 1.
       ! It occurs very rarely but often enough to have brought down
       ! ESM2G four times since the spin-up began. (as of 8/10/2010)
       ! This temporary fix arbitrarily moves the part toward the
       ! center of the current cell.
       xi=(xi-0.5)*(1.-posn_eps)+0.5
       yj=(yj-0.5)*(1.-posn_eps)+0.5
       call error_mesg('diamonds, adjust', 'Part did not move or bounce during iterations AND was not in cell. Adjusting!', WARNING)
    else
      call error_mesg('diamonds, adjust', 'Part iterated many times without bouncing!', WARNING)
    endif
  endif
  if (xi>1.) xi=1.-posn_eps
  if (xi<0.) xi=posn_eps
  if (yj>1.) yj=1.-posn_eps
  if (yj<0.) yj=posn_eps
  lon=bilin(grd, grd%lon, i, j, xi, yj)
  lat=bilin(grd, grd%lat, i, j, xi, yj)
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj

  if (.not. lret) then
    write(0,*) 'i0, j0,=', i0,j0
    write(0,*) 'xi0, yj0,=', xi0,yj0
    write(0,*) 'grd%msk(i0, j0)=', grd%msk(i0,j0)
    write(0,*) 'lon0, lat0,=', lon0,lat0
    write(0,*) 'i,j,lon, lat,grd%msk(i,j)=', i,j,lon,lat,grd%msk(i,j)
    
    write(stderrunit,*) 'diamonds, adjust: Should not get here! Part is not in cell after adjustment'
    if (debug) error=.true.
  endif

end subroutine adjust_index_and_ground

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
