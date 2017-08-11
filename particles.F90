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

use MOM, only : MOM_control_struct

use particles_framework, only: particles_framework_init
use particles_framework, only: particles_gridded, xyt, particle, particles, buffer
use particles_framework, only: verbose, really_debug,debug,use_roundoff_fix
use particles_framework, only: find_cell,find_cell_by_search,count_parts,is_point_in_cell,pos_within_cell
use particles_framework, only: nclasses, budget, old_bug_rotated_weights
use particles_framework, only: bilin,yearday,count_parts,parts_chksum
use particles_framework, only: checksum_gridded,add_new_part_to_list
use particles_framework, only: send_parts_to_other_pes,move_trajectory,move_all_trajectories
use particles_framework, only: record_posn,check_position,print_part,print_parts,print_fld
use particles_framework, only: add_new_part_to_list,delete_particle_from_list,destroy_particle
use particles_framework, only: grd_chksum2,grd_chksum3
use particles_framework, only: offset_part_dates
use particles_framework, only: count_parts_in_list,list_chksum
use particles_framework, only: sum_mass, sum_heat
use particles_framework, only: monitor_a_part,move_part_between_cells, update_halo_particles
use particles_framework, only: is_point_within_xi_yj_bounds

use particles_io,        only: particles_io_init,write_restart,write_trajectory
!use particles_io,        only: read_restart_parts,read_restart_parts_orig,read_restart_calving

implicit none ; private

public particles_init !, particles_end, particles_run, particles_stock_pe, particles
public particles_end, particles_run, particles
public particles_save_restart


real, parameter :: pi_180=pi/180.  !< Converts degrees to radians
real, parameter :: r180_pi=180./pi !< Converts radians to degrees
real, parameter :: Rearth=6360000. !< Radius of earth (m)
real, parameter :: rho_ice=916.7 !< Density of fresh ice @ 0oC (kg/m^3)
real, parameter :: rho_water=999.8 !< Density of fresh water @ 0oC (kg/m^3)
real, parameter :: rho_air=1.1 !< Density of air @ 0oC (kg/m^3) ???
real, parameter :: rho_seawater=1025. !< Approx. density of surface sea water @ 0oC (kg/m^3)
real, parameter :: gravity=9.8 !< Gravitational acceleratio (m/s^2)
real, parameter :: Cd_av=1.3 !< (Vertical) Drag coefficient between parts and atmos (?)
real, parameter :: Cd_ah=0.0055 !< (Horizontal) Drag coefficient between parts and atmos (?)
real, parameter :: Cd_wv=0.9 !< (Vertical) Drag coefficient between parts and ocean (?)
real, parameter :: Cd_wh=0.0012 !< (Horizontal) Drag coefficient between parts and ocean (?)
real, parameter :: Cd_iv=0.9 !< (Vertical) Drag coefficient between parts and sea-ice (?)


#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains

! ##############################################################################
subroutine particles_init(parts, Grid, Time, MOM_CS, dt, axes)

 use particles_io, only: read_restart_parts, particles_io_init

 type(particles), pointer :: parts
 type(ocean_grid_type), pointer :: Grid !< Grid type from parent model
 type(time_type), intent(in) :: Time !< Time type from parent model
 type(MOM_control_struct), pointer, intent(in) :: MOM_CS !< Pointer to a MOM control structure
 real, intent(in)            :: dt !< particle timestep in seconds
 integer, dimension(2), intent(in) :: axes !< diagnostic axis ids

 integer :: io_layout(2)
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


 call mpp_clock_begin(parts%clock_ior)
 call particles_io_init(parts,io_layout)
 call read_restart_parts(parts,Time, MOM_CS)
 call parts_chksum(parts, 'read_restart_particles')
 call mpp_clock_end(parts%clock_ior)

 if (really_debug) call print_parts(stderrunit,parts,'particles_init, initial status')


end subroutine particles_init

! ###############################################################################
!> Returns the axis intercept of a line AB
!!
!! This routine returns the position (x0,y0) at which a line AB intercepts the x or y axis.
!! The value No_intercept_val is returned when the line does not intercept the axis.
subroutine intercept_of_a_line(Ax, Ay, Bx, By, axes1, x0, y0)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  character, intent(in) :: axes1 !< Either 'x' or 'y'
  real, intent(out) :: x0 !< x-position of intercept
  real, intent(out) :: y0 !< y-position of intercept
  ! Local variables
  real :: No_intercept_val ! Huge value used to make sure that the intercept is outside the triangle in the parallel case.

  No_intercept_val=100000000000. ! Huge value used to make sure that the intercept is outside the triangle in the parallel case.
  x0=No_intercept_val
  y0=No_intercept_val

  if (axes1=='x') then ! x intercept
    if (Ay.ne.By) then
      x0=Ax -(((Ax-Bx)/(Ay-By))*Ay)
      y0=0.
    endif
  endif

  if (axes1=='y') then ! y intercept
    if (Ax.ne.Bx) then
      x0=0.
      y0=-(((Ay-By)/(Ax-Bx))*Ax)+Ay
    endif
  endif
end subroutine intercept_of_a_line

! ###############################################################################
!> Calculates the two areas of a triangle divided by an axis line
!!
!! This function calculates the area of a triangle on opposite sides of an axis when the
!! triangle is split with two points on one side, and one point on the other.
!! In this function, A is the point on one side of the axis, and B,C are on the opposite sides.
!! \todo You should change this name a little, so that it not similar the other routine.
subroutine Area_of_triangle_across_axes(Ax, Ay, Bx, By, Cx, Cy, axis1, Area_positive, Area_negative)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  character, intent(in) :: axis1 !< Either 'x' or 'y'
  real, intent(out) :: Area_positive !< Area on negative side of axis line
  real, intent(out) :: Area_negative !< Area on positive side of axis line
  ! Local variables
  real :: pABx, pABy, pACx, pACy, A0
  real :: A_half_triangle, A_triangle

  A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy)

  call intercept_of_a_line(Ax,Ay,Bx,By,axis1,pABx, pABy)
  call intercept_of_a_line(Ax,Ay,Cx,Cy,axis1,pACx, pACy)

  if (axis1=='x')  A0=Ay; !Value used for if statements (deciding up/down vs left/right)
  if (axis1=='y')  A0=Ax; !Value used for if statements (deciding up/down vs left/right)

  A_half_triangle=Area_of_triangle(Ax,Ay,pABx,pABy,pACx,pACy)
  if (A0>=0.) then
    Area_positive= A_half_triangle
    Area_negative= A_triangle-A_half_triangle
  else
    Area_positive= A_triangle-A_half_triangle
    Area_negative= A_half_triangle
  endif

end subroutine Area_of_triangle_across_axes


! ###############################################################################

!> Calculates the area of a triangle on either side of an axis, if any.
!!
!! This routine gives you the area of a triangle on opposite sides of the axis specified.
!! It also takes care of the special case where the triangle is totally on one side.
!! This routine calls Area_of_triangle_across_axes to calculate the areas when the triangles are split.
subroutine divding_triangle_across_axes(Ax, Ay, Bx, By, Cx, Cy, axes1, Area_positive, Area_negative)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  character, intent(in) :: axes1 !< Either 'x' or 'y'
  real, intent(out) :: Area_positive !< Area on negative side of axis line
  real, intent(out) :: Area_negative !< Area on positive side of axis line
  ! Local variables
  real :: A0,B0,C0
  real A_triangle

  if (axes1=='x') then ! Use the y-coordinates for if statements to see which side of the line you are on
    A0=Ay
    B0=By
    C0=Cy
  endif
  if (axes1=='y') then ! Use the y-coordinates for if statements to see which side of the line you are on
    A0=Ax
    B0=Bx
    C0=Cx
  endif

  A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy)
  if ((B0*C0)>0.) then ! B and C are on the same side  (and non-zero)
    if ((A0*B0).ge.0.) then ! all three on the same side (if it equals zero, then A0=0 and the others are not)
      if ((A0>0.)  .or.  ((A0==0.) .and.  (B0>0.))) then
        Area_positive= A_triangle
        Area_negative= 0.
      else
        Area_positive= 0.
        Area_negative= A_triangle
      endif
    else  !A is on the opposite side to B and C
      call Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1,Area_positive, Area_negative)
    endif

  elseif ((B0*C0)<0.) then !B and C are on the opposite sides
    if ((A0*B0).ge. 0.) then !C is all alone
      call Area_of_triangle_across_axes(Cx,Cy,Bx,By,Ax,Ay,axes1,Area_positive, Area_negative)
    else !B is all alone
      call Area_of_triangle_across_axes(Bx,By,Cx,Cy,Ax,Ay,axes1,Area_positive, Area_negative)
    endif

  else  !This is the case when either B or C is equal to zero (or both), A0 could be zero too.
    if (((A0.eq.0.) .and. (B0.eq.0.)) .and. (C0.eq.0.)) then
      Area_positive= 0.
      Area_negative= 0.
    elseif ((A0*B0<0.)  .or.  (A0*C0<0.)) then    !A, B are on opposite sides, and C is zero.  OR  A, C are on opposite sides, and B is zero.
      call Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1,Area_positive, Area_negative)
    elseif (((A0*B0>0.) .or. (A0*C0>0.)) .or. (((abs(A0)>0.) .and. (B0==0.)) .and. (C0==0.))) then
      if (A0>0.) then
        Area_positive= A_triangle
        Area_negative= 0.
      else
        Area_positive= 0.
        Area_negative= A_triangle
      endif

    elseif (A0.eq. 0.) then   !(one of B,C is zero too)
      if ((B0>0.) .or. (C0>0.)) then
        Area_positive= A_triangle
        Area_negative= 0.
      elseif ((B0<0.) .or. (C0<0.)) then
        Area_positive= 0.
        Area_negative= A_triangle
      else
        call error_mesg('diamonds, particle_run', 'Logical error inside triangle dividing routine', FATAL)
      endif
    else
      call error_mesg('diamonds, particle_run', 'Another logical error inside triangle dividing routine', FATAL)
    endif
  endif
end subroutine divding_triangle_across_axes

! ############################################################################################

!> Areas of a triangle divided into quadrants
!!
!! This routine takes a triangle, and finds the intersection with the four quadrants.
subroutine Triangle_divided_into_four_quadrants(Ax, Ay, Bx, By, Cx, Cy, Area_triangle, Area_Q1, Area_Q2 ,Area_Q3 ,Area_Q4)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  real, intent(out) :: Area_triangle !< Are of triangle
  real, intent(out) :: Area_Q1 !< Are in quadrant 1
  real, intent(out) :: Area_Q2 !< Are in quadrant 2
  real, intent(out) :: Area_Q3 !< Are in quadrant 2
  real, intent(out) :: Area_Q4 !< Are in quadrant 4
  ! Local variables
  real :: Area_Upper, Area_Lower, Area_Right, Area_Left
  real :: px, py , qx , qy
  real :: Area_key_quadrant,Error
  real :: tol
  integer :: Key_quadrant
  integer ::sig_fig
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()
  tol=1.e-10

  Area_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy)

  ! Calculating area across axes
  call divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'x',Area_Upper ,Area_Lower)
  call divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'y',Area_Right ,Area_Left)

  ! Decide if the origin is in the triangle. If so, then you have to divide the area 4 ways
  ! This is done by finding a quadrant where the intersection between the triangle and quadrant forms a new triangle
  ! (This occurs when on of the sides of the triangle  intersects both the x and y axis)
  if (point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.)) then
    ! Find a line in the triangle that cuts both axes in/on the triangle
    call intercept_of_a_line(Ax,Ay,Bx,By,'x',px,py); !x_intercept
    call intercept_of_a_line(Ax,Ay,Bx,By,'y',qx,qy); !y_intercept
    ! Note that the 1. here means that we include points on the boundary of the triangle.
    if (.not.((point_in_interval(Ax,Ay,Bx,By,px,py)) .and. (point_in_interval(Ax,Ay,Bx,By,qx,qy)))) then
      call intercept_of_a_line(Ax,Ay,Cx,Cy,'x',px,py); !x_intercept
      call intercept_of_a_line(Ax,Ay,Cx,Cy,'y',qx,qy); !y_intercept
      if (.not.((point_in_interval(Ax,Ay,Cx,Cy,px,py)) .and. (point_in_interval(Ax,Ay,Cx,Cy,qx,qy)))) then
        call intercept_of_a_line(Bx,By,Cx,Cy,'x',px,py); !x_intercept
        call intercept_of_a_line(Bx,By,Cx,Cy,'y',qx,qy); !y_intercept
        if (.not.((point_in_interval(Bx,By,Cx,Cy,px,py)) .and. (point_in_interval(Bx,By,Cx,Cy,qx,qy)))) then
          ! You should not get here, but there might be some bugs in the code to do with points exactly falling on axes.
          !if (mpp_pe().eq.12) then
            write(stderrunit,*) 'diamonds,corners', Ax,Ay,Bx,By,Cx,Cy
          !endif
          call error_mesg('diamonds, particle_run', 'Something went wrong with Triangle_divide_into_four_quadrants', FATAL)
        endif
      endif
    endif

    ! Assigning quadrants. Key_quadrant is the quadrant with the baby triangle in it.
    Area_key_quadrant=Area_of_triangle(px,py,qx,qy,0.,0.)
    if ((px.ge. 0.) .and. (qy.ge. 0.)) then  !First quadrant
      Key_quadrant=1
    elseif ((px.lt.0.) .and. (qy.ge. 0.)) then  !Second quadrant
      Key_quadrant=2
    elseif ((px.lt. 0.) .and. (qy.lt. 0.)) then !Third quadrant
      Key_quadrant=3
    elseif ((px.ge. 0.) .and. (qy.lt. 0.)) then !Forth quadrant
      Key_quadrant=4
    else  !
      call error_mesg('diamonds, particle_run', 'None of the quadrants are Key', WARNING)
      write(stderrunit,*) 'diamonds, Triangle, px,qy', px,qy
    endif

  else ! At least one quadrant is empty, and this can be used to find the areas in the other quadrant.  Assigning quadrants. Key_quadrant is the empty quadrant.
    Area_key_quadrant=0
    if      ( (.not. ((((Ax>0.) .and. (Ay>0.)) .or. ((Bx>0.) .and. (By> 0.))) .or. ((Cx>0.) .and. (Cy> 0.)))) .and. ((Area_Upper+Area_Right).le.Area_triangle) ) then
      ! No points land in this quadrant and triangle does not cross the quadrant
      Key_quadrant=1
    elseif  ( (.not. ((((Ax<0.) .and. (Ay>0)) .or. ((Bx<0.) .and. (By>0.))) .or. ((Cx<0.) .and. (Cy>0.)))) .and. ((Area_Upper+Area_Left).le. Area_triangle) ) then
      Key_quadrant=2
    elseif  ( (.not. ((((Ax<0.) .and. (Ay<0.)) .or. ((Bx<0.) .and. (By< 0.))) .or. ((Cx<0.) .and. (Cy< 0.)))) .and. ((Area_Lower+Area_Left) .le.Area_triangle) ) then
      Key_quadrant=3
    else
      Key_quadrant=4
    endif
  endif

  ! Assign values to quadrants
  if (Key_quadrant .eq. 1) then
    Area_Q1=Area_key_quadrant
    Area_Q2=Area_Upper-Area_Q1
    Area_Q4=Area_Right-Area_Q1
    !Area_Q3=Area_Left-Area_Q2 ! These lines have been changes so that the sum of the 4 quadrants exactly matches the triangle area.
    Area_Q3=Area_triangle-(Area_Q1+Area_Q2+Area_Q4)
  elseif (Key_quadrant .eq. 2) then
    Area_Q2=Area_key_quadrant
    Area_Q1=Area_Upper-Area_Q2
    Area_Q4=Area_Right-Area_Q1
    !Area_Q3=Area_Left-Area_Q2
    Area_Q3=Area_triangle-(Area_Q1+Area_Q2+Area_Q4)
  elseif (Key_quadrant==3) then
    Area_Q3=Area_key_quadrant
    Area_Q2=Area_Left-Area_Q3
    Area_Q1=Area_Upper-Area_Q2
    !Area_Q4=Area_Right-Area_Q1
    Area_Q4=Area_triangle-(Area_Q1+Area_Q2+Area_Q3)
  elseif (Key_quadrant==4) then
    Area_Q4=Area_key_quadrant
    Area_Q1=Area_Right-Area_Q4
    Area_Q2=Area_Upper-Area_Q1
    !Area_Q3=Area_Left-Area_Q2
    Area_Q3=Area_triangle-(Area_Q1+Area_Q2+Area_Q4)
  else
    call error_mesg('diamonds, particle_run', 'Logical error inside triangle into four quadrants. Should not get here.', FATAL)
  endif

  Area_Q1=max(Area_Q1,0.)
  Area_Q2=max(Area_Q2,0.)
  Area_Q3=max(Area_Q3,0.)
  Area_Q4=max(Area_Q4,0.)


  Error=abs(Area_Q1+Area_Q2+Area_Q3+Area_Q4-Area_triangle)
  if (Error>tol) then
    call error_mesg('diamonds, triangle spreading', 'Triangle not evaluated accurately!!', WARNING)
    !if (mpp_pe().eq.mpp_root_pe()) then
    if (mpp_pe().eq. 20) then
      write(stderrunit,*) 'diamonds, Triangle corners:',Ax,Ay,Bx,By,Cx,Cy
      write(stderrunit,*) 'diamonds, Triangle, Full Area', Area_Q1+ Area_Q2+ Area_Q3+ Area_Q4
      write(stderrunit,*) 'diamonds, Triangle, Areas', Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
      write(stderrunit,*) 'diamonds, Triangle, Areas', Error
      write(stderrunit,*) 'diamonds, Key quadrant',Key_quadrant,Area_key_quadrant
      write(stderrunit,*) 'diamonds, point in triangle',(point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.))
      write(stderrunit,*) 'diamonds, halves',Area_Upper,Area_Lower,Area_Right,Area_Left
    endif
  endif

end subroutine Triangle_divided_into_four_quadrants


! ###############################################################################

!> Rotates a point clockwise about origin and then translates by x0,y0
subroutine rotate_and_translate(px, py, theta, x0, y0)
  ! Arguments
  real, intent(in) :: x0 !< x-direction shift
  real, intent(in) :: y0 !< y-direction shift
  real, intent(in) :: theta !< Angle to rotate (degrees)
  real, intent(inout) :: px !< x-coordinate of point
  real, intent(inout) :: py !< y-coordinate of point
  ! Local variables
  real :: px_temp,py_temp

  ! Rotation
  px_temp = ( cos(theta*pi/180)*px) + (sin(theta*pi/180)*py)
  py_temp = (-sin(theta*pi/180)*px) + (cos(theta*pi/180)*py)

  ! Translation
  px= px_temp + x0
  py= py_temp + y0
end subroutine rotate_and_translate

! ##############################################################################

subroutine Hexagon_into_quadrants_using_triangles(x0, y0, H, theta, Area_hex ,Area_Q1, Area_Q2, Area_Q3, Area_Q4)
  ! Arguments
  real, intent(in) :: x0 !< x-coordinate of center of hexagon
  real, intent(in) :: y0 !< y-coordinate of center of hexagon
  real, intent(in) :: H !< Apothem (inner radius of hexagon)
  real, intent(in) :: theta !< Orientation angle of hexagon
  real, intent(out) :: Area_hex !< Area of hexagon
  real, intent(out) :: Area_Q1 !< Are in quadrant 1
  real, intent(out) :: Area_Q2 !< Are in quadrant 2
  real, intent(out) :: Area_Q3 !< Are in quadrant 2
  real, intent(out) :: Area_Q4 !< Are in quadrant 4
  ! Local variables
  real :: C1x, C2x, C3x, C4x, C5x, C6x
  real :: C1y, C2y, C3y, C4y, C5y, C6y
  real :: T12_Area, T12_Q1, T12_Q2, T12_Q3, T12_Q4
  real :: T23_Area, T23_Q1, T23_Q2, T23_Q3, T23_Q4
  real :: T34_Area, T34_Q1, T34_Q2, T34_Q3, T34_Q4
  real :: T45_Area, T45_Q1, T45_Q2, T45_Q3, T45_Q4
  real :: T56_Area, T56_Q1, T56_Q2, T56_Q3, T56_Q4
  real :: T61_Area, T61_Q1, T61_Q2, T61_Q3, T61_Q4
  real :: S, exact_hex_area, Error
  real :: tol
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()
  tol=1.e-10

  ! Length of side of Hexagon
  S=(2/sqrt(3.))*H

  ! Finding positions of corners
  C1x=S           ; C1y=0.  !Corner 1 (right)
  C2x=H/sqrt(3.)  ; C2y=H;  !Corner 2 (top right)
  C3x=-H/sqrt(3.) ; C3y=H;  !Corner 3 (top left)
  C4x=-S          ; C4y=0.; !Corner 4 (left)
  C5x=-H/sqrt(3.) ; C5y=-H; !Corner 5 (bottom left)
  C6x=H/sqrt(3.)  ; C6y=-H; !Corner 6 (bottom right)

  ! Finding positions of corners
  call rotate_and_translate(C1x,C1y,theta,x0,y0)
  call rotate_and_translate(C2x,C2y,theta,x0,y0)
  call rotate_and_translate(C3x,C3y,theta,x0,y0)
  call rotate_and_translate(C4x,C4y,theta,x0,y0)
  call rotate_and_translate(C5x,C5y,theta,x0,y0)
  call rotate_and_translate(C6x,C6y,theta,x0,y0)

  ! Area of Hexagon is the sum of the triangles
  call Triangle_divided_into_four_quadrants(x0,y0,C1x,C1y,C2x,C2y,T12_Area,T12_Q1,T12_Q2,T12_Q3,T12_Q4); !Triangle 012
  call Triangle_divided_into_four_quadrants(x0,y0,C2x,C2y,C3x,C3y,T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4); !Triangle 023
  call Triangle_divided_into_four_quadrants(x0,y0,C3x,C3y,C4x,C4y,T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4); !Triangle 034
  call Triangle_divided_into_four_quadrants(x0,y0,C4x,C4y,C5x,C5y,T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4); !Triangle 045
  call Triangle_divided_into_four_quadrants(x0,y0,C5x,C5y,C6x,C6y,T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4); !Triangle 056
  call Triangle_divided_into_four_quadrants(x0,y0,C6x,C6y,C1x,C1y,T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4); !Triangle 061

  ! Summing up the triangles
  Area_hex=T12_Area+T23_Area+T34_Area+T45_Area+T56_Area+T61_Area
  Area_Q1=T12_Q1+T23_Q1+T34_Q1+T45_Q1+T56_Q1+T61_Q1
  Area_Q2=T12_Q2+T23_Q2+T34_Q2+T45_Q2+T56_Q2+T61_Q2
  Area_Q3=T12_Q3+T23_Q3+T34_Q3+T45_Q3+T56_Q3+T61_Q3
  Area_Q4=T12_Q4+T23_Q4+T34_Q4+T45_Q4+T56_Q4+T61_Q4

  Area_Q1=max(Area_Q1,0.)
  Area_Q2=max(Area_Q2,0.)
  Area_Q3=max(Area_Q3,0.)
  Area_Q4=max(Area_Q4,0.)

  Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
  if ((abs(Error)>tol))then
    if (mpp_pe().eq.mpp_root_pe()) then
      call error_mesg('diamonds, hexagonal spreading', 'Hexagon error is large!!', WARNING)
      write(stderrunit,*) 'diamonds, hex error, H,x0,y0, Error', H, x0 , y0, Error
      write(stderrunit,*) 'diamonds, hex error, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
      write(stderrunit,*) 'diamonds, Triangle1',C1x,C1y,C2x,C2y,T12_Area,T12_Q1,T12_Q2,T12_Q3,T12_Q4,(T12_Q1+T12_Q2+T12_Q3+T12_Q4-T12_Area)
      write(stderrunit,*) 'diamonds, Triangle2',C2x,C2y,C3x,C3y,T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4,(T23_Q1+T23_Q2+T23_Q3+T23_Q4-T23_Area)
      write(stderrunit,*) 'diamonds, Triangle3',C3x,C3y,C4x,C4y,T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4,(T34_Q1+T34_Q2+T34_Q3+T34_Q4-T34_Area)
      write(stderrunit,*) 'diamonds, Triangle4',C4x,C4y,C5x,C5y,T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4,(T45_Q1+T45_Q2+T45_Q3+T45_Q4-T45_Area)
      write(stderrunit,*) 'diamonds, Triangle5',C5x,C5y,C6x,C6y,T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4,(T56_Q1+T56_Q2+T56_Q3+T56_Q4-T56_Area)
      write(stderrunit,*) 'diamonds, Triangle6',C6x,C6y,C1x,C1y,T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4,(T61_Q1+T61_Q2+T61_Q3+T61_Q4-T61_Area)
    endif
  endif

  exact_hex_area=((3.*sqrt(3.)/2)*(S*S))
  if (abs(Area_hex-exact_hex_area)>tol) then
    call error_mesg('diamonds, hexagonal spreading', 'Hexagon not evaluated accurately!!', WARNING)
    if (mpp_pe().eq.mpp_root_pe()) then
      write(stderrunit,*) 'diamonds, hex calculations, H,x0,y0', H, x0 , y0
      write(stderrunit,*) 'diamonds, hex calculations, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
    endif
  endif

  ! Adjust Areas so that the error is zero by subtracting the error from the largest sector.
   if  (((Area_Q1>=Area_Q2) .and. (Area_Q1>=Area_Q3)) .and. (Area_Q1>=Area_Q4)) then
     Area_Q1=Area_Q1+Error
   elseif  (((Area_Q2>=Area_Q1) .and. (Area_Q2>=Area_Q3)) .and. (Area_Q2>=Area_Q4)) then
     Area_Q2=Area_Q2+Error
   elseif  (((Area_Q3>=Area_Q1) .and. (Area_Q3>=Area_Q2)) .and. (Area_Q3>=Area_Q4)) then
     Area_Q3=Area_Q3+Error
   elseif  (((Area_Q4>=Area_Q1) .and. (Area_Q4>=Area_Q2)) .and. (Area_Q4>=Area_Q3)) then
     Area_Q4=Area_Q4+Error
   else
     call error_mesg('diamonds, hexagonal spreading', 'Error in hexagon is larger than any quadrant!!', WARNING)
     if (mpp_pe().eq.mpp_root_pe()) then
      write(stderrunit,*) 'diamonds, hex quadrants, H,x0,y0', H, x0 , y0, Error
      write(stderrunit,*) 'diamonds, hex quadrants, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
     endif
   endif

end subroutine Hexagon_into_quadrants_using_triangles

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

! ##############################################################################


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

!> The main driver the steps updates particles
subroutine particles_run(parts, time, uo, vo, stagger)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  type(time_type), intent(in) :: time !< Model time
  real, dimension(:,:),intent(in) :: uo !< Ocean zonal velocity (m/s)
  real, dimension(:,:),intent(in) :: vo !< Ocean meridional velocity (m/s)
  integer, optional, intent(in) :: stagger

  ! Local variables
  integer :: iyr, imon, iday, ihr, imin, isec, k
  type(particles_gridded), pointer :: grd
  logical :: lerr, sample_traj, write_traj, lverbose
  real :: grdd_u_particle, grdd_v_particle
  integer :: i, j, Iu, ju, iv, Jv, Iu_off, ju_off, iv_off, Jv_off
  real :: mask
  real, dimension(:,:), allocatable :: uC_tmp, vC_tmp, uA_tmp, vA_tmp
  integer :: vel_stagger
  real, dimension(:,:), allocatable :: iCount
  integer :: stderrunit


  print *,' In particles_run 0001'

  ! Get the stderr unit number
  stderrunit = stderr()

  vel_stagger = CGRID_NE ; if (present(stagger)) vel_stagger = stagger

  ! For convenience
  grd=>parts%grd
  grd%u_particle(:,:)=0.
  grd%v_particle(:,:)=0.

  !Initializing _on_ocean_fields
  grd%Uvel_on_ocean(:,:,:)=0. ;   grd%Vvel_on_ocean(:,:,:)=0.

  print *,' In particles_run 0002'

  ! Manage time
  call get_date(time, iyr, imon, iday, ihr, imin, isec)
  parts%current_year=iyr
  parts%current_yearday=yearday(imon, iday, ihr, imin, isec)
  ! Turn on sampling of trajectories, verbosity, budgets
  sample_traj=.false.
  if ( (parts%traj_sample_hrs>0)  .and. (.not. parts%ignore_traj) ) then
    if (mod(60*60*24*iday+ 60*60*ihr + 60*imin + isec ,60*60*parts%traj_sample_hrs).eq.0) &
        sample_traj=.true.
  endif
  write_traj=.false.
  if ((parts%traj_write_hrs>0) .and. (.not. parts%ignore_traj))  then
     if (mod(60*60*24*iday+ 60*60*ihr + 60*imin + isec ,60*60*parts%traj_write_hrs).eq.0) &
         write_traj=.true.
  endif
  lverbose=.false.
  if (parts%verbose_hrs>0) then
     if (mod(24*iday+ihr+(imin/60.),float(parts%verbose_hrs)).eq.0) lverbose=verbose
  endif

  if (mpp_pe()==mpp_root_pe().and.lverbose) write(*,'(a,3i5,a,3i5,a,i5,f8.3)') &
       'diamonds: y,m,d=',iyr, imon, iday,' h,m,s=', ihr, imin, isec, &
       ' yr,yrdy=', parts%current_year, parts%current_yearday

  print *,' In particles_run 0003'

 !call sanitize_field(grd%calving,1.e20)

  if (vel_stagger == BGRID_NE) then
    ! Copy ocean and ice velocities. They are already on B-grid u-points.
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
      ! Interpolate ocean and ice velocities from C-grid velocity points.
      Iu = i + Iu_off ; ju = j + ju_off ; iv = i + iv_off ; Jv = j + Jv_off
      ! This masking is needed for now to prevent particles from running up on to land.
      mask = min(grd%msk(i,j), grd%msk(i+1,j), grd%msk(i,j+1), grd%msk(i+1,j+1))
      grd%uo(I,J) = mask * 0.5*(uo(Iu,ju)+uo(Iu,ju+1))
      grd%vo(I,J) = mask * 0.5*(vo(iv,Jv)+vo(iv+1,Jv))
    enddo ; enddo
  else
    call error_mesg('diamonds, particle_run', 'Unrecognized value of stagger!', FATAL)
  endif

  print *,' In particles_run 0004'
  call mpp_update_domains(grd%uo, grd%vo, grd%domain, gridtype=BGRID_NE)

  ! Make sure that gridded values agree with mask  (to get ride of NaN values)
  do i=grd%isd,grd%ied ; do j=grd%jsd,grd%jed
    ! Initializing all gridded values to zero
    if (grd%msk(i,j).lt. 0.5) then
      grd%uo(i,j) = 0.0 ;  grd%vo(i,j) = 0.0
    endif
    if (grd%uo(i,j) .ne. grd%uo(i,j)) grd%uo(i,j)=0.
    if (grd%vo(i,j) .ne. grd%vo(i,j)) grd%vo(i,j)=0.
  enddo; enddo

  if (debug) call parts_chksum(parts, 'run parts (top)')
  if (debug) call checksum_gridded(parts%grd, 'top of s/r run')


  if (.not.parts%Static_particles) then
    call evolve_particles(parts)
    if (parts%debug_particle_with_id>0) call monitor_a_part(parts, 'particles_run, after evolve()     ')
  endif
  call move_part_between_cells(parts)  !Markpoint6
  if (parts%debug_particle_with_id>0) call monitor_a_part(parts, 'particles_run, after move_lists() ')
  if (debug) call parts_chksum(parts, 'run parts (evolved)',ignore_halo_violation=.true.)
  if (debug) call checksum_gridded(parts%grd, 's/r run after evolve')

  call send_parts_to_other_pes(parts)
  if (parts%debug_particle_with_id>0) call monitor_a_part(parts, 'particles_run, after send_parts() ')

  !Creating gridded fields from new particles
  call create_gridded_particles_fields(parts)

  ! For each part, record
  if (sample_traj) call record_posn(parts)
  if (write_traj) then
    call move_all_trajectories(parts)
    call write_trajectory(parts%trajectories, parts%save_short_traj)
  endif

  ! Gridded diagnostics
  if (grd%id_uo>0) &
    lerr=send_data(grd%id_uo, grd%uo(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_vo>0) &
    lerr=send_data(grd%id_vo, grd%vo(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_u_particle>0) &
    lerr=send_data(grd%id_u_particle, grd%u_particle(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_v_particle>0) &
    lerr=send_data(grd%id_v_particle, grd%v_particle(grd%isc:grd%iec,grd%jsc:grd%jec), Time)
  if (grd%id_count>0) then
    allocate( iCount(grd%isc:grd%iec,grd%jsc:grd%jec) ); iCount(:,:)=0
    do j = grd%jsc, grd%jec ; do i = grd%isc, grd%iec
      iCount(i,j) = count_parts_in_list(parts%list(i,j)%first)
    enddo ; enddo
    lerr=send_data(grd%id_count, iCount(:,:), Time)
    deallocate( iCount )
  endif
  if (grd%id_chksum>0) then
    allocate( iCount(grd%isc:grd%iec,grd%jsc:grd%jec) ); iCount(:,:)=0
    do j = grd%jsc, grd%jec ; do i = grd%isc, grd%iec
      iCount(i,j) = list_chksum(parts%list(i,j)%first)
    enddo ; enddo
    lerr=send_data(grd%id_chksum, iCount(:,:), Time)
    deallocate( iCount )
  endif

  ! Dump particles to screen
  if (really_debug) call print_parts(stderrunit,parts,'particles_run, status')

  if (debug) call parts_chksum(parts, 'run parts (bot)')
  if (debug) call checksum_gridded(parts%grd, 'end of s/r run')
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

  ! Get the stderr unit number
  stderrunit = stderr()

  ! For convenience
  grd=>parts%grd

  interactive_particles_on=parts%interactive_particles_on
  Runge_not_Verlet=parts%Runge_not_Verlet

  do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
    part=>parts%list(grdi,grdj)%first
    do while (associated(part)) ! loop over all parts
      !if (part%static_part .lt. 0.5) then  !Only allow non-static particles to evolve !LUYU: not needed, because all the particles are non-static.

        !Checking it everything is ok:
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

	! Interpolate gridded velocity fields to part and generate uvel and vvel
	call interp_flds(grd,part%ine,part%jne,part%xi,part%yj,part%uvel, part%vvel) 

          !Time stepping schemes:
          if (Runge_not_Verlet) then
            call Runge_Kutta_stepping(parts,part, axn, ayn, bxn, byn, uveln, vveln,lonn, latn, i, j, xi, yj)
          endif
          if (.not.Runge_not_Verlet) then 
            !call verlet_stepping(parts,part, axn, ayn, bxn, byn, uveln, vveln) ! LUYU: consider RK for now
          endif

          !Used for testing the ocean response to fixed particle motion.
          if (parts%override_particle_velocities) then
            uveln  = parts%u_override
            vveln  = parts%v_override
          endif

        ! Saving all the particle variables.
        part%axn=axn
        part%ayn=ayn
        part%bxn=bxn
        part%byn=byn
        part%uvel=uveln
        part%vvel=vveln

        if (Runge_not_Verlet) then
          part%lon=lonn  ;   part%lat=latn
          part%ine=i     ;   part%jne=j
          part%xi=xi     ;   part%yj=yj
        else
          !if (.not. interactive_particles_on) call update_verlet_position(parts,part)
        endif

        !call interp_flds(grd, i, j, xi, yj, part%uo, part%vo, part%ui, part%vi, part%ua, part%va, part%ssh_x, part%ssh_y, part%sst)
        !if (debug) call print_part(stderr(), part, 'evolve_particle, final posn.')
        if (debug) call check_position(grd, part, 'evolve_particle (bot)')
      !endif
      part=>part%next
    enddo ! loop over all parts
  enddo ; enddo

  ! When we are using interactive particles, we update the (old) particle positions and velocities in a second loop, all together (to make code order invarient)
  if (interactive_particles_on) then
    do grdj = grd%jsc,grd%jec ; do grdi = grd%isc,grd%iec
      part=>parts%list(grdi,grdj)%first
      do while (associated(part)) ! loop over all parts
        if (part%static_part .lt. 0.5) then  !Only allow non-static particles to evolve
         !if (.not. Runge_not_Verlet) call update_verlet_position(parts,part)

         !Updating old velocities (for use in particle interactions)
          part%uvel_old=part%uvel
          part%vvel_old=part%vvel
          part%lon_old=part%lon  ! lon_old, lat_old are not really needed for Verlet. But are needed for RK
          part%lat_old=part%lat
        endif
        part=>part%next
      enddo ! loop over all parts
    enddo ; enddo
  endif

end subroutine evolve_particles

!> Calculate explicit and implicit accelerations, new velocity, and new position, using the fourth order Runge-Kutta  method
subroutine Runge_Kutta_stepping(parts, part, axn, ayn, bxn, byn, uveln, vveln, lonn, latn, i, j, xi, yj)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  type(particle), pointer, intent(inout) :: part !< particle
  real, intent(out) :: axn !< Explicit zonal acceleration (m/s2)
  real, intent(out) :: ayn !< Explicit meridional acceleration (m/s2)
  real, intent(out) :: bxn !< Implicit zonal acceleration (m/s2)
  real, intent(out) :: byn !< Implicit meridional acceleration (m/s2)
  real, intent(out) :: uveln !< New zonal velocity (m/s)
  real, intent(out) :: vveln !< New meridional velocity (m/s)
  real, intent(out) :: lonn !< New longitude (degree E)
  real, intent(out) :: latn !< New latitude (degree N)
  integer, intent(out) :: i !< New i-index of containing cell
  integer, intent(out) :: j !< New i-index of containing cell
  real, intent(out) :: xi !< New non-dimensional x-position
  real, intent(out) :: yj !< New non-dimensional y-position
  ! Local variables
  type(particles_gridded), pointer :: grd
  real :: uvel1, vvel1, lon1, lat1, u1, v1, dxdl1, ax1, ay1, axn1, ayn1
  real :: uvel2, vvel2, lon2, lat2, u2, v2, dxdl2, ax2, ay2, axn2, ayn2
  real :: uvel3, vvel3, lon3, lat3, u3, v3, dxdl3, ax3, ay3, axn3, ayn3
  real :: uvel4, vvel4, lon4, lat4, u4, v4, dxdl4, ax4, ay4, axn4, ayn4
  real :: x1, xdot1, xddot1, y1, ydot1, yddot1, xddot1n, yddot1n
  real :: x2, xdot2, xddot2, y2, ydot2, yddot2, xddot2n, yddot2n
  real :: x3, xdot3, xddot3, y3, ydot3, yddot3, xddot3n, yddot3n
  real :: x4, xdot4, xddot4, y4, ydot4, yddot4, xddot4n, yddot4n
  real :: xn, xdotn, xddotn, yn, ydotn, yddotn, xddotnn, yddotnn
  real :: dt, dt_2, dt_6, dydl
  integer :: i1,j1,i2,j2,i3,j3,i4,j4
  integer :: stderrunit
  logical :: bounced, on_tangential_plane, error_flag
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
  dt=parts%dt
  dt_2=0.5*dt
  dt_6=dt/6.

  i=part%ine
  j=part%jne
  xi=part%xi
  yj=part%yj
  bounced=.false.
  on_tangential_plane=.false.
  if ((part%lat>89.) .and. (parts%grd%grid_is_latlon)) on_tangential_plane=.true.
  i1=i;j1=j
  if (parts%add_weight_to_ocean .and. parts%time_average_weight) &
    call spread_mass_across_ocean_cells(parts, part, i, j, xi, yj, part%mass, part%mass_of_bits, 0.25*part%mass_scaling,part%length*part%width, part%thickness)

  ! Loading past accelerations - Alon
  axn=part%axn; ayn=part%ayn !Alon
  axn1=axn; axn2=axn; axn3=axn; axn4=axn
  ayn1=ayn; ayn2=ayn; ayn3=ayn; ayn4=ayn

  ! A1 = A(X1)
  lon1=part%lon; lat1=part%lat
  if (on_tangential_plane) call rotpos_to_tang(lon1,lat1,x1,y1)

  call  convert_from_meters_to_grid(lat1,parts%grd%grid_is_latlon ,dxdl1,dydl)
  !dxdl1=r180_pi/(Rearth*cos(lat1*pi_180))
  !dydl=r180_pi/Rearth
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
  call adjust_index_and_ground(grd, lon2, lat2, uvel2, vvel2, i, j, xi, yj, bounced, error_flag, part%id)
  i2=i; j2=j
  if (parts%add_weight_to_ocean .and. parts%time_average_weight) &
    call spread_mass_across_ocean_cells(parts, part, i, j, xi, yj, part%mass, part%mass_of_bits, 0.25*part%mass_scaling,part%length*part%width, part%thickness)
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon2,lat2,x2,y2)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(parts%grd, lon2, lat2, i, j)) error_flag=.true.
  endif
  if (error_flag) then
    call print_fld(grd, grd%msk, 'msk')
    call print_fld(grd, grd%ssh, 'ssh')
    call print_fld(grd, grd%sst, 'sst')
    call print_fld(grd, grd%sss, 'sss')
    call print_fld(grd, grd%hi, 'hi')
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
    bounced=is_point_in_cell(parts%grd, lon2, lat2, i, j,explain=.true.)
    call error_mesg('diamonds, evolve_particle','part is out of posn at 2!',FATAL)
  endif
  call  convert_from_meters_to_grid(lat2,parts%grd%grid_is_latlon ,dxdl2,dydl)
  !dxdl2=r180_pi/(Rearth*cos(lat2*pi_180))
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
  call adjust_index_and_ground(grd, lon3, lat3, uvel3, vvel3, i, j, xi, yj, bounced, error_flag, part%id)
  i3=i; j3=j
  if (parts%add_weight_to_ocean .and. parts%time_average_weight) &
    call spread_mass_across_ocean_cells(parts, part, i, j, xi, yj, part%mass, part%mass_of_bits, 0.25*part%mass_scaling,part%length*part%width, part%thickness)
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon3,lat3,x3,y3)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(parts%grd, lon3, lat3, i, j)) error_flag=.true.
  endif
  if (error_flag) then
    call print_fld(grd, grd%msk, 'msk')
    call print_fld(grd, grd%ssh, 'ssh')
    call print_fld(grd, grd%sst, 'sst')
    call print_fld(grd, grd%sss, 'sss')
    call print_fld(grd, grd%hi, 'hi')
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
    bounced=is_point_in_cell(parts%grd, lon2, lat2, i, j,explain=.true.)
    call error_mesg('diamonds, evolve_particle','part is out of posn at 3!',FATAL)
  endif
  call  convert_from_meters_to_grid(lat3,parts%grd%grid_is_latlon ,dxdl3,dydl)
  !dxdl3=r180_pi/(Rearth*cos(lat3*pi_180))
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
  call adjust_index_and_ground(grd, lon4, lat4, uvel4, vvel4, i, j, xi, yj, bounced, error_flag, part%id)
  i4=i; j4=j
  ! if (bounced.and.on_tangential_plane) call rotpos_to_tang(lon4,lat4,x4,y4)
  if (.not.error_flag) then
    if (debug .and. .not. is_point_in_cell(parts%grd, lon4, lat4, i, j)) error_flag=.true.
  endif
  if (error_flag) then
    call print_fld(grd, grd%msk, 'msk')
    call print_fld(grd, grd%ssh, 'ssh')
    call print_fld(grd, grd%sst, 'sst')
    call print_fld(grd, grd%sss, 'sss')
    call print_fld(grd, grd%hi, 'hi')
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
  call  convert_from_meters_to_grid(lat4,parts%grd%grid_is_latlon ,dxdl4,dydl)
  !dxdl4=r180_pi/(Rearth*cos(lat4*pi_180))
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
  call adjust_index_and_ground(grd, lonn, latn, uveln, vveln, i, j, xi, yj, bounced, error_flag, part%id)
  if (parts%add_weight_to_ocean .and. parts%time_average_weight) &
    call spread_mass_across_ocean_cells(parts, part, i, j, xi, yj, part%mass, part%mass_of_bits, 0.25*part%mass_scaling,part%length*part%width, part%thickness)

  if (.not.error_flag) then
    if (.not. is_point_in_cell(parts%grd, lonn, latn, i, j)) error_flag=.true.
  endif
  if (error_flag) then
    call print_fld(grd, grd%msk, 'msk')
    call print_fld(grd, grd%ssh, 'ssh')
    call print_fld(grd, grd%sst, 'sst')
    call print_fld(grd, grd%sss, 'sss')
    call print_fld(grd, grd%hi, 'hi')
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
end subroutine Runge_Kutta_stepping

! ###############################################################################
!> Calculate longitude-latitude from tangent plane coordinates
subroutine rotpos_to_tang(lon, lat, x, y, id_in)
  ! Arguments
  real, intent(in) :: lon !< Longitude (degree E)
  real, intent(in) :: lat !< Latitude (degree N)
  real, intent(out) :: x !< x-coordinate in tangent plane
  real, intent(out) :: y !< y-coordinate in tangent plane
  integer(kind=8), intent(in), optional :: id_in !< part identifier
  ! Local variables
  real :: r,colat,clon,slon
  integer :: stderrunit, id

  stderrunit = stderr()
  id=0
  if (present(id_in)) then
        id=id_in
  endif

  if (lat>90.) then
      write(stderrunit,*) 'diamonds, rotpos_to_tang: lat>90 already!',lat, lon, id
      call error_mesg('diamonds, rotpos_to_tang','Something went very wrong!',FATAL)
  endif
  if (lat==90.) then
    write(stderrunit,*) 'diamonds, rotpos_to_tang: lat==90 already!',lat, lon
    call error_mesg('diamonds, rotpos_to_tang','Something went wrong!',FATAL)
  endif

  colat=90.-lat
  r=Rearth*(colat*pi_180)
  clon=cos(lon*pi_180)
  slon=sin(lon*pi_180)
  x=r*clon
  y=r*slon

end subroutine rotpos_to_tang

! ###################################################################################
!> Calculate longitude-latitude from tangent plane coordinates
subroutine rotpos_from_tang(x, y, lon, lat)
  ! Arguments
  real, intent(in) :: x !< x-coordinate in tangent plane
  real, intent(in) :: y !< y-coordinate in tangent plane
  real, intent(out) :: lon !< Longitude (degree E)
  real, intent(out) :: lat !< Latitude (degree N)
  ! Local variables
  real :: r

  r=sqrt(x**2+y**2)
  lat=90.-(r180_pi*r/Rearth)
  lon=r180_pi*acos(x/r)*sign(1.,y)

end subroutine rotpos_from_tang

! ###################################################################################
!> Calculate velocity oriented in geographic coordinates from tangent plane velocity
subroutine rotvec_from_tang(lon, xdot, ydot, uvel, vvel)
  ! Arguments
  real, intent(in) :: lon !< Longitude (degree E)
  real, intent(in) :: xdot !< x-component of velocity in tangent plane (m/s)
  real, intent(in) :: ydot !< y-component of velocity in tangent plane (m/s)
  real, intent(out) :: uvel !< Zonal velocity (m/s)
  real, intent(out) :: vvel !< Meridional velocity (m/s)
  ! Local variables
  real :: clon,slon

  clon=cos(lon*pi_180)
  slon=sin(lon*pi_180)
  uvel=-slon*xdot+clon*ydot
  vvel=-clon*xdot-slon*ydot

end subroutine rotvec_from_tang

! ###################################################################################
!> Calculates tangent plane velocity from velocity in velocity oriented in geographic coordinates
subroutine rotvec_to_tang(lon, uvel, vvel, xdot, ydot)
  ! Arguments
  real, intent(in) :: lon !< Longitude (degree E)
  real, intent(in) :: uvel !< Zonal velocity (m/s)
  real, intent(in) :: vvel !< Meridional velocity (m/s)
  real, intent(out) :: xdot !< x-component of velocity in tangent plane (m/s)
  real, intent(out) :: ydot !< y-component of velocity in tangent plane (m/s)
  ! Local variables
  real :: clon,slon

  clon=cos(lon*pi_180)
  slon=sin(lon*pi_180)
  xdot=-slon*uvel-clon*vvel
  ydot=clon*uvel-slon*vvel

end subroutine rotvec_to_tang

! ####################################################################################
!> Returns metric converting distance in meters to grid distance
subroutine  convert_from_meters_to_grid(lat_ref,grid_is_latlon ,dlon_dx,dlat_dy)
  ! Arguments
  real, intent(in) :: lat_ref !< Latitude at which to make metric conversion (degree N)
  logical, intent(in) :: grid_is_latlon !< True if grid model grid is in lat-lon coordinates
  real, intent(out) :: dlon_dx !< Metric dlon/dx
  real, intent(out) :: dlat_dy !< Metric dlat/dy

  if (grid_is_latlon) then
    dlon_dx=(180./pi)/(Rearth*cos((lat_ref)*(pi/180.)))
    dlat_dy=(180./pi)/Rearth
  else
    dlon_dx=1.
    dlat_dy=1.
  endif

end subroutine convert_from_meters_to_grid

! ###################################################################################
!> Returns metric converting grid distances to meters
subroutine convert_from_grid_to_meters(lat_ref, grid_is_latlon, dx_dlon, dy_dlat)
  ! Arguments
  real, intent(in) :: lat_ref !< Latitude at which to make metric conversion (degree N)
  logical, intent(in) :: grid_is_latlon !< True if grid model grid is in lat-lon coordinates
  real, intent(out) :: dx_dlon !< Metric dx/dlon
  real, intent(out) :: dy_dlat !< Metric dy/dlat

  if (grid_is_latlon) then
    dx_dlon=(pi/180.)*Rearth*cos((lat_ref)*(pi/180.))
    dy_dlat=(pi/180.)*Rearth
  else
    dx_dlon=1.
    dy_dlat=1.
  endif

end subroutine convert_from_grid_to_meters

! ###################################################################################
!> Moves part's cell indexes,(i,j), checking for collisional with coasts
subroutine adjust_index_and_ground(grd, lon, lat, uvel, vvel, i, j, xi, yj, bounced, error, id)
  ! Arguments
  type(particles_gridded), pointer :: grd !< Container for gridded fields
  real, intent(inout) :: lon !< Longitude (degree E)
  real, intent(inout) :: lat !< Latitude (degree N)
  real, intent(inout) :: uvel !< Zonal velocity (m/s)
  real, intent(inout) :: vvel !< Meridional velocity (m/s)
  real, intent(inout) :: xi !< Non-dimension x-position within cell
  real, intent(inout) :: yj !< Non-dimension y-position within cell
  integer, intent(inout) :: i !< i-index of cell
  integer, intent(inout) :: j !< j-index of cell
  logical, intent(out) :: bounced !< True if part collided with coast
  logical, intent(out) :: error !< True if adjustments could not be made consistently
  integer(kind=8), intent(in) :: id !< part identifier
  ! Local variables
  logical lret, lpos
  real, parameter :: posn_eps=0.05
  integer :: icount, i0, j0, inm, jnm
  real :: xi0, yj0, lon0, lat0
  integer :: stderrunit
  logical :: point_in_cell_using_xi_yj

  ! Get the stderr unit number
  stderrunit = stderr()

  bounced=.false.
  error=.false.
  lon0=lon; lat0=lat ! original position
  i0=i; j0=j ! original i,j
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
!  print *, 'Alon:', lon, lat, i, j, xi, yj, lret
  xi0=xi; yj0=yj ! original xi,yj


  !Removing this while debuggin
  if (debug) then
    !Sanity check lret, xi and yj
    lret=is_point_in_cell(grd, lon, lat, i, j)
    point_in_cell_using_xi_yj=is_point_within_xi_yj_bounds(xi,yj)
    if (.not. point_in_cell_using_xi_yj) then

      if (lret) then
        write(stderrunit,*) 'diamonds, adjust: WARNING!!! lret=T but |xi,yj|>1',mpp_pe()
        write(stderrunit,*) 'diamonds, adjust: xi=',xi,' lon=',lon
        write(stderrunit,*) 'diamonds, adjust: x3 x2=',grd%lon(i-1,j),grd%lon(i,j)
        write(stderrunit,*) 'diamonds, adjust: x0 x1=',grd%lon(i-1,j-1),grd%lon(i,j-1)
        write(stderrunit,*) 'diamonds, adjust: yi=',yj,' lat=',lat
        write(stderrunit,*) 'diamonds, adjust: y3 y2=',grd%lat(i-1,j),grd%lat(i,j)
        write(stderrunit,*) 'diamonds, adjust: y0 y1=',grd%lat(i-1,j-1),grd%lat(i,j-1)
        lret=is_point_in_cell(grd, lon, lat, i, j,explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj,explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        call error_mesg('adjust index, ','particle is_point_in_cell=True but xi, yi are out of cell',FATAL)
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
        lret=is_point_in_cell(grd, lon, lat, i, j,  explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn is_point_in_cell=',lret
        lret=pos_within_cell(grd, lon, lat, i, j, xi, yj, explain=.true.)
        write(stderrunit,*) 'diamonds, adjust: fn pos_within_cell=',lret
        write(0,*) 'This should never happen!'
        call error_mesg('adjust index, ','particle is_point_in_cell=False but xi, yi are out of cell',FATAL)
        error=.true.; return
      endif
    endif
    lret=pos_within_cell(grd, lon, lat, i, j, xi, yj)
  endif ! debug

  if (lret) return ! part was already in cell

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
!    elseif (xi.ge.1.) then   !Alon: maybe it should be .ge.
      if (inm<grd%ied) then
        inm=inm+1
      endif
    endif
    if (yj.lt.0.) then
      if (jnm>grd%jsd) then
        jnm=jnm-1
      endif
    elseif (yj.gt.1.) then
!    elseif (yj.ge.1.) then   !Alon:maybe it should be .ge.
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
    elseif (xi.ge.1.) then    !Alon!!!!
!    elseif (xi.gt.1.) then
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
    elseif (yj.ge.1.) then     !Alon.
!    elseif (yj.gt.1.) then
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
      if (xi>=1.) xi=1.-posn_eps   !Alon.
!      if (xi>1.) xi=1.-posn_eps   !
      if (xi<0.) xi=posn_eps
      if (yj>=1.) yj=1.-posn_eps  !Alon.
!      if (yj>1.) yj=1.-posn_eps
      if (yj<0.) yj=posn_eps
      lon=bilin(grd, grd%lon, i, j, xi, yj)
      lat=bilin(grd, grd%lat, i, j, xi, yj)
    endif
    if (debug) then
      if (grd%msk(i,j)==0.) stop 'diamonds, adjust: part is in land! This should not happen...'
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
    if (debug) then
      write(stderrunit,*) 'diamonds, adjust: lon0, lat0=',lon0,lat0
      write(stderrunit,*) 'diamonds, adjust: xi0, yj0=',xi0,yj0
      write(stderrunit,*) 'diamonds, adjust: i0,j0=',i0,j0
      write(stderrunit,*) 'diamonds, adjust: lon, lat=',lon,lat
      write(stderrunit,*) 'diamonds, adjust: xi,yj=',xi,yj
      write(stderrunit,*) 'diamonds, adjust: i,j=',i,j
      write(stderrunit,*) 'diamonds, adjust: inm,jnm=',inm,jnm
      write(stderrunit,*) 'diamonds, adjust: icount=',icount
      lret=pos_within_cell(grd, lon, lat, i, j, xi, yj,explain=.true.)
      write(stderrunit,*) 'diamonds, adjust: lret=',lret
    endif

    if (abs(i-i0)+abs(j-j0)==0) then
      if (use_roundoff_fix) then
        ! This is a special case due to round off where is_point_in_cell()
        ! returns false but xi and yj are between 0 and 1.
        ! It occurs very rarely but often enough to have brought down
        ! ESM2G four times since the spin-up began. (as of 8/10/2010)
        ! This temporary fix arbitrarily moves the part toward the
        ! center of the current cell.
        xi=(xi-0.5)*(1.-posn_eps)+0.5
        yj=(yj-0.5)*(1.-posn_eps)+0.5
      endif
      call error_mesg('diamonds, adjust', 'part did not move or bounce during iterations AND was not in cell. Adjusting!', WARNING)
      write(stderrunit,*) 'diamonds, adjust: The adjusting particle is: ', id,  mpp_pe()
      write(stderrunit,*) 'diamonds, adjust: The adjusting lon,lat,u,v: ', lon, lat, uvel, vvel
      write(stderrunit,*) 'diamonds, adjust: The adjusting xi,ji: ', xi, yj
      lret=pos_within_cell(grd, lon, lat, inm, jnm, xi, yj,explain=.true.)
    else
      call error_mesg('diamonds, adjust', 'part iterated many times without bouncing!', WARNING)
    endif
  endif
!  if (xi>1.) xi=1.-posn_eps    !Alon
  if (xi>=1.) xi=1.-posn_eps
  if (xi<0.) xi=posn_eps
  if (yj>1.) yj=1.-posn_eps
!  if (yj>1.) yj=1.-posn_eps
  if (yj<=0.) yj=posn_eps        !Alon
  lon=bilin(grd, grd%lon, i, j, xi, yj)
  lat=bilin(grd, grd%lat, i, j, xi, yj)
  lret=pos_within_cell(grd, lon, lat, i, j, xi, yj) ! Update xi and yj

  if (.not. lret) then
    write(0,*) 'i0, j0,=', i0,j0
    write(0,*) 'xi0, yj0,=', xi0,yj0
    write(0,*) 'grd%msk(i0, j0)=', grd%msk(i0,j0)
    write(0,*) 'lon0, lat0,=', lon0,lat0
    write(0,*) 'i,j,lon, lat,grd%msk(i,j)=', i,j,lon,lat,grd%msk(i,j)
    write(stderrunit,*) 'diamonds, adjust: Should not get here! part is not in cell after adjustment', id, mpp_pe()
    if (debug) error=.true.
  endif
 end subroutine adjust_index_and_ground

! ###################################################################################

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

  !call write_trajectory(parts%trajectories)

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

!> Calculates parts%grd%mass_on_ocean
subroutine calculate_mass_on_ocean(parts, with_diagnostics)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  logical, intent(in) :: with_diagnostics
  ! Local variables
  type(particle), pointer :: part
  type(particles_gridded), pointer :: grd
  integer :: grdj, grdi
  integer :: j, i

  ! For convenience
  grd=>parts%grd

  !Initialize fields
  grd%mass_on_ocean(:,:,:)=0.
  grd%area_on_ocean(:,:,:)=0.
  grd%Uvel_on_ocean(:,:,:)=0.
  grd%Vvel_on_ocean(:,:,:)=0.

  do grdj = grd%jsc-1,grd%jec+1 ; do grdi = grd%isc-1,grd%iec+1
    part=>parts%list(grdi,grdj)%first
    do while(associated(part))
      i=part%ine  ;     j=part%jne
      if (grd%area(i,j) > 0.) then

        !Increasing Mass on ocean
        if ((parts%add_weight_to_ocean .and. .not. parts%time_average_weight) .or.(parts%find_melt_using_spread_mass)) then
          call spread_mass_across_ocean_cells(parts, part, part%ine, part%jne, part%xi, part%yj, part%mass,part%mass_of_bits, part%mass_scaling, &
                part%length*part%width, part%thickness)
        endif

        !Calculated some particle diagnositcs
        if (with_diagnostics) call calculate_sum_over_parts_diagnositcs(parts,grd,part,i,j)

      endif
      part=>part%next
    enddo
  enddo ;enddo

end subroutine calculate_mass_on_ocean

!> Projects additional diagnostics of parts on to the grid
subroutine calculate_sum_over_parts_diagnositcs(parts, grd, part, i, j)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  type(particles_gridded), pointer :: grd !< Container for gridded fields
  type(particle), pointer :: part !< An particle
  integer, intent(in) :: i !< i-index of cell containing part
  integer, intent(in) :: j !< j-index of cell containing part
  ! Local variables
  real ::  Abits, Lbits, Mbits

  !Virtual area diagnostic
  if (grd%id_virtual_area>0) then
    if (parts%party_bit_erosion_fraction>0.) then
      Lbits=min(part%length,part%width,part%thickness,40.) ! assume party bits are smallest dimension or 40 meters
      Abits=(part%mass_of_bits/parts%rho_parts)/Lbits ! Effective bottom area (assuming T=Lbits)
    else
      Abits=0.0
    endif
    grd%virtual_area(i,j)=grd%virtual_area(i,j)+(part%width*part%length+Abits)*part%mass_scaling ! m^2
  endif

  !Mass diagnostic (also used in u_particle, v_particle
  if ((grd%id_mass>0 ) .or. ((grd%id_u_particle>0) .or. (grd%id_v_particle>0)))   &
       & grd%mass(i,j)=grd%mass(i,j)+part%mass/grd%area(i,j)*part%mass_scaling ! kg/m2

  !Finding the average particle velocity in a grid cell (mass weighted)
  if (grd%id_u_particle>0) &
  grd%u_particle(i,j)=grd%u_particle(i,j)+((part%mass/grd%area(i,j)*part%mass_scaling)*part%uvel) ! kg/m2
  if (grd%id_v_particle>0) &
  grd%v_particle(i,j)=grd%v_particle(i,j)+((part%mass/grd%area(i,j)*part%mass_scaling)*part%vvel) ! kg/m2

  !Mass of party bits
  if (grd%id_party_mass>0 .or. parts%add_weight_to_ocean)&
       & grd%party_mass(i,j)=grd%party_mass(i,j)+part%mass_of_bits/grd%area(i,j)*part%mass_scaling ! kg/m2
end subroutine calculate_sum_over_parts_diagnositcs

!> Spread mass of a part around cells centered on i,j
subroutine spread_mass_across_ocean_cells(parts, part, i, j, x, y, Mpart, Mbits, scaling, Area, Tn)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  type(particle), pointer :: part !< part whose mass is being considered
  integer, intent(in) :: i !< i-index of cell contained center of part
  integer, intent(in) :: j !< j-index of cell contained center of part
  real, intent(in) :: x !< Longitude of part (degree E)
  real, intent(in) :: y !< Latitude of part (degree N)
  real, intent(in) :: Mpart !< Mass of part (kg)
  real, intent(in) :: Mbits !< Mass of party bits (kg)
  real, intent(in) :: scaling !< Multiplier to scale mass (nondim)
  real, intent(in) :: Area !< Area of part (m2)
  real, intent(in) :: Tn !< Thickness of part (m)
  ! Local variables
  type(particles_gridded), pointer :: grd
  real :: xL, xC, xR, yD, yC, yU, Mass, L
  real :: yDxL, yDxC, yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR
  real :: S, H, origin_x, origin_y, x0, y0
  real :: Area_Q1,Area_Q2 , Area_Q3,Area_Q4, Area_hex
  real :: fraction_used !fraction of particle mass included (part of the mass near the boundary is discarded sometimes)
  real :: I_fraction_used !Inverse of fraction used
  real :: tol
  real :: Dn, Hocean
  real, parameter :: rho_seawater=1035.
  integer :: stderrunit
  logical :: debug
  real :: orientation, Mass_part

  ! Get the stderr unit number
  stderrunit = stderr()

  tol=1.e-10
  grd=>parts%grd
  Mass_part=Mpart

  ! Trimming particles to account for grounded fraction.
  if (parts%grounding_fraction>0.) then
    Hocean=parts%grounding_fraction*(grd%ocean_depth(i,j)+grd%ssh(i,j))
    Dn=(parts%rho_parts/rho_seawater)*Tn ! re-calculate draught (keel depth)
    if (Dn>Hocean) Mass_part=Mpart*min(1.,Hocean/Dn)
  endif

  Mass=(Mass_part+Mbits)*scaling
  ! This line attempts to "clip" the weight felt by the ocean. The concept of
  ! clipping is non-physical and this step should be replaced by grounding.
  if (grd%clipping_depth>0.) Mass=min(Mass,grd%clipping_depth*grd%area(i,j)*rho_seawater)

  !Initialize weights for each cell
  yDxL=0.  ; yDxC=0. ; yDxR=0. ; yCxL=0. ; yCxR=0.
  yUxL=0.  ; yUxC=0. ; yUxR=0. ; yCxC=1.

  if (.not. parts%hexagonal_particles) then ! Treat particles as rectangles of size L: (this is the default)

    ! L is the non dimensional length of the particle [ L=(Area of part/ Area of grid cell)^0.5 ] or something like that.
    if (grd%area(i,j)>0) then
      L=min( sqrt(Area / grd%area(i,j)),1.0)
    else
      L=1.
    endif

    if (parts%use_old_spreading) then
      ! Old version before particles were given size L
      xL=min(0.5, max(0., 0.5-x))
      xR=min(0.5, max(0., x-0.5))
      xC=max(0., 1.-(xL+xR))
      yD=min(0.5, max(0., 0.5-y))
      yU=min(0.5, max(0., y-0.5))
      yC=max(0., 1.-(yD+yU))
    else
      xL=min(0.5, max(0., 0.5-(x/L)))
      xR=min(0.5, max(0., (x/L)+(0.5-(1/L) )))
      xC=max(0., 1.-(xL+xR))
      yD=min(0.5, max(0., 0.5-(y/L)))
      yU=min(0.5, max(0., (y/L)+(0.5-(1/L) )))
      yC=max(0., 1.-(yD+yU))
    endif

    yDxL=yD*xL*grd%msk(i-1,j-1)
    yDxC=yD*xC*grd%msk(i  ,j-1)
    yDxR=yD*xR*grd%msk(i+1,j-1)
    yCxL=yC*xL*grd%msk(i-1,j  )
    yCxR=yC*xR*grd%msk(i+1,j  )
    yUxL=yU*xL*grd%msk(i-1,j+1)
    yUxC=yU*xC*grd%msk(i  ,j+1)
    yUxR=yU*xR*grd%msk(i+1,j+1)
    yCxC=1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )

    fraction_used=1. ! rectangular parts do share mass with boundaries (all mass is included in cells)

  else ! Spread mass as if elements area hexagonal

    orientation=parts%initial_orientation
    !if ((parts%particle_bonds_on) .and. (parts%rotate_particles_for_mass_spreading)) call find_orientation_using_particle_bonds(grd,part,orientation)

    if (grd%area(i,j)>0) then
      H=min(( (sqrt(Area/(2.*sqrt(3.))) / sqrt(grd%area(i,j)))),1.) ! Non-dimensionalize element length by grid area. (This gives the non-dim Apothem of the hexagon)
    else
      H=(sqrt(3.)/2)*(0.49) ! Largest allowable H, since this makes S=0.49, and S has to be less than 0.5 (Not sure what the implications of this are)
    endif
    S=(2/sqrt(3.))*H !Side of the hexagon

    if (S>0.5) then
      ! The width of an particle should not be greater than half the grid cell, or else it can spread over 3 cells  (i.e. S must be less than 0.5 non-dimensionally)
      !print 'Elements must be smaller than a whole grid cell', 'i.e.: S= ' , S , '>=0.5'
      call error_mesg('diamonds, hexagonal spreading', 'Diameter of the particle is larger than a grid cell. Use smaller particles', WARNING)
    endif

    !Subtracting the position of the nearest corner from x,y  (The mass will then be spread over the 4 cells connected to that corner)
    origin_x=1. ; origin_y=1.
    if (x<0.5) origin_x=0.
    if (y<0.5) origin_y=0.

    !Position of the hexagon center, relative to origin at the nearest vertex
    x0=(x-origin_x)
    y0=(y-origin_y)

    call Hexagon_into_quadrants_using_triangles(x0,y0,H,orientation,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)

    if (min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4)) <-tol) then
      call error_mesg('diamonds, hexagonal spreading', 'Intersection with hexagons should not be negative!!!', WARNING)
      write(stderrunit,*) 'diamonds, yU,yC,yD', Area_Q1, Area_Q2, Area_Q3, Area_Q4
    endif

    Area_Q1=Area_Q1/Area_hex
    Area_Q2=Area_Q2/Area_hex
    Area_Q3=Area_Q3/Area_hex
    Area_Q4=Area_Q4/Area_hex

    !Now, you decide which quadrant belongs to which mass on ocean cell.
    if ((x.ge. 0.5) .and. (y.ge. 0.5)) then !Top right vertex
      yUxR=Area_Q1
      yUxC=Area_Q2
      yCxC=Area_Q3
      yCxR=Area_Q4
    elseif ((x .lt. 0.5) .and. (y.ge. 0.5)) then  !Top left vertex
      yUxC=Area_Q1
      yUxL=Area_Q2
      yCxL=Area_Q3
      yCxC=Area_Q4
    elseif ((x.lt.0.5) .and. (y.lt. 0.5)) then !Bottom left vertex
      yCxC=Area_Q1
      yCxL=Area_Q2
      yDxL=Area_Q3
      yDxC=Area_Q4
    elseif ((x.ge.0.5) .and. (y.lt. 0.5)) then!Bottom right vertex
      yCxR=Area_Q1
      yCxC=Area_Q2
      yDxC=Area_Q3
      yDxR=Area_Q4
    endif

    !Temporary for debugging reasons.
    if (mpp_pe()==mpp_root_pe()) then
      !write(stderrunit,*) 'diamonds, You are in the hexagonal domain now!!!'
    endif

    !Double check that all the mass is being used.
    if ((abs(yCxC-(1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )))>tol) .and. (mpp_pe().eq. mpp_root_pe())) then
      !call error_mesg('diamonds, hexagonal spreading', 'All the mass is not being used!!!', WARNING)
      write(stderrunit,*) 'diamonds, hexagonal, H,x0,y0', H, x0 , y0
      write(stderrunit,*) 'diamonds, hexagonal, Areas',(Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,  Area_Q4
      debug=.True.
      !call Hexagon_into_quadrants_using_triangles(x0,y0,H,orientation,Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4, debug)
      call error_mesg('diamonds, hexagonal spreading', 'All the mass is not being used!!!', FATAL)
    endif

    !Scale each cell by (1/fraction_used) in order to redisribute ice mass which landed up on the land, back into the ocean
    !Note that for the square elements, the mass has already been reassigned, so fraction_used shoule be equal to 1 aready
    fraction_used= ((yDxL*grd%msk(i-1,j-1)) + (yDxC*grd%msk(i  ,j-1))  +(yDxR*grd%msk(i+1,j-1)) +(yCxL*grd%msk(i-1,j  )) +  (yCxR*grd%msk(i+1,j  ))&
                   +(yUxL*grd%msk(i-1,j+1)) +(yUxC*grd%msk(i  ,j+1))   +(yUxR*grd%msk(i+1,j+1)) + (yCxC**grd%msk(i,j)))
    if  (part%static_part .eq. 1)  fraction_used=1.  !Static particles do not share their mass with the boundary
                                                ! (this allows us to easily  initialize hexagonal particles in regular arrangements against boundaries)
  endif
  I_fraction_used=1./fraction_used !Invert this so that the arithmatec reprocudes

  !Spreading the particle mass onto the ocean
  call spread_variable_across_cells(grd, grd%mass_on_ocean, Mass, i ,j, &
             yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, I_fraction_used)
  !Spreading the particle area onto the ocean
  call spread_variable_across_cells(grd, grd%area_on_ocean, Area*scaling , i ,j, &
             yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR,I_fraction_used)
  !Spreading the particle x momentum onto the ocean
  call spread_variable_across_cells(grd,grd%Uvel_on_ocean, part%uvel*Area*scaling , i ,j, &
             yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, I_fraction_used)
  !Spreading the particle y momentum onto the ocean
  call spread_variable_across_cells(grd,grd%Vvel_on_ocean, part%vvel*Area*scaling , i ,j, &
             yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, I_fraction_used)

end subroutine spread_mass_across_ocean_cells

!> Distribute a quantity among nine cells on a grid centered at cell i,j
subroutine spread_variable_across_cells(grd, variable_on_ocean, Var, i, j, &
           yDxL, yDxC,yDxR, yCxL, yCxC, yCxR, yUxL, yUxC, yUxR, I_fraction_used)
  ! Arguments
  type(particles_gridded), pointer, intent(in) :: grd !< Container for gridded fields
  real, dimension(grd%isd:grd%ied, grd%jsd:grd%jed, 9), intent(inout) :: variable_on_ocean !< Gridded field to augment
  real, intent(in) :: Var !< Variable to be spread accross cell
  real, intent(in) :: yDxL !< Weight for the cell at i-1,j-1
  real, intent(in) :: yDxC !< Weight for the cell at i-1,j
  real, intent(in) :: yDxR !< Weight for the cell at i-1,j+1
  real, intent(in) :: yCxL !< Weight for the cell at i,j-1
  real, intent(in) :: yCxC !< Weight for the cell at i,j
  real, intent(in) :: yCxR !< Weight for the cell at i,j-1
  real, intent(in) :: yUxL !< Weight for the cell at i+1,j-1
  real, intent(in) :: yUxC !< Weight for the cell at i+1,j
  real, intent(in) :: yUxR !< Weight for the cell at i+1,j+1
  real, intent(in) :: I_fraction_used !< Amount of particle used (inverse)
  integer, intent(in) :: i !< i-index of cell containing center of part
  integer, intent(in) :: j !< j-index of cell containing center of part

  !Spreading the particle mass onto the ocean
  variable_on_ocean(i,j,1)=variable_on_ocean(i,j,1)+(yDxL*Var*I_fraction_used)
  variable_on_ocean(i,j,2)=variable_on_ocean(i,j,2)+(yDxC*Var*I_fraction_used)
  variable_on_ocean(i,j,3)=variable_on_ocean(i,j,3)+(yDxR*Var*I_fraction_used)
  variable_on_ocean(i,j,4)=variable_on_ocean(i,j,4)+(yCxL*Var*I_fraction_used)
  variable_on_ocean(i,j,5)=variable_on_ocean(i,j,5)+(yCxC*Var*I_fraction_used)
  variable_on_ocean(i,j,6)=variable_on_ocean(i,j,6)+(yCxR*Var*I_fraction_used)
  variable_on_ocean(i,j,7)=variable_on_ocean(i,j,7)+(yUxL*Var*I_fraction_used)
  variable_on_ocean(i,j,8)=variable_on_ocean(i,j,8)+(yUxC*Var*I_fraction_used)
  variable_on_ocean(i,j,9)=variable_on_ocean(i,j,9)+(yUxR*Var*I_fraction_used)

end subroutine spread_variable_across_cells


!> Returns true of a point is in or on the rectangle with opposite corners A and B
logical function point_in_interval(Ax, Ay, Bx, By, px, py)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: px !< x-position of point
  real, intent(in) :: py !< y-position of point
  point_in_interval=.False.
  if ((px <= max(Ax,Bx)) .and. (px >= min(Ax,Bx))) then
    if ((py <= max(Ay,By)) .and. (py >= min(Ay,By))) then
      point_in_interval=.True.
    endif
  endif
end function point_in_interval

!> Returns true if point q is on a line through points A and B
logical function point_is_on_the_line(Ax, Ay, Bx, By, qx, qy)
  ! Arguments
  real, intent(in) :: Ax !< x-position of point A
  real, intent(in) :: Ay !< y-position of point A
  real, intent(in) :: Bx !< x-position of point B
  real, intent(in) :: By !< y-position of point B
  real, intent(in) :: qx !< x-position of point q
  real, intent(in) :: qy !< y-position of point q
  ! Local variables
  real :: tol, dxc,dyc,dxl,dyl,cross
  !tol=1.e-12
  tol=0.0
  dxc = qx - Ax
  dyc = qy - Ay
  dxl = Bx - Ax
  dyl = By - Ay
  cross = dxc * dyl - dyc * dxl
  if (abs(cross)<=tol) then
    point_is_on_the_line=.True.
  else
   point_is_on_the_line=.False.
  endif
end function point_is_on_the_line

!> Returns True if a point q is inside a triangle ABC
!!
!! This function decides whether a point (qx,qy) is inside the triangle ABC.
!! There is also the option to include the boundary of the triangle.
logical function point_in_triangle(Ax, Ay, Bx, By, Cx, Cy, qx, qy)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  real, intent(in) :: qx !< x-position of point q
  real, intent(in) :: qy !< y-position of point q
  ! Local variables
  real :: l0,l1,l2,p0,p1,p2
  real :: v0x,v1x,v2x,v0y,v1y,v2y,dot00,dot01,dot02,dot11,dot12

  point_in_triangle = .False.
  if ((Ax==qx .and. Ay==qy) .or. (Bx==qx .and. By==qy) .or. (Cx==qx .and. Cy==qy)) then ! Exclude the pathelogical case
      point_in_triangle = .False.
  else
    if (((point_is_on_the_line(Ax,Ay,Bx,By,qx,qy) .or. (point_is_on_the_line(Ax,Ay,Cx,Cy,qx,qy))) .or. (point_is_on_the_line(Bx,By,Cx,Cy,qx,qy)))) then
      point_in_triangle = .False.
    else
      ! Compute point in triangle using Barycentric coordinates (the same as sum_sign_dot_prod routines)
      l0=(qx-Ax)*(By-Ay)-(qy-Ay)*(Bx-Ax)
      l1=(qx-Bx)*(Cy-By)-(qy-By)*(Cx-Bx)
      l2=(qx-Cx)*(Ay-Cy)-(qy-Cy)*(Ax-Cx)

      p0=sign(1., l0); if (l0==0.)  p0=0.
      p1=sign(1., l1); if (l1==0.)  p1=0.
      p2=sign(1., l2); if (l2==0.)  p2=0.

      if ( (abs(p0)+abs(p2))+(abs(p1)) == abs((p0+p2)+(p1)) )  point_in_triangle = .True.
    endif
  endif
end function point_in_triangle

!> Returns area of a triangle
real function Area_of_triangle(Ax, Ay, Bx, By, Cx, Cy)
  ! Arguments
  real, intent(in) :: Ax !< x-position of corner A
  real, intent(in) :: Ay !< y-position of corner A
  real, intent(in) :: Bx !< x-position of corner B
  real, intent(in) :: By !< y-position of corner B
  real, intent(in) :: Cx !< x-position of corner C
  real, intent(in) :: Cy !< y-position of corner C
  Area_of_triangle    =   abs(    0.5*((Ax*(By-Cy))+(Bx*(Cy-Ay))+(Cx*(Ay-By))) )
end function Area_of_triangle

!> Sums up a part property into a gridded field
subroutine sum_up_spread_fields(parts, field, field_name)
  ! Arguments
  type(particles), pointer :: parts !< Container for all types and memory
  real, dimension(parts%grd%isc:parts%grd%iec,parts%grd%jsc:parts%grd%jec), intent(out) :: field !< Gridded field
  character(len=4), intent(in) :: field_name !< Name of field to grid
  ! Local variables
  integer :: i, j
  type(particles_gridded), pointer :: grd
  real :: dmda
  logical :: lerr
  real, dimension(parts%grd%isd:parts%grd%ied, parts%grd%jsd:parts%grd%jed,9) :: var_on_ocean   !Variable being spread onto the ocean  (mass, area, Uvel, Vvel)
  integer :: stderrunit

  ! Get the stderr unit number
  stderrunit = stderr()
  ! For convenience
  grd=>parts%grd

  field(:,:)=0.

  !Deciding which varibale to spread across cells across grid cells
  if (field_name=='mass') var_on_ocean(:,:,:)=grd%mass_on_ocean(:,:,:)
  if (field_name=='area') var_on_ocean(:,:,:)=grd%area_on_ocean(:,:,:)
  if (field_name=='Uvel') var_on_ocean(:,:,:)=grd%Uvel_on_ocean(:,:,:)
  if (field_name=='Vvel') var_on_ocean(:,:,:)=grd%Vvel_on_ocean(:,:,:)

  !This line has been removed, for that routine can be used for other fields
  !if (.not. parts%add_weight_to_ocean) return

  !Update the halos of the var_on_ocean
  call mpp_update_domains(var_on_ocean, grd%domain)

  !Rotatine when old_bug_rotated_weights is on  - we should remove this.
  if (.not. old_bug_rotated_weights) then
    do j=grd%jsd, grd%jed; do i=grd%isd, grd%ied
      if (grd%parity_x(i,j)<0.) then
        ! This block assumes both parity_x and parity_y are negative
        ! (i.e. a 180 degree rotation). In general, we should handle
        ! +/- 90 degree rotations as well but in CM2*-class models
        ! this is not necessary. -aja
        dmda=var_on_ocean(i,j,9); var_on_ocean(i,j,9)=var_on_ocean(i,j,1); var_on_ocean(i,j,1)=dmda
        dmda=var_on_ocean(i,j,8); var_on_ocean(i,j,8)=var_on_ocean(i,j,2); var_on_ocean(i,j,2)=dmda
        dmda=var_on_ocean(i,j,7); var_on_ocean(i,j,7)=var_on_ocean(i,j,3); var_on_ocean(i,j,3)=dmda
        dmda=var_on_ocean(i,j,6); var_on_ocean(i,j,6)=var_on_ocean(i,j,4); var_on_ocean(i,j,4)=dmda
      endif
    enddo; enddo
  endif

  !Here we add the contribution of the 9 cells. This is the heart of the routine.
  do j=grd%jsc, grd%jec; do i=grd%isc, grd%iec
    dmda=var_on_ocean(i,j,5) &
         + ( ( (var_on_ocean(i-1,j-1,9)+var_on_ocean(i+1,j+1,1))   &
         +     (var_on_ocean(i+1,j-1,7)+var_on_ocean(i-1,j+1,3)) ) &
         +   ( (var_on_ocean(i-1,j  ,6)+var_on_ocean(i+1,j  ,4))   &
         +     (var_on_ocean(i  ,j-1,8)+var_on_ocean(i  ,j+1,2)) ) )
    if (grd%area(i,j)>0) dmda=dmda/grd%area(i,j)*grd%msk(i,j)

    !Make sure that area <=1.0
    if (field_name=='area') dmda=min(dmda,1.0)

    field(i,j)=dmda
  enddo; enddo

  if (debug) then
    grd%tmp(:,:)=0.; grd%tmp(grd%isc:grd%iec,grd%jsc:grd%jec)=field
    if (field_name=='mass') then
      call grd_chksum3(grd, grd%mass_on_ocean, 'mass parts (incr)')
      call grd_chksum2(grd, grd%tmp, 'mass out (incr)')
    elseif (field_name=='area') then
      call grd_chksum3(grd, grd%area_on_ocean, 'area parts (incr)')
      call grd_chksum2(grd, grd%tmp, 'area out (incr)')
    endif
 endif
end subroutine sum_up_spread_fields

!> Approximately convert a wind-stress into a velocity difference
subroutine invert_tau_for_du(u, v)
  ! Arguments
  real, dimension(:,:), intent(inout) :: u !< On entry, zonal wind stress (Pa). On exit, zonal velocity difference (m/s).
  real, dimension(:,:), intent(inout) :: v !< On entry, meridional wind stress (Pa). On exit, meridional velocity difference (m/s).
  ! Local variables
  integer :: i, j
  real :: cd, cddvmod, tau2

  cd=0.0015

  do j=lbound(u,2), ubound(u,2)
    do i=lbound(u,1), ubound(u,1)
      tau2=u(i,j)*u(i,j)+v(i,j)*v(i,j)
      cddvmod=sqrt(cd*sqrt(tau2))
      if (cddvmod.ne.0.) then
        u(i,j)=u(i,j)/cddvmod
        v(i,j)=v(i,j)/cddvmod
      else
        u(i,j)=0.
        v(i,j)=0.
      endif
    enddo
  enddo

end subroutine invert_tau_for_du

!> Prints consistency summary of start and end states
subroutine report_consistant(budgetstr, budgetunits, startstr, startval, endstr, endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: budgetunits !< Units of budgeted quantity
  character*(*), intent(in) :: startstr !< Start label
  real, intent(in) :: startval !< Start value for budget
  character*(*), intent(in) :: endstr !< End label
  real, intent(in) :: endval !< End value for budget
  ! Local variables
  write(*,200) budgetstr//' check:', &
                      startstr,startval,budgetunits, &
                      endstr,endval,budgetunits, &
                      'error',(endval-startval)/((endval+startval)+1e-30),'nd'
  200 format("diamonds: ",a19,10(a18,"=",es14.7,x,a2,:,","))
end subroutine report_consistant

!> Prints summary of start and end states
subroutine report_state(budgetstr, budgetunits, startstr, startval, endstr, endval, delstr, nparts)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: budgetunits !< Units of budgeted quantity
  character*(*), intent(in) :: startstr !< Start label
  real, intent(in) :: startval !< Start value for budget
  character*(*), intent(in) :: endstr !< End label
  real, intent(in) :: endval !< End value for budget
  character*(*), intent(in) :: delstr !< Delta label
  integer, intent(in), optional :: nparts !< Number of parts
  ! Local variables
  if (present(nparts)) then
    write(*,100) budgetstr//' state:', &
                        startstr//' start',startval,budgetunits, &
                        endstr//' end',endval,budgetunits, &
                        'Delta '//delstr,endval-startval,budgetunits, &
                        '# of parts',nparts
  else
    write(*,100) budgetstr//' state:', &
                        startstr//' start',startval,budgetunits, &
                        endstr//' end',endval,budgetunits, &
                        delstr//'Delta',endval-startval,budgetunits
  endif
  100 format("diamonds: ",a19,3(a18,"=",es14.7,x,a2,:,","),a12,i8)
end subroutine report_state


!> Prints a budget
subroutine report_budget(budgetstr, budgetunits, instr, inval, outstr, outval, delstr, startval, endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: budgetunits !< Units of budgeted quantity
  character*(*), intent(in) :: instr !< Incoming label
  real, intent(in) :: inval !< Incoming value
  character*(*), intent(in) :: outstr !< Outgoing label
  real, intent(in) :: outval !< Outgoing value
  character*(*), intent(in) :: delstr !< Delta label
  real, intent(in) :: startval !< Start value for budget
  real, intent(in) :: endval !< End value for budget
  ! Local variables
  write(*,200) budgetstr//' budget:', &
                      instr//' in',inval,budgetunits, &
                      outstr//' out',outval,budgetunits, &
                      'Delta '//delstr,inval-outval,budgetunits, &
                      'error',((endval-startval)-(inval-outval))/max(1.e-30,max(abs(endval-startval),abs(inval-outval))),'nd'
  200 format("diamonds: ",a19,3(a18,"=",es14.7,x,a2,:,","),a8,"=",es10.3,x,a2)
end subroutine report_budget

!> Prints summary of start and end states
subroutine report_istate(budgetstr, startstr, startval, endstr, endval, delstr)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: startstr !< Start label
  integer, intent(in) :: startval !< Start value for budget
  character*(*), intent(in) :: endstr !< End label
  integer, intent(in) :: endval !< End value for budget
  character*(*), intent(in) :: delstr !< Delta label
  ! Local variables
  write(*,100) budgetstr//' state:', &
                        startstr//' start',startval, &
                        endstr//' end',endval, &
                        delstr//'Delta',endval-startval
  100 format("diamonds: ",a19,3(a18,"=",i14,x,:,","))
end subroutine report_istate

!> Prints a budget
subroutine report_ibudget(budgetstr,instr,inval,outstr,outval,delstr,startval,endval)
  ! Arguments
  character*(*), intent(in) :: budgetstr !< Budget title
  character*(*), intent(in) :: instr !< Incoming label
  integer, intent(in) :: inval !< Incoming value
  character*(*), intent(in) :: outstr !< Outgoing label
  integer, intent(in) :: outval !< Outgoing value
  character*(*), intent(in) :: delstr !< Delta label
  integer, intent(in) :: startval !< Start value for budget
  integer, intent(in) :: endval !< End value for budget
  ! Local variables
  write(*,200) budgetstr//' budget:', &
                      instr//' in',inval, &
                      outstr//' out',outval, &
                      'Delta '//delstr,inval-outval, &
                      'error',((endval-startval)-(inval-outval))
  200 format("diamonds: ",a19,10(a18,"=",i14,x,:,","))
end subroutine report_ibudget

subroutine create_gridded_particles_fields(parts)
! Arguments
type(particles), pointer :: parts !< Container for all types and memory
! Local variables
type(particles_gridded), pointer :: grd



  ! For convenience
  grd=>parts%grd

  !Finding the spread fields
  if ((grd%id_spread_uvel>0)  .or. (parts%pass_fields_to_ocean_model)) then
    grd%spread_uvel(:,:)=0.
    call sum_up_spread_fields(parts, grd%spread_uvel(grd%isc:grd%iec,grd%jsc:grd%jec), 'Uvel')
  endif
  if ( (grd%id_spread_vvel>0)  .or. (parts%pass_fields_to_ocean_model)) then
    grd%spread_vvel(:,:)=0.
    call sum_up_spread_fields(parts, grd%spread_vvel(grd%isc:grd%iec,grd%jsc:grd%jec), 'Vvel')
  endif

end subroutine create_gridded_particles_fields

end module particles_mod
