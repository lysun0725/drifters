module particles_io

use constants_mod, only: pi, omega, HLF

use mpp_domains_mod, only: domain2D
use mpp_domains_mod, only: mpp_domain_is_tile_root_pe,mpp_get_domain_tile_root_pe
use mpp_domains_mod, only: mpp_get_tile_pelist,mpp_get_tile_npes,mpp_get_io_domain,mpp_get_tile_id

use mpp_mod, only: mpp_npes, mpp_pe, mpp_root_pe, mpp_sum, mpp_min, mpp_max, NULL_PE
use mpp_mod, only: mpp_send, mpp_recv, mpp_gather, mpp_chksum
use mpp_mod, only: COMM_TAG_11, COMM_TAG_12, COMM_TAG_13, COMM_TAG_14

use fms_mod, only: stdlog, stderr, error_mesg, FATAL, WARNING, NOTE
use fms_mod, only: field_exist, file_exist, read_data, write_data

use fms_io_mod, only: get_instance_filename
use fms_io_mod, only : save_restart, restart_file_type, free_restart_type, set_meta_global
use fms_io_mod, only : register_restart_axis, register_restart_field, set_domain, nullify_domain
use fms_io_mod, only : read_unlimited_axis =>read_compressed, field_exist, get_field_size

use mpp_mod,    only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
use mpp_mod,    only : CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_LOOP
use fms_mod,    only : clock_flag_default

use time_manager_mod, only: time_type, get_date, get_time, set_date, operator(-)

use MOM_grid, only : ocean_grid_type
use MOM, only      : MOM_control_struct

use particles_framework, only: particles_gridded, xyt, particle, particles, buffer
use particles_framework, only: pack_traj_into_buffer2,unpack_traj_from_buffer2
use particles_framework, only: find_cell,find_cell_by_search,count_parts,is_point_in_cell,pos_within_cell,append_posn
!use particles_framework, only: count_bonds, form_a_bond
use particles_framework, only: find_individual_particle
use particles_framework, only: push_posn
use particles_framework, only: add_new_part_to_list,destroy_particle
use particles_framework, only: increase_ibuffer,grd_chksum2,grd_chksum3
use particles_framework, only: sum_mass,sum_heat,bilin
!params !Niki: write a subroutine to get these
use particles_framework, only: nclasses, buffer_width, buffer_width_traj
use particles_framework, only: verbose, really_debug, debug, restart_input_dir,make_calving_reproduce
use particles_framework, only: ignore_ij_restart, use_slow_find,generate_test_particles!,print_part
use particles_framework, only: force_all_pes_traj
!use particles_framework, only: check_for_duplicates_in_parallel
!use particles_framework, only: split_id, id_from_2_ints, generate_id

implicit none ; private

include 'netcdf.inc'

public particles_io_init
public read_restart_parts, write_restart,write_trajectory
public interp_flds

!Local Vars
integer, parameter :: file_format_major_version=0
integer, parameter :: file_format_minor_version=1
!I/O vars
type(domain2d), pointer, save :: io_domain=>NULL()
integer, save :: io_tile_id(1), io_tile_root_pe, io_npes
integer, allocatable,save :: io_tile_pelist(:)
logical :: is_io_tile_root_pe = .true.

integer :: clock_trw,clock_trp

#ifdef _FILE_VERSION
  character(len=128) :: version = _FILE_VERSION
#else
  character(len=128) :: version = 'unknown'
#endif

contains

!> Initialize parallel i/o
subroutine particles_io_init(parts, io_layout)
type(particles), pointer :: parts !< particles container
integer, intent(in) :: io_layout(2) !< Decomposition of i/o processors

integer :: np
integer :: stdlogunit, stderrunit

  ! Get the stderr and stdlog unit numbers
  stderrunit=stderr()
  stdlogunit=stdlog()
  write(stdlogunit,*) "particles_framework: "//trim(version)

  !I/O layout init
  io_tile_id=-1
  io_domain => mpp_get_io_domain(parts%grd%domain)
  if(associated(io_domain)) then
     io_tile_id = mpp_get_tile_id(io_domain)
     is_io_tile_root_pe = mpp_domain_is_tile_root_pe(io_domain)
     io_tile_root_pe = mpp_get_domain_tile_root_pe(io_domain)
     np=mpp_get_tile_npes(io_domain)
     allocate(io_tile_pelist(np))
     call mpp_get_tile_pelist(io_domain,io_tile_pelist)
     io_npes = io_layout(1)*io_layout(2)
  endif

  clock_trw=mpp_clock_id( 'particles-traj write', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )
  clock_trp=mpp_clock_id( 'particles-traj prepare', flags=clock_flag_default, grain=CLOCK_SUBCOMPONENT )

end subroutine particles_io_init

! ##############################################################################

!> Write an particle restart file
subroutine write_restart(parts)
! Arguments
type(particles), pointer :: parts !< particles container
! Local variables
!type(bond), pointer :: current_bond
integer :: i,j,id
character(len=35) :: filename
character(len=35) :: filename_bonds
type(particle), pointer :: this=>NULL()
integer :: stderrunit
!I/O vars
type(restart_file_type) :: parts_restart
type(restart_file_type) :: parts_bond_restart
integer :: nparts, nbonds
integer :: n_static_parts
logical :: check_bond_quality
type(particles_gridded), pointer :: grd
real, allocatable, dimension(:) :: lon,          &
                                   lat,          &
                                   uvel,         &
                                   vvel,         &
                                   mass,         &
                                   axn,          &
                                   ayn,          &
                                   bxn,          &
                                   byn,          &
                                   thickness,    &
                                   width,        &
                                   length,       &
                                   start_lon,    &
                                   start_lat,    &
                                   start_day,    &
                                   start_mass,   &
                                   mass_scaling, &
                                   mass_of_bits, &
                                   static_part,  &
                                   heat_density

integer, allocatable, dimension(:) :: ine,              &
                                      jne,              &
                                      id_cnt,           &
                                      id_ij,            &
                                      start_year,       &
                                      first_id_cnt,     &
                                      other_id_cnt,     &
                                      first_id_ij,      &
                                      other_id_ij,      &
                                      first_part_jne,         &
                                      first_part_ine,         &
                                      other_part_jne,         &
                                      other_part_ine


integer :: grdi, grdj


 
end subroutine write_restart

! ##############################################################################

subroutine read_restart_parts(parts,Time, MOM_CS)
! Arguments
type(particles), pointer :: parts
type(time_type), intent(in) :: Time
type(MOM_control_struct), pointer, intent(in) :: MOM_CS

!Local variables
integer :: k, siz(4), nparts_in_file, nparts_read
logical :: lres, found_restart, found!, replace_particle_num
logical :: explain
logical :: multiPErestart  ! Not needed with new restart read; currently kept for compatibility
real :: lon0, lon1, lat0, lat1
real :: pos_is_good, pos_is_good_all_pe
character(len=33) :: filename, filename_base
type(particles_gridded), pointer :: grd

type(particle) :: localpart
integer :: stderrunit, i, j, cnt, ij

real, allocatable,dimension(:) :: lon,	&
                                  lat,	&
                                  depth, &
                                  id	

  ! Get the stderr unit number
  stderrunit=stderr()

  ! For convenience
  grd=>parts%grd

  ! Zero out nparts_in_file
  nparts_in_file = 0

  filename_base=trim(restart_input_dir)//'drifters.res.nc'

  found_restart = find_restart_file(filename_base, filename, multiPErestart, io_tile_id(1))

  if (found_restart) then
    filename = filename_base
    call get_field_size(filename,'i',siz, field_found=found, domain=grd%domain) 
    nparts_in_file = siz(1)
    print *,'NPARTS= ',nparts_in_file
    allocate(lon(nparts_in_file))
    allocate(lat(nparts_in_file))
    allocate(depth(nparts_in_file))
    allocate(id(nparts_in_file))

    call read_unlimited_axis(filename,'lon',lon,domain=grd%domain)
    call read_unlimited_axis(filename,'lat',lat,domain=grd%domain)
    call read_unlimited_axis(filename,'depth',depth,domain=grd%domain)
    call read_unlimited_axis(filename,'drifter_num',id,domain=grd%domain)
  end if ! found_restart ln 569

  ! Find approx outer bounds for tile
  lon0=minval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lon1=maxval( grd%lon(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat0=minval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )
  lat1=maxval( grd%lat(grd%isc-1:grd%iec,grd%jsc-1:grd%jec) )

  do k=1, nparts_in_file
    localpart%lon=lon(k)
    localpart%lat=lat(k)

    if (use_slow_find) then
      lres=find_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
    else
      lres=find_cell_by_search(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne)
    endif

    if (really_debug) then
      write(stderrunit,'(a,i8,a,2f9.4,a,i8)') 'diamonds, read_restart_part: part ',k,' is at ',localpart%lon,localpart%lat,&
           & ' on PE ',mpp_pe()
      write(stderrunit,*) 'diamonds, read_restart_parts: lres = ',lres
    endif

    if (lres) then ! True if the particle resides on the current processors computational grid
      lres=pos_within_cell(grd, localpart%lon, localpart%lat, localpart%ine, localpart%jne, localpart%xi, localpart%yj)
      call interp_flds(grd,localpart%ine,localpart%jne,localpart%xi,localpart%yj,localpart%uvel, localpart%vvel)
      call add_new_part_to_list(parts%list(localpart%ine,localpart%jne)%first, localpart)
    endif
  end do ! ln 650

  if (found_restart) then
    deallocate(lon,          &
               lat,          &
               depth,        &
               id )
  end if


end subroutine read_restart_parts

! ##############################################################################
!> Write a trajectory-based diagnostics file
subroutine write_trajectory(trajectory, save_short_traj)
! Arguments
type(xyt), pointer :: trajectory !< An iceberg trajectory
logical, intent(in) :: save_short_traj !< If true, record less data
! Local variables
integer :: iret, ncid, i_dim, i
integer :: lonid, latid, yearid, dayid, uvelid, vvelid, idcntid, idijid
integer :: uoid, void, uiid, viid, uaid, vaid, sshxid, sshyid, sstid, sssid
integer :: cnid, hiid
integer :: mid, did, wid, lid, mbid, hdid
character(len=37) :: filename
character(len=7) :: pe_name
type(xyt), pointer :: this, next
integer :: stderrunit, cnt, ij
!I/O vars
type(xyt), pointer :: traj4io=>null()
integer :: ntrajs_sent_io,ntrajs_rcvd_io
integer :: from_pe,np
type(buffer), pointer :: obuffer_io=>null(), ibuffer_io=>null()
logical :: io_is_in_append_mode


  
end subroutine write_trajectory


! ##############################################################################

integer function inq_var(ncid, var, unsafe)
! Arguments
integer, intent(in) :: ncid
character(len=*), intent(in) :: var
logical, optional, intent(in) :: unsafe
! Local variables
integer :: iret
integer :: stderrunit
logical :: unsafely=.false.

if(present(unsafe)) unsafely=unsafe
  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_inq_varid(ncid, var, inq_var)
  if (iret .ne. NF_NOERR) then
    if (.not. unsafely) then
      write(stderrunit,*) 'diamonds, inq_var: nf_inq_varid ',var,' failed'
      call error_mesg('diamonds, inq_var', 'netcdf function returned a failure!', FATAL)
    else
      inq_var=-1
    endif
  endif

end function inq_var

! ##############################################################################

integer function def_var(ncid, var, ntype, idim)
! Arguments
integer, intent(in) :: ncid, ntype, idim
character(len=*), intent(in) :: var
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_def_var(ncid, var, ntype, 1, idim, def_var)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, def_var: nf_def_var failed for ',trim(var)
    call error_mesg('diamonds, def_var', 'netcdf function returned a failure!', FATAL)
  endif

end function def_var

! ##############################################################################

integer function inq_varid(ncid, var)
! Arguments
integer, intent(in) :: ncid
character(len=*), intent(in) :: var
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_inq_varid(ncid, var, inq_varid)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, inq_varid: nf_inq_varid failed for ',trim(var)
    call error_mesg('diamonds, inq_varid', 'netcdf function returned a failure!', FATAL)
  endif

end function inq_varid

! ##############################################################################

subroutine put_att(ncid, id, att, attval)
! Arguments
integer, intent(in) :: ncid, id
character(len=*), intent(in) :: att, attval
! Local variables
integer :: vallen, iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  vallen=len_trim(attval)
  iret = nf_put_att_text(ncid, id, att, vallen, attval)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, put_att: nf_put_att_text failed adding', &
      trim(att),' = ',trim(attval)
    call error_mesg('diamonds, put_att', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_att

! ##############################################################################

real function get_double(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_double(ncid, id, i, get_double)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, get_double: nf_get_var1_double failed reading'
    call error_mesg('diamonds, get_double', 'netcdf function returned a failure!', FATAL)
  endif

end function get_double

! ##############################################################################

integer function get_int(ncid, id, i)
! Arguments
integer, intent(in) :: ncid, id, i
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret=nf_get_var1_int(ncid, id, i, get_int)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, get_int: nf_get_var1_int failed reading'
    call error_mesg('diamonds, get_int', 'netcdf function returned a failure!', FATAL)
  endif

end function get_int

! ##############################################################################

subroutine put_double(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid, id, i
real, intent(in) :: val
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_double(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, put_double: nf_put_vara_double failed writing'
    call error_mesg('diamonds, put_double', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_double

! ##############################################################################

subroutine put_int(ncid, id, i, val)
! Arguments
integer, intent(in) :: ncid, id, i, val
! Local variables
integer :: iret
integer :: stderrunit

  ! Get the stderr unit number
  stderrunit=stderr()

  iret = nf_put_vara_int(ncid, id, i, 1, val)
  if (iret .ne. NF_NOERR) then
    write(stderrunit,*) 'diamonds, put_int: nf_put_vara_int failed writing'
    call error_mesg('diamonds, put_int', 'netcdf function returned a failure!', FATAL)
  endif

end subroutine put_int


! ##############################################################################

logical function find_restart_file(filename, actual_file, multiPErestart, tile_id)
  character(len=*), intent(in) :: filename
  character(len=*), intent(out) :: actual_file
  logical, intent(out) :: multiPErestart
  integer, intent(in) :: tile_id

  character(len=6) :: pe_name

  find_restart_file = .false.

  ! If running as ensemble, add the ensemble id string to the filename
  call get_instance_filename(filename, actual_file)
    
  ! Prefer combined restart files.
  inquire(file=actual_file,exist=find_restart_file)
  if (find_restart_file) return
    
  ! Uncombined restart
  if(tile_id .ge. 0) then
    write(actual_file,'(A,".",I4.4)') trim(actual_file), tile_id
  else
  if (mpp_npes()>10000) then
     write(pe_name,'(a,i6.6)' )'.', mpp_pe()    
  else
     write(pe_name,'(a,i4.4)' )'.', mpp_pe()    
  endif
  actual_file=trim(actual_file)//trim(pe_name)
  endif
  inquire(file=actual_file,exist=find_restart_file)
  if (find_restart_file) then
     multiPErestart=.true.
     return
  endif

  ! No file found, Reset all return parameters
  find_restart_file=.false.
  actual_file = ''
  multiPErestart=.false.

end function find_restart_file

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

end module
