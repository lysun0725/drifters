program main



use fms_mod, only: read_data


real, dimension(:,:) : u,v


!allocatate u,v
call read_data(filename,varname,u)
call read_data(filename,varname,v)

call particles_init(...)
call particles_init(...)
call particles_run(...)

end program main
