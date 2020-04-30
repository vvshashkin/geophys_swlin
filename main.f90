!linear shallow-water equations solver on a sphere
!for practical excersises of geophys.hydrodyn. course
!V. Shashkin, INM RAS, March 2019

!Horizontal discretization C-grid, second order
!Time discretization: matrix exponent method

program swlin
use const
use prmt
use fld
use gem
implicit none
integer irec, it
character(:), allocatable :: namfname
integer nargs, arg_len

nargs = command_argument_count()
if(nargs >= 1) then
    call get_command_argument(1, length=arg_len)
    allocate(character(arg_len) :: namfname)
    call get_command_argument(1, namfname)
else
    namfname = "namelist"
end if
print *, "namelist: ", namfname

call init_const(namfname)
call init_prmt(namfname)
call init_fld(namfname)
call init_gem

open(117,file="f.dat", access="direct",recl=4*NLON*NLAT)
call wrfld(h,u,v, 1, 117)
irec=2

do it=1, NSTEP
  call step(h,u,v,dt)
  if(mod(it,NZAP)==0) then
    call wrfld(h,u,v,irec,117)
    print *, "fields record ",irec, "written"
    print *, "h_max=", maxval(h), "h_min=", minval(h)
    irec = irec+1
  end if
end do


close(117)
end program swlin

subroutine wrfld(h,u,v,irec,ifd)
use prmt
implicit none
real(8) h(0:NLON,NLAT)
real(8) u(0:NLON,NLAT)
real(8) v(0:NLON,0:NLAT)
integer irec, ifd

write(ifd, rec=3*(irec-1)+1) real(h(0:NLON-1,1:NLAT),4)
u(0,:) = u(NLON,:)
write(ifd, rec=3*(irec-1)+2) real(.5_8*(u(0:NLON-1,1:NLAT)+u(1:NLON,1:NLAT)),4)
v(:,0) = v(:,1)
v(:,NLAT) = v(:,NLAT-1)
write(ifd, rec=3*(irec-1)+3) real(.5_8*(v(0:NLON-1,0:NLAT-1)+v(0:NLON-1,1:NLAT)),4)
end subroutine wrfld

