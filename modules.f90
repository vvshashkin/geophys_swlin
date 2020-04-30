module const
implicit none
real(8), parameter :: pi = acos(-1._8)
real(8), parameter :: radz = 0.6371229E7_8
real(8), parameter :: grav = 9.80616_8
real(8) href
real(8) omega
logical lbeta

contains
subroutine init_const
namelist /cst/ href, omega, lbeta
implicit none

HREF = 3000._8
omega= 0._8
lbeta = .false.

open(13, file="namelist", form = "formatted")
read(13, cst)
close(13)

end subroutine init_const

end module const

module prmt
implicit none
integer NLAT,NLON
integer NSTEP
integer NZAP
real(8) dt
real(8) rdfi, rdlam

contains
subroutine init_prmt
use const
namelist /prm/ NLON,NLAT,NSTEP,NZAP,dt
implicit none

NLAT = 120
NLON = 240
dt = 1800._8
NSTEP = 96
NZAP = 2

open(13, file="namelist", form = "formatted")
read(13, prm)
close(13)

rdfi = pi/NLAT
rdlam = 2._8*pi/NLON
end subroutine init_prmt
end module prmt

module gem
implicit none
real(8), allocatable :: ds(:), rf(:)
real(8), allocatable :: fcoriu(:), fcoriv(:)

contains
subroutine init_gem
use prmt, only: NLON,NLAT, rdfi, rdlam
use const, only: pi, omega, lbeta
integer j
allocate(ds(1:NLAT))
allocate(rf(0:NLAT))
allocate(fcoriu(1:NLAT))
allocate(fcoriv(0:NLAT))

do j=1,NLAT
  ds(j) = cos(-.5_8*pi+(j-.5_8)*rdfi)
  if(lbeta) then
    fcoriu(j) = 2._8*omega*sin(-.5_8*pi+(j-.5_8)*rdfi)
  else
    fcoriu(j) = 2._8*omega
  endif
end do

do j=0,NLAT
  rf(j) = cos(-.5_8*pi+j*rdfi)
  if(lbeta) then
    fcoriv(j) = 2._8*omega*sin(-.5_8*pi+j*rdfi)
  else
    fcoriv(j) = 2._8*omega
  endif
end do

end subroutine init_gem
end module gem

module fld
implicit none
real(8), allocatable :: h(:,:)
real(8), allocatable :: u(:,:)
real(8), allocatable :: v(:,:)

contains
subroutine init_fld
use prmt
allocate(h(0:NLON,NLAT))
allocate(u(0:NLON,NLAT))
allocate(v(0:NLON,0:NLAT))

call initial_cond(h, v, u)

end subroutine init_fld

subroutine initial_cond(h, u, v)
use const
use prmt
namelist /ini/ phi0, lam0, rhill, hill_type
implicit none
real(8) h(0:NLON,NLAT)
real(8) u(0:NLON,NLAT)
real(8) v(0:NLON,0:NLAT)
!test parameters:
real(8) phi0, lam0
real(8), parameter :: hmax = 1._8
real(8) rhill
character(256) :: hill_type
!local vars
integer i, j
real(8) phi, lam
real(8) zx, zy, zz, zx0, zy0, zz0, zr

phi0 = 45._8; lam0 = 90._8
rhill = 1000.e3_8
hill_type="gauss"

open(13, file="namelist", form = "formatted")
read(13, ini)
close(13)
lam0 = lam0/180.*pi
phi0 = phi0/180.*pi

u = 0._8
v = 0._8

zx0 = cos(phi0)*cos(lam0)
zy0 = cos(phi0)*sin(lam0)
zz0 = sin(phi0)
do j = 1, NLAT
  phi = -.5_8*pi+(j-.5_8)*rdfi
  do i=0,NLON-1
    lam = i*rdlam
    zx = cos(phi)*cos(lam)
    zy = cos(phi)*sin(lam)
    zz = sin(phi)
    if(trim(hill_type)=="gauss") then
      h(i,j) = hmax*exp(-(radz*acos(zx0*zx+zy0*zy+zz0*zz)/rhill)**2)
    else if(trim(hill_type)=="cone") then
      h(i,j) = hmax*max(0._8,1._8-radz*acos(zx0*zx+zy0*zy+zz0*zz)/rhill)
    else if(trim(hill_type)=="cylinder") then
      if(acos(zx0*zx+zy0*zy+zz0*zz)<=rhill/radz) then
        h(i,j) = hmax
      else
        h(i,j) = 0._8
      end if
    else if(trim(hill_type)=="cosbell") then
      zr = radz*acos(zx0*zx+zy0*zy+zz0*zz)
      if(zr<=rhill) then
        h(i,j) = hmax*.5_8*(1._8+cos(pi*zr/rhill))
      else
        h(i,j) = 0._8
      end if
    else
      print *, "unknown type of initial height distribution ", trim(hill_type)
    end if
  end do
end do

end subroutine initial_cond
end module fld
