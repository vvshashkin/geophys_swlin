subroutine step(h,u,v,dt)
use prmt, only:NLON,NLAT
use const
implicit none
real(8) h(0:NLON,NLAT)
real(8) u(0:NLON,NLAT)
real(8) v(0:NLON,0:NLAT)
real(8) dt
integer iter
real(8) zeps0
real(8) zeps
real(8) h1(0:NLON,NLAT)
real(8) u1(0:NLON,NLAT)
real(8) v1(0:NLON,0:NLAT)

zeps0 = max(sqrt(sum(u**2)+sum(v**2)+sum(h**2)),1.e-12_8)
zeps = 1._8

u(0,:) = u(NLON,:)
h(NLON,:) = h(0,:)
v(:,NLAT) = 0._8; v(:,0) = 0._8
v(NLON,:) = v(0,:)

h1 = h
u1 = u
v1 = v
iter = 1

do while(zeps>1.e-12_8.and.iter<=1000)
  call model_op(h1,u1,v1)
  h1 = h1*dt/dble(iter)
  u1 = u1*dt/dble(iter)
  v1 = v1*dt/dble(iter)
  
  h = h+h1
  u = u+u1
  v = v+v1

  zeps = sqrt(sum(u1**2)+sum(v1**2)+sum(h1**2))/zeps0
  !print *, iter, zeps
  iter = iter+1
end do

u(0,:) = u(NLON,:)
h(NLON,:) = h(0,:)
v(:,NLAT) = 0._8; v(:,0) = 0._8
v(NLON,:) = v(0,:)

print *, "exp solver:", iter, " iterations, residual ~", zeps!, zeps0

end subroutine step

subroutine model_op(h, u, v)
use prmt, only:NLON,NLAT,rdlam,rdfi
use const
use gem
implicit none
real(8) h(0:NLON,NLAT)
real(8) u(0:NLON,NLAT)
real(8) v(0:NLON,0:NLAT)
real(8) dt
integer i, j
real(8) zeps0
real(8) zeps
real(8) h1(0:NLON,NLAT)
real(8) u1(0:NLON,NLAT)
real(8) v1(0:NLON,0:NLAT)
real(8) zdiv1


do j=1, 20
  call polar_filter(h(0:NLON-1,j),21-j)
  call polar_filter(u(0:NLON-1,j),21-j)
  call polar_filter(v(0:NLON-1,j),21-j)
  call polar_filter(h(0:NLON-1,NLAT-j+1),21-j)
  call polar_filter(u(0:NLON-1,NLAT-j+1),21-j)
  call polar_filter(u(0:NLON-1,NLAT-j),21-j)
end do


h1 = h
u1 = u
v1 = v

u1(0,:) = u1(NLON,:)
h1(NLON,:) = h1(0,:)
v1(:,NLAT) = 0._8; v1(:,0) = 0._8
v1(NLON,:) = v1(0,:)

do j=1, NLAT
  do i=0,NLON-1
    zdiv1  = (u1(i+1,j)-u1(i,j))/(rdlam*RADZ*ds(j))+(v1(i,j)*rf(j)-v1(i,j-1)*rf(j-1))/(rdfi*RADZ*ds(j))
    h(i,j) =-HREF*zdiv1
  end do
end do
v(:,0) = 0._8
do j=1, NLAT-1
  do i=0,NLON-1
    v(i,j) =-grav*(h1(i,j+1)-h1(i,j))/(rdfi*RADZ)-.25_8*(fcoriu(j+1)*(u1(i+1,j+1)+u1(i,j+1))+fcoriu(j)*(u1(i+1,j)+u1(i,j)))
  end do
end do
v(:,NLAT) = 0._8
!v-boundary conditions for coriollis
v1(0:NLON/2-1,0) = .5_8*(v1(0:NLON/2-1,1)-v1(NLON/2:NLON-1,1))
v1(NLON/2:NLON-1,0) = -v1(0:NLON/2-1,0)
v1(NLON,0) = v1(0,0)
v1(0:NLON/2-1,NLAT) = .5_8*(v1(0:NLON/2-1,NLAT-1)-v1(NLON/2:NLON-1,NLAT-1))
v1(NLON/2:NLON-1,NLAT) = -v1(0:NLON/2-1,NLAT)
v1(NLON,NLAT) = v1(0,NLAT)
do j=1, NLAT
  do i=1,NLON
    u(i,j) =-grav*(h1(i,j)-h1(i-1,j))/(rdlam*RADZ*ds(j))+.25_8*(fcoriv(j)*(v1(i,j)+v1(i-1,j))+fcoriv(j-1)*(v1(i,j-1)+v1(i-1,j-1)))
  end do
end do

u(0,:) = u(NLON,:)
h(NLON,:) = h(0,:)
v(:,NLAT) = 0._8; v(:,0) = 0._8
v(NLON,:) = v(0,:)
end subroutine model_op

subroutine polar_filter(f,niter)
use prmt, only: NLON
implicit none
real(8) f(1:NLON)
integer iter, niter
real(8) f1(0:NLON+1)

f1(1:NLON) = f(1:NLON)
do iter = 1, niter
  f1(0) = f1(NLON)
  f1(NLON+1) = f1(1)
  f1(1:NLON) = .25_8*(f1(0:NLON-1)+2._8*f1(1:NLON)+f1(2:NLON+1))    
end do
f(1:NLON) = f1(1:NLON)
end subroutine polar_filter
