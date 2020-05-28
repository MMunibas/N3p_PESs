module surface
implicit none
real*8, parameter :: dk26f1 = 1.0d0/14.0d0, dk26f2 = 1.0d0/18.0d0, &
dk24f1 = 2.0d0/15.0d0, dk24f2 = 2.0d0/21.0d0, dk25f1 = 2.0d0/21.0d0,&
dk25f2 = 1.0d0/14.0d0, akf1=2.0d0/3.0d0
real*8, allocatable, dimension(:,:) :: asy_array1
integer :: na1

contains

function drker24(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker24, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker24 = dk24f1/xl**5 - dk24f2*xs/xl**6

end function drker24

function drker25(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker25, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker25 = dk25f1/xl**6 - dk25f2*xs/xl**7

end function drker25

function ddrker25(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker25, xl, xs

if (x .lt. xi) then
  ddrker25 = -dk25f2/xi**7
else
  ddrker25 = -6.0d0*dk25f1/x**7 + 7.0d0*dk25f2*xi/x**8
end if

end function ddrker25

function drker26(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker26, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker26 = dk26f1/xl**7 - dk26f2*xs/xl**8

end function drker26

function atker23(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: atker23, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

atker23 = 1.0d0 + xs*xl + 2.0d0*xs**2*xl  - akf1*xs**3

end function atker23

subroutine calcener(capr,smlr,theta, ener)
use RKHS            ! This module needs to be used by your code
implicit none
real*8 :: lambda
real*8, intent(out) :: ener
real*8, intent(in) :: capr, smlr, theta
real*8,parameter :: pi = acos(-1.0d0), piby180 = pi/180.0d0
real*8 :: asener, anener
real*8, dimension(:) :: x(3)
integer :: kk, ii
type(kernel), save  :: pes1           ! The kernel type is needed to set up and evaluate a RKHS model
logical, save :: stored = .false., kread = .false.
logical, save :: ker1 = .false.

if (.not. ker1) then
  inquire(file="pes1.kernel", exist=ker1)   ! file_exists will be true if the file exists and false otherwise
end if

if (.not. ker1) print*,"pes1.kernel does not exist"

lambda=0.1d-19

if (.not. stored ) then

open(unit=1001,file="asymp.dat", status = "old")

read(1001,*)na1
allocate(asy_array1(na1,2))
do ii = 1, na1
  read(1001,*)asy_array1(ii,1), asy_array1(ii,2)
end do

stored=.true.

end if

if (.not. kread) then
    call pes1%load_from_file("pes1.kernel")
    kread = .true.
end if

x(1)=(1.0d0-cos(theta))/2.0d0
x(2)=capr
x(3)=smlr

  asener = 0.0d0
  do kk = 1,na1
    asener = asener + drker26(smlr,asy_array1(kk,1))*asy_array1(kk,2)
  end do

  anener = 0.0d0
  call pes1%evaluate_fast(x,anener)

  ener = anener+asener

return

end subroutine calcener

end module

!=============================================================================
!potv is a subroutine for N3+ (N+ + N2) potential energy surface
!          N
!         /
!        /
!    R  /
!      /theta
!N____/____N
!     r
!v=total energy in hartree (output)
!smlr (r) = N2 diatomic distance in bohr (input)
!capr (R) = distance from center of mass of N2 to third N in bohr (input)
!ct = cos(theta), theta = angle between r and R (input)
!r range  covers 1.5 to 3.2 bohr)
!R range covers 1.6 to infinite bohr)
!theta range covers 0 to 180 degree
!=============================================================================

subroutine potv(v, smlr, capr, ct)
use surface
implicit none
real*8, dimension ( : ) :: r ( 3 )
real*8 :: capr, smlr, ct, theta, bigr, tinyr, rn1n2, rn2n3, rn3n1, intangle, v
real*8, parameter :: pi=acos(-1.0d0)

theta=acos(ct)
if (theta < pi/2.0d0)  theta=pi-theta

!Symmetry
rn1n2=smlr
rn2n3=sqrt((smlr/2.0d0)**2+capr**2-capr*smlr*cos(pi-theta))
rn3n1=sqrt((smlr/2.0d0)**2+capr**2-capr*smlr*cos(theta))

if ( rn1n2 >rn2n3 ) then
  rn1n2=rn2n3
  rn2n3=smlr
end if

tinyr=rn1n2
intangle=(rn1n2**2+rn2n3**2-rn3n1**2)/2.0d0/rn1n2/rn2n3
intangle=min(1.0d0,max(-1.0d0,intangle))
intangle=acos(intangle)

bigr=sqrt((tinyr/2.0d0)**2+rn2n3**2.0-tinyr*rn2n3*cos(intangle))
theta=(bigr**2+(tinyr/2.0d0)**2-rn2n3**2)/bigr/tinyr
theta=min(1.0d0,max(-1.0d0,theta))
theta=acos(theta)
!

r(1)=rn1n2
r(2)=rn2n3
r(3)=rn3n1

if (theta < pi/2.0d0)  theta=pi-theta

if  ( minval ( r ) < 1.1d0 .or. sum(r) < 4.5d0 ) then
  v = 0.25d0
else
  call calcener(bigr,tinyr,theta,v)
  if(v/=v)print*,"nan"
  if ( v > 0.25d0 ) then
    v = 0.25d0
  end if
end if
end subroutine
