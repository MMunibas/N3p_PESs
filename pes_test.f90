program testpes
implicit none
real*8 :: capr, smlr, theta, ct, ener
real*8, parameter :: pi = acos(-1.0d0)

!==============================
!          N
!         /
!        /
!    R  /
!      /theta
!N____/____N
!     r
!==============================

smlr=2.3d0 !r in bohr
capr=4.5d0 !R in bohr
theta=175.0d0 ! theta in degree
theta=theta*pi/180.0d0
ct=cos(theta)

call potv(ener, smlr, capr, ct)

write(*,*)"Energy = ", ener, "Hartree"

end program
