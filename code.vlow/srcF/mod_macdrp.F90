module macdrp_mod

! This module contains the variables and subroutines
! used in the DRP/opt MacCormack fd operator
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-14 16:30:51 -0500 (Wed, 14 Jan 2009) $
! $Revision: 510 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"

use constants_mod, only : SEIS_GEO
use math_mod
use para_mod
use mpi
use mpi_mod
use media_mod
use grid_mod

implicit none
private
public ::                                      &
    Txx, Tyy, Txy, Vx, Vy, Tzz, Txz, Tyz, Vz,  &
    hTxx,hTyy,hTxy,hVx,hVy,hTzz,hTxz,hTyz,hVz, &
    macdrp_init,                               &
    macdrp_syn,                                &
    macdrp_mesg_init,                          &
    macdrp_destroy,                            &
    macdrp_LxF_LyF_LzF,                        &
    macdrp_LxB_LyB_LzB,                        &
    macdrp_LxB_LyB_LzF,                        &
    macdrp_LxF_LyF_LzB,                        &
    macdrp_LxB_LyF_LzF,                        &
    macdrp_LxF_LyB_LzB,                        &
    macdrp_LxF_LyB_LzF,                        &
    macdrp_LxB_LyF_LzB,                        &
    macdrp_RK_beg,                             &
    macdrp_RK_inn,                             &
    macdrp_RK_fin,                             &
    macdrp_check,                              &
    atten_graves

interface macdrp_LxF_LyF_LzF
  module procedure in_LxF_LyF_LzF
end interface
interface macdrp_LxB_LyB_LzB
  module procedure in_LxB_LyB_LzB
end interface
interface macdrp_LxF_LyF_LzB
  module procedure in_LxF_LyF_LzB
end interface
interface macdrp_LxB_LyB_LzF
  module procedure in_LxB_LyB_LzF
end interface
interface macdrp_LxB_LyF_LzF
  module procedure in_LxB_LyF_LzF
end interface
interface macdrp_LxF_LyB_LzB
  module procedure in_LxF_LyB_LzB
end interface
interface macdrp_LxF_LyB_LzF
  module procedure in_LxF_LyB_LzF
end interface
interface macdrp_LxB_LyF_LzB
  module procedure in_LxB_LyF_LzB
end interface

DEFFDWET
DEFFDWET24
DEFFDWET22
DEFLDDRK2A
DEFLDDRK2B
DEFLDDRK4A
DEFLDDRK4B
HOCWETL
HOCWETR

real(SP),dimension(:,:,:),allocatable ::        &
      Txx, Tyy, Txy, Vx, Vy, Tzz, Txz, Tyz, Vz, &
     hTxx,hTyy,hTxy,hVx,hVy,hTzz,hTxz,hTyz,hVz, &
     mTxx,mTyy,mTxy,mVx,mVy,mTzz,mTxz,mTyz,mVz, &
     tTxx,tTyy,tTxy,tVx,tVy,tTzz,tTxz,tTyz,tVz
real(SP),dimension(:,:,:,:),allocatable,public :: &
     matVx2Vz,matVy2Vz,matF2Vz
real(SP),dimension(:,:),allocatable,public :: &
     TxSrc,TySrc,TzSrc,                       &
     VxSrc,VySrc,VzSrc
real(SP),dimension(4),public :: firRKa,firRKb, secRKa,secRKb
integer,dimension(SEIS_GEO*2,SEIS_GEO*2+1),public :: indx
integer ierr
integer,dimension(MPI_STATUS_SIZE) :: istatus
integer,dimension(36) :: reqXB, reqXF, reqYB, reqYF, reqZB, reqZF
integer,dimension(MPI_STATUS_SIZE,36) :: reqstat
#ifdef VERBOSE
integer fid_out
#endif

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine macdrp_init
integer ierr
allocate( Txx(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Txx=0.0_SP
allocate( Tyy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Tyy=0.0_SP
allocate( Txy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Txy=0.0_SP
allocate( Vx (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Vx =0.0_SP
allocate( Vy (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Vy =0.0_SP
allocate( Tzz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Tzz=0.0_SP
allocate( Txz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Txz=0.0_SP
allocate( Tyz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Tyz=0.0_SP
allocate( Vz (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Vz =0.0_SP
allocate(hTxx(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTxx=0.0_SP
allocate(hTyy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTyy=0.0_SP
allocate(hTxy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTxy=0.0_SP
allocate(hVx (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hVx =0.0_SP
allocate(hVy (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hVy =0.0_SP
allocate(hTzz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTzz=0.0_SP
allocate(hTxz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTxz=0.0_SP
allocate(hTyz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTyz=0.0_SP
allocate(hVz (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hVz =0.0_SP
allocate(mTxx(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTxx=0.0_SP
allocate(mTyy(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTyy=0.0_SP
allocate(mTxy(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTxy=0.0_SP
allocate(mVx (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mVx =0.0_SP
allocate(mVy (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mVy =0.0_SP
allocate(mTzz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTzz=0.0_SP
allocate(mTxz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTxz=0.0_SP
allocate(mTyz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTyz=0.0_SP
allocate(mVz (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mVz =0.0_SP
allocate(tTxx(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTxx=0.0_SP
allocate(tTyy(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTyy=0.0_SP
allocate(tTxy(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTxy=0.0_SP
allocate(tVx (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tVx =0.0_SP
allocate(tVy (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tVy =0.0_SP
allocate(tTzz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTzz=0.0_SP
allocate(tTxz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTxz=0.0_SP
allocate(tTyz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTyz=0.0_SP
allocate(tVz (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tVz =0.0_SP
if (ierr>0) then
   print *, "can't allocate variable in macdrp_init"
   stop 1
end if
allocate(matVx2Vz(SEIS_GEO,SEIS_GEO,nx,ny),stat=ierr); matVx2Vz=0.0_SP
allocate(matVy2Vz(SEIS_GEO,SEIS_GEO,nx,ny),stat=ierr); matVy2Vz=0.0_SP
allocate(matF2Vz(SEIS_GEO,SEIS_GEO,nx,ny),stat=ierr); matF2Vz=0.0_SP
allocate(TxSrc(nx1:nx2,ny1:ny2),stat=ierr); TxSrc=0.0_SP
allocate(TySrc(nx1:nx2,ny1:ny2),stat=ierr); TySrc=0.0_SP
allocate(TzSrc(nx1:nx2,ny1:ny2),stat=ierr); TzSrc=0.0_SP
allocate(VxSrc(nx1:nx2,ny1:ny2),stat=ierr); VxSrc=0.0_SP
allocate(VySrc(nx1:nx2,ny1:ny2),stat=ierr); VySrc=0.0_SP
allocate(VzSrc(nx1:nx2,ny1:ny2),stat=ierr); VzSrc=0.0_SP
! main
indx(:,SEIS_GEO*2+1)=(/ ni1+LenFD,ni2-LenFD, &
                        nj1+LenFD,nj2-LenFD, &
                        nk1+LenFD,nk2-LenFD /)
indx(:,SEIS_GEO*2  )=(/ ni1,ni2,nj1,nj2,nk2-LenFD+1,nk2 /) ! z2
indx(:,SEIS_GEO*2-1)=(/ ni1,ni2,nj1,nj2,nk1,nk1+LenFD-1 /) ! z1
indx(:,1)=(/ ni1,ni2,nj1,nj1+LenFD-1,nk1+LenFD,nk2-LenFD /) ! y1
indx(:,2)=(/ ni1,ni2,nj2-LenFD+1,nj2,nk1+LenFD,nk2-LenFD /) ! y2
indx(:,3)=(/ ni1,ni1+LenFD-1,nj1+LenFD,nj2-LenFD,nk1+LenFD,nk2-LenFD /) ! x1
indx(:,4)=(/ ni2-LenFD+1,ni2,nj1+LenFD,nj2-LenFD,nk1+LenFD,nk2-LenFD /) ! x2
! rk coefficient
firRKa=(/ RK4a2, RK4a3, RK4a4, 0.0 /)
firRKb=(/ RK4b1, RK4b2, RK4b3, RK4b4 /)
!secRKa=(/ RK4a2, RK4a3, RK4a4, 0.0 /)
!secRKb=(/ RK4b1, RK4b2, RK4b3, RK4b4 /)
secRKa=(/ RK2a2, 0.0, 0.0, 0.0 /)
secRKb=(/ RK2b1, RK2b2, 0.0, 0.0 /)
! mat to convert V,z
call coef_fdxy2fdz
#ifdef VERBOSE
  fid_out=9050
  open(fid_out,                                                                      &
       file='log_maxval_'//trim(set_mpi_subfix(thisid(1),thisid(2),thisid(3)))//'.dat', &
       status='unknown')
#endif
end subroutine macdrp_init
subroutine macdrp_destroy
deallocate( Txx, Tyy, Txy, Vx, Vy)
deallocate(hTxx,hTyy,hTxy,hVx,hVy)
deallocate(mTxx,mTyy,mTxy,mVx,mVy)
deallocate(tTxx,tTyy,tTxy,tVx,tVy)
#ifdef VERBOSE
  close(fid_out)
#endif
end subroutine macdrp_destroy
subroutine macdrp_check(ntime)
integer,intent(in) :: ntime
real(SP) :: V1,V2,V3,T11,T22,T33,T12,T13,T23,W
integer ierr
#ifndef CheckOverFlow
 return
#endif
if (mod(ntime,1)==0) then
    V1=maxval(abs(Vx))
    V2=maxval(abs(Vy))
    V3=maxval(abs(Vz))
   T11=maxval(abs(Txx))
   T22=maxval(abs(Tyy))
   T33=maxval(abs(Tzz))
   T12=maxval(abs(Txy))
   T13=maxval(abs(Txz))
   T23=maxval(abs(Tyz))
#ifdef VERBOSE
   write(fid_out,'(i5,9es12.5)') ntime, V1,V2,V3,T11,T22,T33,T12,T13,T23
#endif
   W=max(V1,V2,V3,T11,T22,T33,T12,T13,T23)
   if (W>=huge(1.0)) then
      print *, "Overflow error: "
      write(*,"(i5,i3.2,2(i2.2),9(es12.5,3i5))") ntime,thisid(1),thisid(2),thisid(3), &
         V1, maxloc(abs(Vx)), V2, maxloc(abs(Vy)), V3, maxloc(abs(Vz)),               &
        T11, maxloc(abs(Txx)), T22, maxloc(abs(Tyy)), T33, maxloc(abs(Tzz)),          &
        T12, maxloc(abs(Txy)), T13, maxloc(abs(Txz)), T23, maxloc(abs(Tyz))
#ifdef VERBOSE
      write(fid_out,"(i5,9(es12.5,3i5))") ntime,                                      &
         V1, maxloc(abs(Vx)), V2, maxloc(abs(Vy)), V3, maxloc(abs(Vz)),               &
        T11, maxloc(abs(Txx)), T22, maxloc(abs(Tyy)), T33, maxloc(abs(Tzz)),          &
        T12, maxloc(abs(Txy)), T13, maxloc(abs(Txz)), T23, maxloc(abs(Tyz))
#endif
      call MPI_ABORT(SWMPI_COMM,1,ierr)
   end if
end if
end subroutine macdrp_check

subroutine coef_fdxy2fdz
  use math_mod, only : invert
  real(SP),dimension(SEIS_GEO,SEIS_GEO) :: A,B,C
  real(SP) :: lam,miu,lam2mu
  real(SP) :: e11,e12,e13,e21,e22,e23,e31,e32,e33
  integer i,j,k
  k=nk2
  do j=nj1,nj2
  do i=ni1,ni2
     lam=lambda(i,j,k); miu =mu(i,j,k); lam2mu=lam+2.0*miu
       e11=  xi_x(i,j,k)
       e12=  xi_y(i,j,k)
       e13=  xi_z(i,j,k)
       e21= eta_x(i,j,k)
       e22= eta_y(i,j,k)
       e23= eta_z(i,j,k)
       e31=zeta_x(i,j,k)
       e32=zeta_y(i,j,k)
       e33=zeta_z(i,j,k)
       !----
       A(1,1)=lam2mu*e31*e31+miu*(e32*e32+e33*e33)
       A(1,2)=lam*e31*e32+miu*e32*e31
       A(1,3)=lam*e31*e33+miu*e33*e31
       A(2,1)=lam*e32*e31+miu*e31*e32
       A(2,2)=lam2mu*e32*e32+miu*(e31*e31+e33*e33)
       A(2,3)=lam*e32*e33+miu*e33*e32
       A(3,1)=lam*e33*e31+miu*e31*e33
       A(3,2)=lam*e33*e32+miu*e32*e33
       A(3,3)=lam2mu*e33*e33+miu*(e31*e31+e32*e32)
       call invert(A)
       !----
       B(1,1)=lam2mu*e31*e11+miu*(e32*e12+e33*e13)
       B(1,2)=lam*e31*e12+miu*e32*e11
       B(1,3)=lam*e31*e13+miu*e33*e11
       B(2,1)=lam*e32*e11+miu*e31*e12
       B(2,2)=lam2mu*e32*e12+miu*(e31*e11+e33*e13)
       B(2,3)=lam*e32*e13+miu*e33*e12
       B(3,1)=lam*e33*e11+miu*e31*e13
       B(3,2)=lam*e33*e12+miu*e32*e13
       B(3,3)=lam2mu*e33*e13+miu*(e31*e11+e32*e12)
       B=-B
       !----
       C(1,1)=lam2mu*e31*e21+miu*(e32*e22+e33*e23)
       C(1,2)=lam*e31*e22+miu*e32*e21
       C(1,3)=lam*e31*e23+miu*e33*e21
       C(2,1)=lam*e32*e21+miu*e31*e22
       C(2,2)=lam2mu*e32*e22+miu*(e31*e21+e33*e23)
       C(2,3)=lam*e32*e23+miu*e33*e22
       C(3,1)=lam*e33*e21+miu*e31*e23
       C(3,2)=lam*e33*e22+miu*e32*e23
       C(3,3)=lam2mu*e33*e23+miu*(e31*e21+e32*e22)
       C=-C
       !----
       matVx2Vz(1:SEIS_GEO,1:SEIS_GEO,i,j)=matmul(A,B)
       matVy2Vz(1:SEIS_GEO,1:SEIS_GEO,i,j)=matmul(A,C)
       matF2Vz (1:SEIS_GEO,1:SEIS_GEO,i,j)=A
  end do
  end do
end subroutine coef_fdxy2fdz

!{----------- 4-6 LDDRK stages ----------------

subroutine macdrp_syn

 integer i,j,k
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
  mTxx(i,j,k)=Txx(i,j,k)
  mTyy(i,j,k)=Tyy(i,j,k)
  mTxy(i,j,k)=Txy(i,j,k)
  mVx (i,j,k)=Vx (i,j,k)
  mVy (i,j,k)=Vy (i,j,k)
  mTzz(i,j,k)=Tzz (i,j,k)
  mTxz(i,j,k)=Txz (i,j,k)
  mTyz(i,j,k)=Tyz (i,j,k)
  mVz (i,j,k)=Vz  (i,j,k)
 end do
 end do
 end do
!$OMP END PARALLEL DO
end subroutine macdrp_syn

!-- Generic RK --
subroutine macdrp_RK_beg(rka,rkb)
 real(SP),intent(in) :: rka,rkb
 real(SP) :: a,b
 integer i,j,k
 a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
    Txx(i,j,k)=mTxx(i,j,k)+a*hTxx(i,j,k)
    Tyy(i,j,k)=mTyy(i,j,k)+a*hTyy(i,j,k)
    Tzz(i,j,k)=mTzz(i,j,k)+a*hTzz(i,j,k)
    Txy(i,j,k)=mTxy(i,j,k)+a*hTxy(i,j,k)
    Txz(i,j,k)=mTxz(i,j,k)+a*hTxz(i,j,k)
    Tyz(i,j,k)=mTyz(i,j,k)+a*hTyz(i,j,k)
    Vx (i,j,k)=mVx (i,j,k)+a*hVx (i,j,k)
    Vy (i,j,k)=mVy (i,j,k)+a*hVy (i,j,k)
    Vz (i,j,k)=mVz (i,j,k)+a*hVz (i,j,k)

    tTxx(i,j,k)=mTxx(i,j,k)+b*hTxx(i,j,k)
    tTyy(i,j,k)=mTyy(i,j,k)+b*hTyy(i,j,k)
    tTzz(i,j,k)=mTzz(i,j,k)+b*hTzz(i,j,k)
    tTxy(i,j,k)=mTxy(i,j,k)+b*hTxy(i,j,k)
    tTxz(i,j,k)=mTxz(i,j,k)+b*hTxz(i,j,k)
    tTyz(i,j,k)=mTyz(i,j,k)+b*hTyz(i,j,k)
    tVx (i,j,k)=mVx (i,j,k)+b*hVx (i,j,k)
    tVy (i,j,k)=mVy (i,j,k)+b*hVy (i,j,k)
    tVz (i,j,k)=mVz (i,j,k)+b*hVz (i,j,k)
 end do
 end do
 end do
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
call free_charac
call free_extrap
#endif
end subroutine macdrp_RK_beg

subroutine macdrp_RK_inn(rka,rkb)
 real(SP),intent(in) :: rka,rkb
 real(SP) :: a,b
 integer i,j,k
 a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
    Txx(i,j,k)=mTxx(i,j,k)+a*hTxx(i,j,k)
    Tyy(i,j,k)=mTyy(i,j,k)+a*hTyy(i,j,k)
    Tzz(i,j,k)=mTzz(i,j,k)+a*hTzz(i,j,k)
    Txy(i,j,k)=mTxy(i,j,k)+a*hTxy(i,j,k)
    Txz(i,j,k)=mTxz(i,j,k)+a*hTxz(i,j,k)
    Tyz(i,j,k)=mTyz(i,j,k)+a*hTyz(i,j,k)
    Vx (i,j,k)=mVx (i,j,k)+a*hVx (i,j,k)
    Vy (i,j,k)=mVy (i,j,k)+a*hVy (i,j,k)
    Vz (i,j,k)=mVz (i,j,k)+a*hVz (i,j,k)

    tTxx(i,j,k)=tTxx(i,j,k)+b*hTxx(i,j,k)
    tTyy(i,j,k)=tTyy(i,j,k)+b*hTyy(i,j,k)
    tTzz(i,j,k)=tTzz(i,j,k)+b*hTzz(i,j,k)
    tTxy(i,j,k)=tTxy(i,j,k)+b*hTxy(i,j,k)
    tTxz(i,j,k)=tTxz(i,j,k)+b*hTxz(i,j,k)
    tTyz(i,j,k)=tTyz(i,j,k)+b*hTyz(i,j,k)
    tVx (i,j,k)=tVx (i,j,k)+b*hVx (i,j,k)
    tVy (i,j,k)=tVy (i,j,k)+b*hVy (i,j,k)
    tVz (i,j,k)=tVz (i,j,k)+b*hVz (i,j,k)
 end do
 end do
 end do
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
call free_charac
call free_extrap
#endif
end subroutine macdrp_RK_inn

subroutine macdrp_RK_fin(rkb)
 real(SP),intent(in) :: rkb
 real(SP) :: b
 integer i,j,k
 b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
    Txx(i,j,k)=tTxx(i,j,k)+b*hTxx(i,j,k)
    Tyy(i,j,k)=tTyy(i,j,k)+b*hTyy(i,j,k)
    Tzz(i,j,k)=tTzz(i,j,k)+b*hTzz(i,j,k)
    Txy(i,j,k)=tTxy(i,j,k)+b*hTxy(i,j,k)
    Txz(i,j,k)=tTxz(i,j,k)+b*hTxz(i,j,k)
    Tyz(i,j,k)=tTyz(i,j,k)+b*hTyz(i,j,k)
    Vx (i,j,k)=tVx (i,j,k)+b*hVx (i,j,k)
    Vy (i,j,k)=tVy (i,j,k)+b*hVy (i,j,k)
    Vz (i,j,k)=tVz (i,j,k)+b*hVz (i,j,k)
 end do
 end do
 end do
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
call free_charac
call free_extrap
#endif
end subroutine macdrp_RK_fin

subroutine atten_graves
 integer :: i,j,k
 real(SP) :: Qatt
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
 do k=nk1,nk2
 do j=nj1,nj2
 do i=ni1,ni2
 if (Qs(i,j,k)<QsINF) then
    Qatt=exp((-PI*QsF0*stept)/Qs(i,j,k))
    Txx(i,j,k)=Txx(i,j,k)*Qatt
    Tyy(i,j,k)=Tyy(i,j,k)*Qatt
    Tzz(i,j,k)=Tzz(i,j,k)*Qatt
    Txy(i,j,k)=Txy(i,j,k)*Qatt
    Txz(i,j,k)=Txz(i,j,k)*Qatt
    Tyz(i,j,k)=Tyz(i,j,k)*Qatt
    Vx (i,j,k)=Vx (i,j,k)*Qatt
    Vy (i,j,k)=Vy (i,j,k)*Qatt
    Vz (i,j,k)=Vz (i,j,k)*Qatt
 end if
 end do
 end do
 end do
!$OMP END PARALLEL DO
end subroutine atten_graves
!---------------------------------------------------}

!{------------------------ DRP/opt macdrp --------------------------------

!{----- wrapper of LxF_LyF_LzF ---------
subroutine in_LxF_LyF_LzF
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
  n=SEIS_GEO*2+1
  call MPI_STARTALL(36,reqXF,ierr)
  call MPI_STARTALL(36,reqYF,ierr)
  call MPI_STARTALL(36,reqZF,ierr)
  call LxF_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(36,reqXF,reqstat,ierr)
  call MPI_WAITALL(36,reqYF,reqstat,ierr)
  call MPI_WAITALL(36,reqZF,reqstat,ierr)

  do n=1,SEIS_GEO*2
  call LxF_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeTIMG
  if (freenode) call LxF_LyF_LzF_TIMG
#endif
#ifdef CondFreeVHOC
   if (freenode)  call LxF_LyF_LzF_VHOC
#endif
end subroutine in_LxF_LyF_LzF
subroutine in_LxB_LyB_LzB
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif

  !Txx=real(myid)
  n=SEIS_GEO*2+1
  !if (myid==0) then
  !  Vx=1.0;Vy=2.0;Vz=3.0;Txx=4.0;Tyy=5.0;Tzz=6.0;Txy=7.0;Txz=8.0;Tyz=9.0
  !end if

  call MPI_STARTALL(36,reqXB,ierr)
  call MPI_STARTALL(36,reqYB,ierr)
  call MPI_STARTALL(36,reqZB,ierr)
  call LxB_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(36,reqXB,reqstat,ierr)
  call MPI_WAITALL(36,reqYB,reqstat,ierr)
  call MPI_WAITALL(36,reqZB,reqstat,ierr)

  do n=1,SEIS_GEO*2
  call LxB_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeTIMG
  if (freenode) call LxB_LyB_LzB_TIMG
#endif
#ifdef CondFreeVHOC
   if (freenode)  call LxB_LyB_LzB_VHOC
#endif
end subroutine in_LxB_LyB_LzB
subroutine in_LxB_LyB_LzF
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif

  n=SEIS_GEO*2+1
  call MPI_STARTALL(36,reqXB,ierr)
  call MPI_STARTALL(36,reqYB,ierr)
  call MPI_STARTALL(36,reqZF,ierr)
  call LxB_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(36,reqXB,reqstat,ierr)
  call MPI_WAITALL(36,reqYB,reqstat,ierr)
  call MPI_WAITALL(36,reqZF,reqstat,ierr)

  do n=1,SEIS_GEO*2
  call LxB_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeTIMG
  if (freenode) call LxB_LyB_LzF_TIMG
#endif
#ifdef CondFreeVHOC
   if (freenode)  call LxB_LyB_LzF_VHOC
#endif
end subroutine in_LxB_LyB_LzF
subroutine in_LxF_LyF_LzB
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif

  n=SEIS_GEO*2+1
  call MPI_STARTALL(36,reqXF,ierr)
  call MPI_STARTALL(36,reqYF,ierr)
  call MPI_STARTALL(36,reqZB,ierr)
  call LxF_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(36,reqXF,reqstat,ierr)
  call MPI_WAITALL(36,reqYF,reqstat,ierr)
  call MPI_WAITALL(36,reqZB,reqstat,ierr)

  do n=1,SEIS_GEO*2
  call LxF_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeTIMG
  if (freenode) call LxF_LyF_LzB_TIMG
#endif
#ifdef CondFreeVHOC
   if (freenode)  call LxF_LyF_LzB_VHOC
#endif
end subroutine in_LxF_LyF_LzB
subroutine in_LxB_LyF_LzF
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
  n=SEIS_GEO*2+1
  call MPI_STARTALL(36,reqXB,ierr)
  call MPI_STARTALL(36,reqYF,ierr)
  call MPI_STARTALL(36,reqZF,ierr)
  call LxB_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(36,reqXB,reqstat,ierr)
  call MPI_WAITALL(36,reqYF,reqstat,ierr)
  call MPI_WAITALL(36,reqZF,reqstat,ierr)

  do n=1,SEIS_GEO*2
  call LxB_LyF_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeTIMG
  if (freenode) call LxB_LyF_LzF_TIMG
#endif
#ifdef CondFreeVHOC
   if (freenode)  call LxB_LyF_LzF_VHOC
#endif
end subroutine in_LxB_LyF_LzF
subroutine in_LxF_LyB_LzB
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
  n=SEIS_GEO*2+1
  call MPI_STARTALL(36,reqXF,ierr)
  call MPI_STARTALL(36,reqYB,ierr)
  call MPI_STARTALL(36,reqZB,ierr)
  call LxF_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(36,reqXF,reqstat,ierr)
  call MPI_WAITALL(36,reqYB,reqstat,ierr)
  call MPI_WAITALL(36,reqZB,reqstat,ierr)

  do n=1,SEIS_GEO*2
  call LxF_LyB_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeTIMG
  if (freenode) call LxF_LyB_LzB_TIMG
#endif
#ifdef CondFreeVHOC
   if (freenode)  call LxF_LyB_LzB_VHOC
#endif
end subroutine in_LxF_LyB_LzB
subroutine in_LxF_LyB_LzF
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
  n=SEIS_GEO*2+1
  call MPI_STARTALL(36,reqXF,ierr)
  call MPI_STARTALL(36,reqYB,ierr)
  call MPI_STARTALL(36,reqZF,ierr)
  call LxF_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(36,reqXF,reqstat,ierr)
  call MPI_WAITALL(36,reqYB,reqstat,ierr)
  call MPI_WAITALL(36,reqZF,reqstat,ierr)

  do n=1,SEIS_GEO*2
  call LxF_LyB_LzF( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeTIMG
  if (freenode) call LxF_LyB_LzF_TIMG
#endif
#ifdef CondFreeVHOC
   if (freenode)  call LxF_LyB_LzF_VHOC
#endif
end subroutine in_LxF_LyB_LzF
subroutine in_LxB_LyF_LzB
  integer n
#ifdef MPIBARRIER
  integer ierr
  call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
  !Txx=real(myid)
  n=SEIS_GEO*2+1
  call MPI_STARTALL(36,reqXB,ierr)
  call MPI_STARTALL(36,reqYF,ierr)
  call MPI_STARTALL(36,reqZB,ierr)
  call LxB_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  call MPI_WAITALL(36,reqXB,reqstat,ierr)
  call MPI_WAITALL(36,reqYF,reqstat,ierr)
  call MPI_WAITALL(36,reqZB,reqstat,ierr)

  do n=1,SEIS_GEO*2
  call LxB_LyF_LzB( indx(1,n),indx(2,n), &
                    indx(3,n),indx(4,n), &
                    indx(5,n),indx(6,n) )
  end do

#ifdef CondFreeTIMG
  if (freenode) call LxB_LyF_LzB_TIMG
#endif
#ifdef CondFreeVHOC
   if (freenode)  call LxB_LyF_LzB_VHOC
#endif
end subroutine in_LxB_LyF_LzB
!-- private subroutine --
subroutine LxF_LyF_LzF(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz, &
!$OMP   DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )/steph
   DxTyy = (              &
     m3d_FDxF1(Tyy,i,j,k) &
     m3d_FDxF2(Tyy,i,j,k) &
     )/steph
   DxTzz = (              &
     m3d_FDxF1(Tzz,i,j,k) &
     m3d_FDxF2(Tzz,i,j,k) &
     )/steph
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )/steph
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )/steph
   DxTyz = (              &
     m3d_FDxF1(Tyz,i,j,k) &
     m3d_FDxF2(Tyz,i,j,k) &
     )/steph
   DxVx =  (              &
     m3d_FDxF1(Vx,i,j,k)  &
     m3d_FDxF2(Vx,i,j,k)  &
     )/steph
   DxVy = (               &
     m3d_FDxF1(Vy,i,j,k)  &
     m3d_FDxF2(Vy,i,j,k)  &
     )/steph
   DxVz = (               &
     m3d_FDxF1(Vz,i,j,k)  &
     m3d_FDxF2(Vz,i,j,k)  &
     )/steph

   DyTxx = (              &
     m3d_FDyF1(Txx,i,j,k) &
     m3d_FDyF2(Txx,i,j,k) &
     )/steph
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )/steph
   DyTzz = (              &
     m3d_FDyF1(Tzz,i,j,k) &
     m3d_FDyF2(Tzz,i,j,k) &
     )/steph
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )/steph
   DyTxz = (              &
     m3d_FDyF1(Txz,i,j,k) &
     m3d_FDyF2(Txz,i,j,k) &
     )/steph
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )/steph
   DyVx = (               &
     m3d_FDyF1(Vx,i,j,k)  &
     m3d_FDyF2(Vx,i,j,k)  &
     )/steph
   DyVy = (               &
     m3d_FDyF1(Vy,i,j,k)  &
     m3d_FDyF2(Vy,i,j,k)  &
     )/steph
   DyVz = (               &
     m3d_FDyF1(Vz,i,j,k)  &
     m3d_FDyF2(Vz,i,j,k)  &
     )/steph

   DzTxx = (              &
     m3d_FDzF1(Txx,i,j,k) &
     m3d_FDzF2(Txx,i,j,k) &
     )/steph
   DzTyy = (              &
     m3d_FDzF1(Tyy,i,j,k) &
     m3d_FDzF2(Tyy,i,j,k) &
     )/steph
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )/steph
   DzTxy = (              &
     m3d_FDzF1(Txy,i,j,k) &
     m3d_FDzF2(Txy,i,j,k) &
     )/steph
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )/steph
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )/steph

#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzF1(Vx,i,j,k)  &
     m22_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m22_FDzF1(Vy,i,j,k)  &
     m22_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m22_FDzF1(Vz,i,j,k)  &
     m22_FDzF2(Vz,i,j,k)  &
     )/steph
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzF1(Vx,i,j,k)  &
     m24_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m24_FDzF1(Vy,i,j,k)  &
     m24_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m24_FDzF1(Vz,i,j,k)  &
     m24_FDzF2(Vz,i,j,k)  &
     )/steph
   else
#endif
   DzVx = (               &
     m3d_FDzF1(Vx,i,j,k)  &
     m3d_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m3d_FDzF1(Vy,i,j,k)  &
     m3d_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m3d_FDzF1(Vz,i,j,k)  &
     m3d_FDzF2(Vz,i,j,k)  &
     )/steph
#ifdef CondFreeVLOW
   end if
#endif

#ifndef CondFreeVHOC
#ifndef CondFreeCharac
if (freenode .and. k==nk2) then
   DzVx = matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz &
         +matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy = matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz &
         +matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz = matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz &
         +matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx=DzVx                          &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy=DzVy                          &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz=DzVz                          &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)
end if
#endif
#endif

   hVx(i,j,k)= rrho*(                          &
           xix*DxTxx + xiy*DxTxy + xiz*DxTxz   &
          +etx*DyTxx + ety*DyTxy + etz*DyTxz   &
          +ztx*DzTxx + zty*DzTxy + ztz*DzTxz )
   hVy(i,j,k)= rrho*(                          &
           xix*DxTxy + xiy*DxTyy + xiz*DxTyz   &
          +etx*DyTxy + ety*DyTyy + etz*DyTyz   &
          +ztx*DzTxy + zty*DzTyy + ztz*DzTyz )
   hVz(i,j,k)= rrho*(                          &
           xix*DxTxz + xiy*DxTyz + xiz*DxTzz   &
          +etx*DyTxz + ety*DyTyz + etz*DyTzz   &
          +ztx*DzTxz + zty*DzTyz + ztz*DzTzz )

   hTxx(i,j,k)= lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz &
               +lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz &
               +lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz
   hTyy(i,j,k)= lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz &
               +lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz &
               +lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz
   hTzz(i,j,k)= lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz &
               +lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz &
               +lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz
   hTxy(i,j,k)= miu *(              &
                xiy*DxVx + xix*DxVy &
               +ety*DyVx + etx*DyVy &
               +zty*DzVx + ztx*DzVy &
               )
   hTxz(i,j,k)= miu *(              &
                xiz*DxVx + xix*DxVz &
               +etz*DyVx + etx*DyVz &
               +ztz*DzVx + ztx*DzVz &
               )
   hTyz(i,j,k)= miu *(              &
                xiz*DxVy + xiy*DxVz &
               +etz*DyVy + ety*DyVz &
               +ztz*DzVy + zty*DzVz &
               )
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxF_LyF_LzF

subroutine LxB_LyB_LzB(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz, &
!$OMP   DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )/steph
   DxTyy = (              &
     m3d_FDxB1(Tyy,i,j,k) &
     m3d_FDxB2(Tyy,i,j,k) &
     )/steph
   DxTzz = (              &
     m3d_FDxB1(Tzz,i,j,k) &
     m3d_FDxB2(Tzz,i,j,k) &
     )/steph
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )/steph
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )/steph
   DxTyz = (              &
     m3d_FDxB1(Tyz,i,j,k) &
     m3d_FDxB2(Tyz,i,j,k) &
     )/steph
   DxVx =  (              &
     m3d_FDxB1(Vx,i,j,k)  &
     m3d_FDxB2(Vx,i,j,k)  &
     )/steph
   DxVy = (               &
     m3d_FDxB1(Vy,i,j,k)  &
     m3d_FDxB2(Vy,i,j,k)  &
     )/steph
   DxVz = (               &
     m3d_FDxB1(Vz,i,j,k)  &
     m3d_FDxB2(Vz,i,j,k)  &
     )/steph

   DyTxx = (              &
     m3d_FDyB1(Txx,i,j,k) &
     m3d_FDyB2(Txx,i,j,k) &
     )/steph
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )/steph
   DyTzz = (              &
     m3d_FDyB1(Tzz,i,j,k) &
     m3d_FDyB2(Tzz,i,j,k) &
     )/steph
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )/steph
   DyTxz = (              &
     m3d_FDyB1(Txz,i,j,k) &
     m3d_FDyB2(Txz,i,j,k) &
     )/steph
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )/steph
   DyVx = (               &
     m3d_FDyB1(Vx,i,j,k)  &
     m3d_FDyB2(Vx,i,j,k)  &
     )/steph
   DyVy = (               &
     m3d_FDyB1(Vy,i,j,k)  &
     m3d_FDyB2(Vy,i,j,k)  &
     )/steph
   DyVz = (               &
     m3d_FDyB1(Vz,i,j,k)  &
     m3d_FDyB2(Vz,i,j,k)  &
     )/steph

   DzTxx = (              &
     m3d_FDzB1(Txx,i,j,k) &
     m3d_FDzB2(Txx,i,j,k) &
     )/steph
   DzTyy = (              &
     m3d_FDzB1(Tyy,i,j,k) &
     m3d_FDzB2(Tyy,i,j,k) &
     )/steph
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )/steph
   DzTxy = (              &
     m3d_FDzB1(Txy,i,j,k) &
     m3d_FDzB2(Txy,i,j,k) &
     )/steph
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )/steph
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )/steph

#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzB1(Vx,i,j,k)  &
     m22_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m22_FDzB1(Vy,i,j,k)  &
     m22_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m22_FDzB1(Vz,i,j,k)  &
     m22_FDzB2(Vz,i,j,k)  &
     )/steph
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzB1(Vx,i,j,k)  &
     m24_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m24_FDzB1(Vy,i,j,k)  &
     m24_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m24_FDzB1(Vz,i,j,k)  &
     m24_FDzB2(Vz,i,j,k)  &
     )/steph
   else
#endif
   DzVx = (               &
     m3d_FDzB1(Vx,i,j,k)  &
     m3d_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m3d_FDzB1(Vy,i,j,k)  &
     m3d_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m3d_FDzB1(Vz,i,j,k)  &
     m3d_FDzB2(Vz,i,j,k)  &
     )/steph
#ifdef CondFreeVLOW
   end if
#endif

#ifndef CondFreeVHOC
#ifndef CondFreeCharac
if (freenode .and. k==nk2) then
   DzVx = matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz &
         +matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy = matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz &
         +matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz = matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz &
         +matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx=DzVx                          &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy=DzVy                          &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz=DzVz                          &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)
end if
#endif
#endif

   hVx(i,j,k)= rrho*(                          &
           xix*DxTxx + xiy*DxTxy + xiz*DxTxz   &
          +etx*DyTxx + ety*DyTxy + etz*DyTxz   &
          +ztx*DzTxx + zty*DzTxy + ztz*DzTxz )
   hVy(i,j,k)= rrho*(                          &
           xix*DxTxy + xiy*DxTyy + xiz*DxTyz   &
          +etx*DyTxy + ety*DyTyy + etz*DyTyz   &
          +ztx*DzTxy + zty*DzTyy + ztz*DzTyz )
   hVz(i,j,k)= rrho*(                          &
           xix*DxTxz + xiy*DxTyz + xiz*DxTzz   &
          +etx*DyTxz + ety*DyTyz + etz*DyTzz   &
          +ztx*DzTxz + zty*DzTyz + ztz*DzTzz )

   hTxx(i,j,k)= lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz &
               +lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz &
               +lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz
   hTyy(i,j,k)= lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz &
               +lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz &
               +lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz
   hTzz(i,j,k)= lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz &
               +lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz &
               +lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz
   hTxy(i,j,k)= miu *(              &
                xiy*DxVx + xix*DxVy &
               +ety*DyVx + etx*DyVy &
               +zty*DzVx + ztx*DzVy &
               )
   hTxz(i,j,k)= miu *(              &
                xiz*DxVx + xix*DxVz &
               +etz*DyVx + etx*DyVz &
               +ztz*DzVx + ztx*DzVz &
               )
   hTyz(i,j,k)= miu *(              &
                xiz*DxVy + xiy*DxVz &
               +etz*DyVy + ety*DyVz &
               +ztz*DzVy + zty*DzVz &
               )
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxB_LyB_LzB

subroutine LxF_LyF_LzB(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz, &
!$OMP   DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )/steph
   DxTyy = (              &
     m3d_FDxF1(Tyy,i,j,k) &
     m3d_FDxF2(Tyy,i,j,k) &
     )/steph
   DxTzz = (              &
     m3d_FDxF1(Tzz,i,j,k) &
     m3d_FDxF2(Tzz,i,j,k) &
     )/steph
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )/steph
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )/steph
   DxTyz = (              &
     m3d_FDxF1(Tyz,i,j,k) &
     m3d_FDxF2(Tyz,i,j,k) &
     )/steph
   DxVx =  (              &
     m3d_FDxF1(Vx,i,j,k)  &
     m3d_FDxF2(Vx,i,j,k)  &
     )/steph
   DxVy = (               &
     m3d_FDxF1(Vy,i,j,k)  &
     m3d_FDxF2(Vy,i,j,k)  &
     )/steph
   DxVz = (               &
     m3d_FDxF1(Vz,i,j,k)  &
     m3d_FDxF2(Vz,i,j,k)  &
     )/steph

   DyTxx = (              &
     m3d_FDyF1(Txx,i,j,k) &
     m3d_FDyF2(Txx,i,j,k) &
     )/steph
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )/steph
   DyTzz = (              &
     m3d_FDyF1(Tzz,i,j,k) &
     m3d_FDyF2(Tzz,i,j,k) &
     )/steph
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )/steph
   DyTxz = (              &
     m3d_FDyF1(Txz,i,j,k) &
     m3d_FDyF2(Txz,i,j,k) &
     )/steph
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )/steph
   DyVx = (               &
     m3d_FDyF1(Vx,i,j,k)  &
     m3d_FDyF2(Vx,i,j,k)  &
     )/steph
   DyVy = (               &
     m3d_FDyF1(Vy,i,j,k)  &
     m3d_FDyF2(Vy,i,j,k)  &
     )/steph
   DyVz = (               &
     m3d_FDyF1(Vz,i,j,k)  &
     m3d_FDyF2(Vz,i,j,k)  &
     )/steph

   DzTxx = (              &
     m3d_FDzB1(Txx,i,j,k) &
     m3d_FDzB2(Txx,i,j,k) &
     )/steph
   DzTyy = (              &
     m3d_FDzB1(Tyy,i,j,k) &
     m3d_FDzB2(Tyy,i,j,k) &
     )/steph
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )/steph
   DzTxy = (              &
     m3d_FDzB1(Txy,i,j,k) &
     m3d_FDzB2(Txy,i,j,k) &
     )/steph
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )/steph
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )/steph

#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzB1(Vx,i,j,k)  &
     m22_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m22_FDzB1(Vy,i,j,k)  &
     m22_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m22_FDzB1(Vz,i,j,k)  &
     m22_FDzB2(Vz,i,j,k)  &
     )/steph
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzB1(Vx,i,j,k)  &
     m24_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m24_FDzB1(Vy,i,j,k)  &
     m24_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m24_FDzB1(Vz,i,j,k)  &
     m24_FDzB2(Vz,i,j,k)  &
     )/steph
   else
#endif
   DzVx = (               &
     m3d_FDzB1(Vx,i,j,k)  &
     m3d_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m3d_FDzB1(Vy,i,j,k)  &
     m3d_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m3d_FDzB1(Vz,i,j,k)  &
     m3d_FDzB2(Vz,i,j,k)  &
     )/steph
#ifdef CondFreeVLOW
   end if
#endif

#ifndef CondFreeVHOC
#ifndef CondFreeCharac
if (freenode .and. k==nk2) then
   DzVx = matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz &
         +matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy = matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz &
         +matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz = matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz &
         +matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx=DzVx                          &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy=DzVy                          &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz=DzVz                          &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)
end if
#endif
#endif

   hVx(i,j,k)= rrho*(                          &
           xix*DxTxx + xiy*DxTxy + xiz*DxTxz   &
          +etx*DyTxx + ety*DyTxy + etz*DyTxz   &
          +ztx*DzTxx + zty*DzTxy + ztz*DzTxz )
   hVy(i,j,k)= rrho*(                          &
           xix*DxTxy + xiy*DxTyy + xiz*DxTyz   &
          +etx*DyTxy + ety*DyTyy + etz*DyTyz   &
          +ztx*DzTxy + zty*DzTyy + ztz*DzTyz )
   hVz(i,j,k)= rrho*(                          &
           xix*DxTxz + xiy*DxTyz + xiz*DxTzz   &
          +etx*DyTxz + ety*DyTyz + etz*DyTzz   &
          +ztx*DzTxz + zty*DzTyz + ztz*DzTzz )

   hTxx(i,j,k)= lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz &
               +lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz &
               +lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz
   hTyy(i,j,k)= lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz &
               +lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz &
               +lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz
   hTzz(i,j,k)= lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz &
               +lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz &
               +lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz
   hTxy(i,j,k)= miu *(              &
                xiy*DxVx + xix*DxVy &
               +ety*DyVx + etx*DyVy &
               +zty*DzVx + ztx*DzVy &
               )
   hTxz(i,j,k)= miu *(              &
                xiz*DxVx + xix*DxVz &
               +etz*DyVx + etx*DyVz &
               +ztz*DzVx + ztx*DzVz &
               )
   hTyz(i,j,k)= miu *(              &
                xiz*DxVy + xiy*DxVz &
               +etz*DyVy + ety*DyVz &
               +ztz*DzVy + zty*DzVz &
               )
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxF_LyF_LzB

subroutine LxB_LyB_LzF(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz, &
!$OMP   DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )/steph
   DxTyy = (              &
     m3d_FDxB1(Tyy,i,j,k) &
     m3d_FDxB2(Tyy,i,j,k) &
     )/steph
   DxTzz = (              &
     m3d_FDxB1(Tzz,i,j,k) &
     m3d_FDxB2(Tzz,i,j,k) &
     )/steph
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )/steph
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )/steph
   DxTyz = (              &
     m3d_FDxB1(Tyz,i,j,k) &
     m3d_FDxB2(Tyz,i,j,k) &
     )/steph
   DxVx =  (              &
     m3d_FDxB1(Vx,i,j,k)  &
     m3d_FDxB2(Vx,i,j,k)  &
     )/steph
   DxVy = (               &
     m3d_FDxB1(Vy,i,j,k)  &
     m3d_FDxB2(Vy,i,j,k)  &
     )/steph
   DxVz = (               &
     m3d_FDxB1(Vz,i,j,k)  &
     m3d_FDxB2(Vz,i,j,k)  &
     )/steph

   DyTxx = (              &
     m3d_FDyB1(Txx,i,j,k) &
     m3d_FDyB2(Txx,i,j,k) &
     )/steph
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )/steph
   DyTzz = (              &
     m3d_FDyB1(Tzz,i,j,k) &
     m3d_FDyB2(Tzz,i,j,k) &
     )/steph
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )/steph
   DyTxz = (              &
     m3d_FDyB1(Txz,i,j,k) &
     m3d_FDyB2(Txz,i,j,k) &
     )/steph
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )/steph
   DyVx = (               &
     m3d_FDyB1(Vx,i,j,k)  &
     m3d_FDyB2(Vx,i,j,k)  &
     )/steph
   DyVy = (               &
     m3d_FDyB1(Vy,i,j,k)  &
     m3d_FDyB2(Vy,i,j,k)  &
     )/steph
   DyVz = (               &
     m3d_FDyB1(Vz,i,j,k)  &
     m3d_FDyB2(Vz,i,j,k)  &
     )/steph

   DzTxx = (              &
     m3d_FDzF1(Txx,i,j,k) &
     m3d_FDzF2(Txx,i,j,k) &
     )/steph
   DzTyy = (              &
     m3d_FDzF1(Tyy,i,j,k) &
     m3d_FDzF2(Tyy,i,j,k) &
     )/steph
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )/steph
   DzTxy = (              &
     m3d_FDzF1(Txy,i,j,k) &
     m3d_FDzF2(Txy,i,j,k) &
     )/steph
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )/steph
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )/steph

#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzF1(Vx,i,j,k)  &
     m22_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m22_FDzF1(Vy,i,j,k)  &
     m22_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m22_FDzF1(Vz,i,j,k)  &
     m22_FDzF2(Vz,i,j,k)  &
     )/steph
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzF1(Vx,i,j,k)  &
     m24_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m24_FDzF1(Vy,i,j,k)  &
     m24_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m24_FDzF1(Vz,i,j,k)  &
     m24_FDzF2(Vz,i,j,k)  &
     )/steph
   else
#endif
   DzVx = (               &
     m3d_FDzF1(Vx,i,j,k)  &
     m3d_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m3d_FDzF1(Vy,i,j,k)  &
     m3d_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m3d_FDzF1(Vz,i,j,k)  &
     m3d_FDzF2(Vz,i,j,k)  &
     )/steph
#ifdef CondFreeVLOW
   end if
#endif

#ifndef CondFreeVHOC
#ifndef CondFreeCharac
if (freenode .and. k==nk2) then
   DzVx = matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz &
         +matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy = matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz &
         +matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz = matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz &
         +matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx=DzVx                          &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy=DzVy                          &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz=DzVz                          &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)
end if
#endif
#endif

   hVx(i,j,k)= rrho*(                          &
           xix*DxTxx + xiy*DxTxy + xiz*DxTxz   &
          +etx*DyTxx + ety*DyTxy + etz*DyTxz   &
          +ztx*DzTxx + zty*DzTxy + ztz*DzTxz )
   hVy(i,j,k)= rrho*(                          &
           xix*DxTxy + xiy*DxTyy + xiz*DxTyz   &
          +etx*DyTxy + ety*DyTyy + etz*DyTyz   &
          +ztx*DzTxy + zty*DzTyy + ztz*DzTyz )
   hVz(i,j,k)= rrho*(                          &
           xix*DxTxz + xiy*DxTyz + xiz*DxTzz   &
          +etx*DyTxz + ety*DyTyz + etz*DyTzz   &
          +ztx*DzTxz + zty*DzTyz + ztz*DzTzz )

   hTxx(i,j,k)= lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz &
               +lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz &
               +lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz
   hTyy(i,j,k)= lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz &
               +lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz &
               +lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz
   hTzz(i,j,k)= lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz &
               +lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz &
               +lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz
   hTxy(i,j,k)= miu *(              &
                xiy*DxVx + xix*DxVy &
               +ety*DyVx + etx*DyVy &
               +zty*DzVx + ztx*DzVy &
               )
   hTxz(i,j,k)= miu *(              &
                xiz*DxVx + xix*DxVz &
               +etz*DyVx + etx*DyVz &
               +ztz*DzVx + ztx*DzVz &
               )
   hTyz(i,j,k)= miu *(              &
                xiz*DxVy + xiy*DxVz &
               +etz*DyVy + ety*DyVz &
               +ztz*DzVy + zty*DzVz &
               )
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxB_LyB_LzF

subroutine LxB_LyF_LzF(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz, &
!$OMP   DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )/steph
   DxTyy = (              &
     m3d_FDxB1(Tyy,i,j,k) &
     m3d_FDxB2(Tyy,i,j,k) &
     )/steph
   DxTzz = (              &
     m3d_FDxB1(Tzz,i,j,k) &
     m3d_FDxB2(Tzz,i,j,k) &
     )/steph
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )/steph
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )/steph
   DxTyz = (              &
     m3d_FDxB1(Tyz,i,j,k) &
     m3d_FDxB2(Tyz,i,j,k) &
     )/steph
   DxVx =  (              &
     m3d_FDxB1(Vx,i,j,k)  &
     m3d_FDxB2(Vx,i,j,k)  &
     )/steph
   DxVy = (               &
     m3d_FDxB1(Vy,i,j,k)  &
     m3d_FDxB2(Vy,i,j,k)  &
     )/steph
   DxVz = (               &
     m3d_FDxB1(Vz,i,j,k)  &
     m3d_FDxB2(Vz,i,j,k)  &
     )/steph

   DyTxx = (              &
     m3d_FDyF1(Txx,i,j,k) &
     m3d_FDyF2(Txx,i,j,k) &
     )/steph
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )/steph
   DyTzz = (              &
     m3d_FDyF1(Tzz,i,j,k) &
     m3d_FDyF2(Tzz,i,j,k) &
     )/steph
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )/steph
   DyTxz = (              &
     m3d_FDyF1(Txz,i,j,k) &
     m3d_FDyF2(Txz,i,j,k) &
     )/steph
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )/steph
   DyVx = (               &
     m3d_FDyF1(Vx,i,j,k)  &
     m3d_FDyF2(Vx,i,j,k)  &
     )/steph
   DyVy = (               &
     m3d_FDyF1(Vy,i,j,k)  &
     m3d_FDyF2(Vy,i,j,k)  &
     )/steph
   DyVz = (               &
     m3d_FDyF1(Vz,i,j,k)  &
     m3d_FDyF2(Vz,i,j,k)  &
     )/steph

   DzTxx = (              &
     m3d_FDzF1(Txx,i,j,k) &
     m3d_FDzF2(Txx,i,j,k) &
     )/steph
   DzTyy = (              &
     m3d_FDzF1(Tyy,i,j,k) &
     m3d_FDzF2(Tyy,i,j,k) &
     )/steph
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )/steph
   DzTxy = (              &
     m3d_FDzF1(Txy,i,j,k) &
     m3d_FDzF2(Txy,i,j,k) &
     )/steph
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )/steph
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )/steph

#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzF1(Vx,i,j,k)  &
     m22_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m22_FDzF1(Vy,i,j,k)  &
     m22_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m22_FDzF1(Vz,i,j,k)  &
     m22_FDzF2(Vz,i,j,k)  &
     )/steph
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzF1(Vx,i,j,k)  &
     m24_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m24_FDzF1(Vy,i,j,k)  &
     m24_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m24_FDzF1(Vz,i,j,k)  &
     m24_FDzF2(Vz,i,j,k)  &
     )/steph
   else
#endif
   DzVx = (               &
     m3d_FDzF1(Vx,i,j,k)  &
     m3d_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m3d_FDzF1(Vy,i,j,k)  &
     m3d_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m3d_FDzF1(Vz,i,j,k)  &
     m3d_FDzF2(Vz,i,j,k)  &
     )/steph
#ifdef CondFreeVLOW
   end if
#endif

#ifndef CondFreeVHOC
#ifndef CondFreeCharac
if (freenode .and. k==nk2) then
   DzVx = matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz &
         +matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy = matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz &
         +matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz = matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz &
         +matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx=DzVx                          &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy=DzVy                          &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz=DzVz                          &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)
end if
#endif
#endif

   hVx(i,j,k)= rrho*(                          &
           xix*DxTxx + xiy*DxTxy + xiz*DxTxz   &
          +etx*DyTxx + ety*DyTxy + etz*DyTxz   &
          +ztx*DzTxx + zty*DzTxy + ztz*DzTxz )
   hVy(i,j,k)= rrho*(                          &
           xix*DxTxy + xiy*DxTyy + xiz*DxTyz   &
          +etx*DyTxy + ety*DyTyy + etz*DyTyz   &
          +ztx*DzTxy + zty*DzTyy + ztz*DzTyz )
   hVz(i,j,k)= rrho*(                          &
           xix*DxTxz + xiy*DxTyz + xiz*DxTzz   &
          +etx*DyTxz + ety*DyTyz + etz*DyTzz   &
          +ztx*DzTxz + zty*DzTyz + ztz*DzTzz )

   hTxx(i,j,k)= lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz &
               +lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz &
               +lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz
   hTyy(i,j,k)= lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz &
               +lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz &
               +lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz
   hTzz(i,j,k)= lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz &
               +lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz &
               +lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz
   hTxy(i,j,k)= miu *(              &
                xiy*DxVx + xix*DxVy &
               +ety*DyVx + etx*DyVy &
               +zty*DzVx + ztx*DzVy &
               )
   hTxz(i,j,k)= miu *(              &
                xiz*DxVx + xix*DxVz &
               +etz*DyVx + etx*DyVz &
               +ztz*DzVx + ztx*DzVz &
               )
   hTyz(i,j,k)= miu *(              &
                xiz*DxVy + xiy*DxVz &
               +etz*DyVy + ety*DyVz &
               +ztz*DzVy + zty*DzVz &
               )
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxB_LyF_LzF

subroutine LxF_LyB_LzB(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz, &
!$OMP   DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )/steph
   DxTyy = (              &
     m3d_FDxF1(Tyy,i,j,k) &
     m3d_FDxF2(Tyy,i,j,k) &
     )/steph
   DxTzz = (              &
     m3d_FDxF1(Tzz,i,j,k) &
     m3d_FDxF2(Tzz,i,j,k) &
     )/steph
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )/steph
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )/steph
   DxTyz = (              &
     m3d_FDxF1(Tyz,i,j,k) &
     m3d_FDxF2(Tyz,i,j,k) &
     )/steph
   DxVx =  (              &
     m3d_FDxF1(Vx,i,j,k)  &
     m3d_FDxF2(Vx,i,j,k)  &
     )/steph
   DxVy = (               &
     m3d_FDxF1(Vy,i,j,k)  &
     m3d_FDxF2(Vy,i,j,k)  &
     )/steph
   DxVz = (               &
     m3d_FDxF1(Vz,i,j,k)  &
     m3d_FDxF2(Vz,i,j,k)  &
     )/steph

   DyTxx = (              &
     m3d_FDyB1(Txx,i,j,k) &
     m3d_FDyB2(Txx,i,j,k) &
     )/steph
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )/steph
   DyTzz = (              &
     m3d_FDyB1(Tzz,i,j,k) &
     m3d_FDyB2(Tzz,i,j,k) &
     )/steph
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )/steph
   DyTxz = (              &
     m3d_FDyB1(Txz,i,j,k) &
     m3d_FDyB2(Txz,i,j,k) &
     )/steph
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )/steph
   DyVx = (               &
     m3d_FDyB1(Vx,i,j,k)  &
     m3d_FDyB2(Vx,i,j,k)  &
     )/steph
   DyVy = (               &
     m3d_FDyB1(Vy,i,j,k)  &
     m3d_FDyB2(Vy,i,j,k)  &
     )/steph
   DyVz = (               &
     m3d_FDyB1(Vz,i,j,k)  &
     m3d_FDyB2(Vz,i,j,k)  &
     )/steph

   DzTxx = (              &
     m3d_FDzB1(Txx,i,j,k) &
     m3d_FDzB2(Txx,i,j,k) &
     )/steph
   DzTyy = (              &
     m3d_FDzB1(Tyy,i,j,k) &
     m3d_FDzB2(Tyy,i,j,k) &
     )/steph
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )/steph
   DzTxy = (              &
     m3d_FDzB1(Txy,i,j,k) &
     m3d_FDzB2(Txy,i,j,k) &
     )/steph
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )/steph
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )/steph

#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzB1(Vx,i,j,k)  &
     m22_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m22_FDzB1(Vy,i,j,k)  &
     m22_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m22_FDzB1(Vz,i,j,k)  &
     m22_FDzB2(Vz,i,j,k)  &
     )/steph
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzB1(Vx,i,j,k)  &
     m24_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m24_FDzB1(Vy,i,j,k)  &
     m24_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m24_FDzB1(Vz,i,j,k)  &
     m24_FDzB2(Vz,i,j,k)  &
     )/steph
   else
#endif
   DzVx = (               &
     m3d_FDzB1(Vx,i,j,k)  &
     m3d_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m3d_FDzB1(Vy,i,j,k)  &
     m3d_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m3d_FDzB1(Vz,i,j,k)  &
     m3d_FDzB2(Vz,i,j,k)  &
     )/steph
#ifdef CondFreeVLOW
   end if
#endif

#ifndef CondFreeVHOC
#ifndef CondFreeCharac
if (freenode .and. k==nk2) then
   DzVx = matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz &
         +matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy = matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz &
         +matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz = matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz &
         +matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx=DzVx                          &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy=DzVy                          &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz=DzVz                          &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)
end if
#endif
#endif

   hVx(i,j,k)= rrho*(                          &
           xix*DxTxx + xiy*DxTxy + xiz*DxTxz   &
          +etx*DyTxx + ety*DyTxy + etz*DyTxz   &
          +ztx*DzTxx + zty*DzTxy + ztz*DzTxz )
   hVy(i,j,k)= rrho*(                          &
           xix*DxTxy + xiy*DxTyy + xiz*DxTyz   &
          +etx*DyTxy + ety*DyTyy + etz*DyTyz   &
          +ztx*DzTxy + zty*DzTyy + ztz*DzTyz )
   hVz(i,j,k)= rrho*(                          &
           xix*DxTxz + xiy*DxTyz + xiz*DxTzz   &
          +etx*DyTxz + ety*DyTyz + etz*DyTzz   &
          +ztx*DzTxz + zty*DzTyz + ztz*DzTzz )

   hTxx(i,j,k)= lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz &
               +lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz &
               +lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz
   hTyy(i,j,k)= lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz &
               +lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz &
               +lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz
   hTzz(i,j,k)= lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz &
               +lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz &
               +lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz
   hTxy(i,j,k)= miu *(              &
                xiy*DxVx + xix*DxVy &
               +ety*DyVx + etx*DyVy &
               +zty*DzVx + ztx*DzVy &
               )
   hTxz(i,j,k)= miu *(              &
                xiz*DxVx + xix*DxVz &
               +etz*DyVx + etx*DyVz &
               +ztz*DzVx + ztx*DzVz &
               )
   hTyz(i,j,k)= miu *(              &
                xiz*DxVy + xiy*DxVz &
               +etz*DyVy + ety*DyVz &
               +ztz*DzVy + zty*DzVz &
               )
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxF_LyB_LzB

subroutine LxF_LyB_LzF(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz, &
!$OMP   DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   DxTxx = (              &
     m3d_FDxF1(Txx,i,j,k) &
     m3d_FDxF2(Txx,i,j,k) &
     )/steph
   DxTyy = (              &
     m3d_FDxF1(Tyy,i,j,k) &
     m3d_FDxF2(Tyy,i,j,k) &
     )/steph
   DxTzz = (              &
     m3d_FDxF1(Tzz,i,j,k) &
     m3d_FDxF2(Tzz,i,j,k) &
     )/steph
   DxTxy = (              &
     m3d_FDxF1(Txy,i,j,k) &
     m3d_FDxF2(Txy,i,j,k) &
     )/steph
   DxTxz = (              &
     m3d_FDxF1(Txz,i,j,k) &
     m3d_FDxF2(Txz,i,j,k) &
     )/steph
   DxTyz = (              &
     m3d_FDxF1(Tyz,i,j,k) &
     m3d_FDxF2(Tyz,i,j,k) &
     )/steph
   DxVx =  (              &
     m3d_FDxF1(Vx,i,j,k)  &
     m3d_FDxF2(Vx,i,j,k)  &
     )/steph
   DxVy = (               &
     m3d_FDxF1(Vy,i,j,k)  &
     m3d_FDxF2(Vy,i,j,k)  &
     )/steph
   DxVz = (               &
     m3d_FDxF1(Vz,i,j,k)  &
     m3d_FDxF2(Vz,i,j,k)  &
     )/steph

   DyTxx = (              &
     m3d_FDyB1(Txx,i,j,k) &
     m3d_FDyB2(Txx,i,j,k) &
     )/steph
   DyTyy = (              &
     m3d_FDyB1(Tyy,i,j,k) &
     m3d_FDyB2(Tyy,i,j,k) &
     )/steph
   DyTzz = (              &
     m3d_FDyB1(Tzz,i,j,k) &
     m3d_FDyB2(Tzz,i,j,k) &
     )/steph
   DyTxy = (              &
     m3d_FDyB1(Txy,i,j,k) &
     m3d_FDyB2(Txy,i,j,k) &
     )/steph
   DyTxz = (              &
     m3d_FDyB1(Txz,i,j,k) &
     m3d_FDyB2(Txz,i,j,k) &
     )/steph
   DyTyz = (              &
     m3d_FDyB1(Tyz,i,j,k) &
     m3d_FDyB2(Tyz,i,j,k) &
     )/steph
   DyVx = (               &
     m3d_FDyB1(Vx,i,j,k)  &
     m3d_FDyB2(Vx,i,j,k)  &
     )/steph
   DyVy = (               &
     m3d_FDyB1(Vy,i,j,k)  &
     m3d_FDyB2(Vy,i,j,k)  &
     )/steph
   DyVz = (               &
     m3d_FDyB1(Vz,i,j,k)  &
     m3d_FDyB2(Vz,i,j,k)  &
     )/steph

   DzTxx = (              &
     m3d_FDzF1(Txx,i,j,k) &
     m3d_FDzF2(Txx,i,j,k) &
     )/steph
   DzTyy = (              &
     m3d_FDzF1(Tyy,i,j,k) &
     m3d_FDzF2(Tyy,i,j,k) &
     )/steph
   DzTzz = (              &
     m3d_FDzF1(Tzz,i,j,k) &
     m3d_FDzF2(Tzz,i,j,k) &
     )/steph
   DzTxy = (              &
     m3d_FDzF1(Txy,i,j,k) &
     m3d_FDzF2(Txy,i,j,k) &
     )/steph
   DzTxz = (              &
     m3d_FDzF1(Txz,i,j,k) &
     m3d_FDzF2(Txz,i,j,k) &
     )/steph
   DzTyz = (              &
     m3d_FDzF1(Tyz,i,j,k) &
     m3d_FDzF2(Tyz,i,j,k) &
     )/steph

#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzF1(Vx,i,j,k)  &
     m22_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m22_FDzF1(Vy,i,j,k)  &
     m22_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m22_FDzF1(Vz,i,j,k)  &
     m22_FDzF2(Vz,i,j,k)  &
     )/steph
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzF1(Vx,i,j,k)  &
     m24_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m24_FDzF1(Vy,i,j,k)  &
     m24_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m24_FDzF1(Vz,i,j,k)  &
     m24_FDzF2(Vz,i,j,k)  &
     )/steph
   else
#endif
   DzVx = (               &
     m3d_FDzF1(Vx,i,j,k)  &
     m3d_FDzF2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m3d_FDzF1(Vy,i,j,k)  &
     m3d_FDzF2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m3d_FDzF1(Vz,i,j,k)  &
     m3d_FDzF2(Vz,i,j,k)  &
     )/steph
#ifdef CondFreeVLOW
   end if
#endif

#ifndef CondFreeVHOC
#ifndef CondFreeCharac
if (freenode .and. k==nk2) then
   DzVx = matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz &
         +matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy = matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz &
         +matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz = matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz &
         +matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx=DzVx                          &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy=DzVy                          &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz=DzVz                          &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)
end if
#endif
#endif

   hVx(i,j,k)= rrho*(                          &
           xix*DxTxx + xiy*DxTxy + xiz*DxTxz   &
          +etx*DyTxx + ety*DyTxy + etz*DyTxz   &
          +ztx*DzTxx + zty*DzTxy + ztz*DzTxz )
   hVy(i,j,k)= rrho*(                          &
           xix*DxTxy + xiy*DxTyy + xiz*DxTyz   &
          +etx*DyTxy + ety*DyTyy + etz*DyTyz   &
          +ztx*DzTxy + zty*DzTyy + ztz*DzTyz )
   hVz(i,j,k)= rrho*(                          &
           xix*DxTxz + xiy*DxTyz + xiz*DxTzz   &
          +etx*DyTxz + ety*DyTyz + etz*DyTzz   &
          +ztx*DzTxz + zty*DzTyz + ztz*DzTzz )

   hTxx(i,j,k)= lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz &
               +lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz &
               +lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz
   hTyy(i,j,k)= lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz &
               +lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz &
               +lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz
   hTzz(i,j,k)= lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz &
               +lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz &
               +lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz
   hTxy(i,j,k)= miu *(              &
                xiy*DxVx + xix*DxVy &
               +ety*DyVx + etx*DyVy &
               +zty*DzVx + ztx*DzVy &
               )
   hTxz(i,j,k)= miu *(              &
                xiz*DxVx + xix*DxVz &
               +etz*DyVx + etx*DyVz &
               +ztz*DzVx + ztx*DzVz &
               )
   hTyz(i,j,k)= miu *(              &
                xiz*DxVy + xiy*DxVz &
               +etz*DyVy + ety*DyVz &
               +ztz*DzVy + zty*DzVz &
               )
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxF_LyB_LzF

subroutine LxB_LyF_LzB(I1,I2,J1,J2,K1,K2)
integer,intent(in) :: I1,I2,J1,J2,K1,K2
integer :: i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz, &
!$OMP   DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz )
do k=K1,K2
do j=J1,J2
do i=I1,I2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   DxTxx = (              &
     m3d_FDxB1(Txx,i,j,k) &
     m3d_FDxB2(Txx,i,j,k) &
     )/steph
   DxTyy = (              &
     m3d_FDxB1(Tyy,i,j,k) &
     m3d_FDxB2(Tyy,i,j,k) &
     )/steph
   DxTzz = (              &
     m3d_FDxB1(Tzz,i,j,k) &
     m3d_FDxB2(Tzz,i,j,k) &
     )/steph
   DxTxy = (              &
     m3d_FDxB1(Txy,i,j,k) &
     m3d_FDxB2(Txy,i,j,k) &
     )/steph
   DxTxz = (              &
     m3d_FDxB1(Txz,i,j,k) &
     m3d_FDxB2(Txz,i,j,k) &
     )/steph
   DxTyz = (              &
     m3d_FDxB1(Tyz,i,j,k) &
     m3d_FDxB2(Tyz,i,j,k) &
     )/steph
   DxVx =  (              &
     m3d_FDxB1(Vx,i,j,k)  &
     m3d_FDxB2(Vx,i,j,k)  &
     )/steph
   DxVy = (               &
     m3d_FDxB1(Vy,i,j,k)  &
     m3d_FDxB2(Vy,i,j,k)  &
     )/steph
   DxVz = (               &
     m3d_FDxB1(Vz,i,j,k)  &
     m3d_FDxB2(Vz,i,j,k)  &
     )/steph

   DyTxx = (              &
     m3d_FDyF1(Txx,i,j,k) &
     m3d_FDyF2(Txx,i,j,k) &
     )/steph
   DyTyy = (              &
     m3d_FDyF1(Tyy,i,j,k) &
     m3d_FDyF2(Tyy,i,j,k) &
     )/steph
   DyTzz = (              &
     m3d_FDyF1(Tzz,i,j,k) &
     m3d_FDyF2(Tzz,i,j,k) &
     )/steph
   DyTxy = (              &
     m3d_FDyF1(Txy,i,j,k) &
     m3d_FDyF2(Txy,i,j,k) &
     )/steph
   DyTxz = (              &
     m3d_FDyF1(Txz,i,j,k) &
     m3d_FDyF2(Txz,i,j,k) &
     )/steph
   DyTyz = (              &
     m3d_FDyF1(Tyz,i,j,k) &
     m3d_FDyF2(Tyz,i,j,k) &
     )/steph
   DyVx = (               &
     m3d_FDyF1(Vx,i,j,k)  &
     m3d_FDyF2(Vx,i,j,k)  &
     )/steph
   DyVy = (               &
     m3d_FDyF1(Vy,i,j,k)  &
     m3d_FDyF2(Vy,i,j,k)  &
     )/steph
   DyVz = (               &
     m3d_FDyF1(Vz,i,j,k)  &
     m3d_FDyF2(Vz,i,j,k)  &
     )/steph

   DzTxx = (              &
     m3d_FDzB1(Txx,i,j,k) &
     m3d_FDzB2(Txx,i,j,k) &
     )/steph
   DzTyy = (              &
     m3d_FDzB1(Tyy,i,j,k) &
     m3d_FDzB2(Tyy,i,j,k) &
     )/steph
   DzTzz = (              &
     m3d_FDzB1(Tzz,i,j,k) &
     m3d_FDzB2(Tzz,i,j,k) &
     )/steph
   DzTxy = (              &
     m3d_FDzB1(Txy,i,j,k) &
     m3d_FDzB2(Txy,i,j,k) &
     )/steph
   DzTxz = (              &
     m3d_FDzB1(Txz,i,j,k) &
     m3d_FDzB2(Txz,i,j,k) &
     )/steph
   DzTyz = (              &
     m3d_FDzB1(Tyz,i,j,k) &
     m3d_FDzB2(Tyz,i,j,k) &
     )/steph

#ifdef CondFreeVLOW
   if (freenode .and. k==nk2-1) then
   DzVx = (               &
     m22_FDzB1(Vx,i,j,k)  &
     m22_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m22_FDzB1(Vy,i,j,k)  &
     m22_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m22_FDzB1(Vz,i,j,k)  &
     m22_FDzB2(Vz,i,j,k)  &
     )/steph
   elseif (freenode .and. k==nk2-2) then
   DzVx = (               &
     m24_FDzB1(Vx,i,j,k)  &
     m24_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m24_FDzB1(Vy,i,j,k)  &
     m24_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m24_FDzB1(Vz,i,j,k)  &
     m24_FDzB2(Vz,i,j,k)  &
     )/steph
   else
#endif
   DzVx = (               &
     m3d_FDzB1(Vx,i,j,k)  &
     m3d_FDzB2(Vx,i,j,k)  &
     )/steph
   DzVy = (               &
     m3d_FDzB1(Vy,i,j,k)  &
     m3d_FDzB2(Vy,i,j,k)  &
     )/steph
   DzVz = (               &
     m3d_FDzB1(Vz,i,j,k)  &
     m3d_FDzB2(Vz,i,j,k)  &
     )/steph
#ifdef CondFreeVLOW
   end if
#endif

#ifndef CondFreeVHOC
#ifndef CondFreeCharac
if (freenode .and. k==nk2) then
   DzVx = matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz &
         +matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy = matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz &
         +matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz = matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz &
         +matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx=DzVx                          &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy=DzVy                          &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz=DzVz                          &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)
end if
#endif
#endif

   hVx(i,j,k)= rrho*(                          &
           xix*DxTxx + xiy*DxTxy + xiz*DxTxz   &
          +etx*DyTxx + ety*DyTxy + etz*DyTxz   &
          +ztx*DzTxx + zty*DzTxy + ztz*DzTxz )
   hVy(i,j,k)= rrho*(                          &
           xix*DxTxy + xiy*DxTyy + xiz*DxTyz   &
          +etx*DyTxy + ety*DyTyy + etz*DyTyz   &
          +ztx*DzTxy + zty*DzTyy + ztz*DzTyz )
   hVz(i,j,k)= rrho*(                          &
           xix*DxTxz + xiy*DxTyz + xiz*DxTzz   &
          +etx*DyTxz + ety*DyTyz + etz*DyTzz   &
          +ztx*DzTxz + zty*DzTyz + ztz*DzTzz )

   hTxx(i,j,k)= lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz &
               +lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz &
               +lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz
   hTyy(i,j,k)= lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz &
               +lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz &
               +lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz
   hTzz(i,j,k)= lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz &
               +lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz &
               +lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz
   hTxy(i,j,k)= miu *(              &
                xiy*DxVx + xix*DxVy &
               +ety*DyVx + etx*DyVy &
               +zty*DzVx + ztx*DzVy &
               )
   hTxz(i,j,k)= miu *(              &
                xiz*DxVx + xix*DxVz &
               +etz*DyVx + etx*DyVz &
               +ztz*DzVx + ztx*DzVz &
               )
   hTyz(i,j,k)= miu *(              &
                xiz*DxVy + xiy*DxVz &
               +etz*DyVy + ety*DyVz &
               +ztz*DzVy + zty*DzVz &
               )
end do
end do
end do
!$OMP END PARALLEL DO
end subroutine LxB_LyF_LzB

!*************************************************************************
!* use Traction Image method to calculate fd of stress compoents,        *
!* then assemble the right hand side to update velocities                *
!*************************************************************************
subroutine symmetric_cond
   Txx(:,:,nk2+1:nk2+LenFD)=Txx(:,:,nk2-1:nk2-LenFD:-1)
   Tyy(:,:,nk2+1:nk2+LenFD)=Tyy(:,:,nk2-1:nk2-LenFD:-1)
   Tzz(:,:,nk2+1:nk2+LenFD)=Tzz(:,:,nk2-1:nk2-LenFD:-1)
   Txy(:,:,nk2+1:nk2+LenFD)=Txy(:,:,nk2-1:nk2-LenFD:-1)
   Txz(:,:,nk2+1:nk2+LenFD)=Txz(:,:,nk2-1:nk2-LenFD:-1)
   Tyz(:,:,nk2+1:nk2+LenFD)=Tyz(:,:,nk2-1:nk2-LenFD:-1)
   Vx (:,:,nk2+1:nk2+LenFD)=Vx (:,:,nk2-1:nk2-LenFD:-1)
   Vy (:,:,nk2+1:nk2+LenFD)=Vy (:,:,nk2-1:nk2-LenFD:-1)
   Vz (:,:,nk2+1:nk2+LenFD)=Vz (:,:,nk2-1:nk2-LenFD:-1)
end subroutine symmetric_cond
subroutine LxF_LyF_LzF_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n

do n=1,LenFD ! from surface to inner domain
   k=nk2-n+1
do j=nj1,nj2
do i=ni1,ni2
   rrhojac=1.0/rho(i,j,k)/jac(i,j,k)
   !-- hVx --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txx(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txx(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txx(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TxSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TxSrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVx(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVy --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TySrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TySrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVy(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVz --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tzz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tzz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tzz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TzSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TzSrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVz(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
end do
end do
end do

end subroutine LxF_LyF_LzF_TIMG
subroutine LxB_LyB_LzB_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n

do n=1,LenFD ! from surface to inner domain
   k=nk2-n+1
do j=nj1,nj2
do i=ni1,ni2
   rrhojac=1.0/rho(i,j,k)/jac(i,j,k)
   !-- hVx --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txx(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txx(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txx(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TxSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TxSrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVx(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVy --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TySrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TySrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVy(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVz --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tzz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tzz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tzz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TzSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TzSrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVz(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
end do
end do
end do

end subroutine LxB_LyB_LzB_TIMG

subroutine LxF_LyF_LzB_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n

do n=1,LenFD ! from surface to inner domain
   k=nk2-n+1
do j=nj1,nj2
do i=ni1,ni2
   rrhojac=1.0/rho(i,j,k)/jac(i,j,k)
   !-- hVx --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txx(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txx(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txx(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TxSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TxSrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVx(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVy --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TySrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TySrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVy(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVz --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tzz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tzz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tzz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TzSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TzSrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVz(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
end do
end do
end do
end subroutine LxF_LyF_LzB_TIMG
subroutine LxB_LyB_LzF_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n

do n=1,LenFD ! from surface to inner domain
   k=nk2-n+1
do j=nj1,nj2
do i=ni1,ni2
   rrhojac=1.0/rho(i,j,k)/jac(i,j,k)
   !-- hVx --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txx(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txx(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txx(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TxSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TxSrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVx(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVy --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TySrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TySrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVy(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVz --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tzz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tzz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tzz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TzSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TzSrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVz(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
end do
end do
end do

end subroutine LxB_LyB_LzF_TIMG

subroutine LxB_LyF_LzF_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n

do n=1,LenFD ! from surface to inner domain
   k=nk2-n+1
do j=nj1,nj2
do i=ni1,ni2
   rrhojac=1.0/rho(i,j,k)/jac(i,j,k)
   !-- hVx --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txx(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txx(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txx(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TxSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TxSrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVx(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVy --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TySrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TySrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVy(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVz --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tzz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tzz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tzz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TzSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TzSrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVz(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
end do
end do
end do

end subroutine LxB_LyF_LzF_TIMG
subroutine LxF_LyB_LzB_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n

do n=1,LenFD ! from surface to inner domain
   k=nk2-n+1
do j=nj1,nj2
do i=ni1,ni2
   rrhojac=1.0/rho(i,j,k)/jac(i,j,k)
   !-- hVx --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txx(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txx(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txx(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TxSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TxSrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVx(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVy --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TySrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TySrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVy(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVz --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tzz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tzz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tzz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TzSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TzSrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVz(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
end do
end do
end do

end subroutine LxF_LyB_LzB_TIMG

subroutine LxF_LyB_LzF_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n

do n=1,LenFD ! from surface to inner domain
   k=nk2-n+1
do j=nj1,nj2
do i=ni1,ni2
   rrhojac=1.0/rho(i,j,k)/jac(i,j,k)
   !-- hVx --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txx(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txx(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txx(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TxSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TxSrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVx(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVy --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TySrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TySrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVy(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVz --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tzz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tzz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tzz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TzSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TzSrc(i,j)
   DxTx = (             &
     vec_FD_F1(vecTx,0) &
     vec_FD_F2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_B1(vecTy,0) &
     vec_FD_B2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_F1(vecTz,0) &
     vec_FD_F2(vecTz,0) &
     )/steph
   hVz(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
end do
end do
end do

end subroutine LxF_LyB_LzF_TIMG
subroutine LxB_LyF_LzB_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n

do n=1,LenFD ! from surface to inner domain
   k=nk2-n+1
do j=nj1,nj2
do i=ni1,ni2
   rrhojac=1.0/rho(i,j,k)/jac(i,j,k)
   !-- hVx --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txx(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txx(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txx(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TxSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TxSrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVx(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVy --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txy(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyy(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txy(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyy(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txy(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyy(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TySrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TySrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVy(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
   !-- hVz --
   vecTx(-LenFD:LenFD)=jac(i-LenFD:i+LenFD,j,k)*(                   &
               xi_x(i-LenFD:i+LenFD,j,k)*Txz(i-LenFD:i+LenFD,j,k)   &
              +xi_y(i-LenFD:i+LenFD,j,k)*Tyz(i-LenFD:i+LenFD,j,k)   &
              +xi_z(i-LenFD:i+LenFD,j,k)*Tzz(i-LenFD:i+LenFD,j,k)   &
              )
   vecTy(-LenFD:LenFD)=jac(i,j-LenFD:j+LenFD,k)*(                   &
               eta_x(i,j-LenFD:j+LenFD,k)*Txz(i,j-LenFD:j+LenFD,k)  &
              +eta_y(i,j-LenFD:j+LenFD,k)*Tyz(i,j-LenFD:j+LenFD,k)  &
              +eta_z(i,j-LenFD:j+LenFD,k)*Tzz(i,j-LenFD:j+LenFD,k)  &
              )
   vecTz(-LenFD:LenFD)=jac(i,j,k-LenFD:k+LenFD)*(                   &
               zeta_x(i,j,k-LenFD:k+LenFD)*Txz(i,j,k-LenFD:k+LenFD) &
              +zeta_y(i,j,k-LenFD:k+LenFD)*Tyz(i,j,k-LenFD:k+LenFD) &
              +zeta_z(i,j,k-LenFD:k+LenFD)*Tzz(i,j,k-LenFD:k+LenFD) &
              )
   vecTz(n:LenFD)=2.0_SP*TzSrc(i,j)-vecTz(n-2:n-2 -(LenFD-n):-1 )
   vecTz(n-1)=TzSrc(i,j)
   DxTx = (             &
     vec_FD_B1(vecTx,0) &
     vec_FD_B2(vecTx,0) &
     )/steph
   DyTy = (             &
     vec_FD_F1(vecTy,0) &
     vec_FD_F2(vecTy,0) &
     )/steph
   DzTz = (             &
     vec_FD_B1(vecTz,0) &
     vec_FD_B2(vecTz,0) &
     )/steph
   hVz(i,j,k)=( DxTx+DyTy+DzTz )*rrhojac
end do
end do
end do

end subroutine LxB_LyF_LzB_TIMG

!*************************************************************************
!* use Compact MacCormack scheme to calculate velocities fd with         *
!* respect to eta and assemble the right hand side to update stresses    *
!*************************************************************************
subroutine LxF_LyF_LzF_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz

#ifdef WithoutVHOC
   return
#endif
!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )/steph

   DzVx(n) = matVx2Vz(1,1,i,j)*DxVx(n) &
            +matVx2Vz(1,2,i,j)*DxVy(n) &
            +matVx2Vz(1,3,i,j)*DxVz(n) &
            +matVy2Vz(1,1,i,j)*DyVx(n) &
            +matVy2Vz(1,2,i,j)*DyVy(n) &
            +matVy2Vz(1,3,i,j)*DyVz(n)
   DzVy(n) = matVx2Vz(2,1,i,j)*DxVx(n) &
            +matVx2Vz(2,2,i,j)*DxVy(n) &
            +matVx2Vz(2,3,i,j)*DxVz(n) &
            +matVy2Vz(2,1,i,j)*DyVx(n) &
            +matVy2Vz(2,2,i,j)*DyVy(n) &
            +matVy2Vz(2,3,i,j)*DyVz(n)
   DzVz(n) = matVx2Vz(3,1,i,j)*DxVx(n) &
            +matVx2Vz(3,2,i,j)*DxVy(n) &
            +matVx2Vz(3,3,i,j)*DxVz(n) &
            +matVy2Vz(3,1,i,j)*DyVx(n) &
            +matVy2Vz(3,2,i,j)*DyVy(n) &
            +matVy2Vz(3,3,i,j)*DyVz(n)
   DzVx(n)=DzVx(n)                    &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy(n)=DzVy(n)                    &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz(n)=DzVz(n)                    &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)

   !-- j=ny-LenFD+1:ny-1 --
do n=LenFD,2,-1
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )/steph

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vx,i,j,k) &
     m3d_HOCzF2_RHS(Vx,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vy,i,j,k) &
     m3d_HOCzF2_RHS(Vy,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vz,i,j,k) &
     m3d_HOCzF2_RHS(Vz,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz

end do
do n=2,LenFD+1
   k=nk2-LenFD +n-1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   hTxx(i,j,k)= lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTyy(i,j,k)= lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTzz(i,j,k)= lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n)
   hTxy(i,j,k)= miu *(                    &
                xiy*DxVx(n) + xix*DxVy(n) &
               +ety*DyVx(n) + etx*DyVy(n) &
               +zty*DzVx(n) + ztx*DzVy(n) &
               )
   hTxz(i,j,k)= miu *(                    &
                xiz*DxVx(n) + xix*DxVz(n) &
               +etz*DyVx(n) + etx*DyVz(n) &
               +ztz*DzVx(n) + ztx*DzVz(n) &
               )
   hTyz(i,j,k)= miu *(                    &
                xiz*DxVy(n) + xiy*DxVz(n) &
               +etz*DyVy(n) + ety*DyVz(n) &
               +ztz*DzVy(n) + zty*DzVz(n) &
               )
end do
end do loop_xi
end do loop_eta
end subroutine LxF_LyF_LzF_VHOC
subroutine LxB_LyB_LzB_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz

#ifdef WithoutVHOC
   return
#endif
!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz-LenFD --
   n=1; k=nk2-LenFD
   DzVx(n) =  (          &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )/steph
   DzVy(n) = (           &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )/steph
   DzVz(n) = (           &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )/steph

   !-- k=nk2 --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )/steph

   DzVx(n) = matVx2Vz(1,1,i,j)*DxVx(n) &
            +matVx2Vz(1,2,i,j)*DxVy(n) &
            +matVx2Vz(1,3,i,j)*DxVz(n) &
            +matVy2Vz(1,1,i,j)*DyVx(n) &
            +matVy2Vz(1,2,i,j)*DyVy(n) &
            +matVy2Vz(1,3,i,j)*DyVz(n)
   DzVy(n) = matVx2Vz(2,1,i,j)*DxVx(n) &
            +matVx2Vz(2,2,i,j)*DxVy(n) &
            +matVx2Vz(2,3,i,j)*DxVz(n) &
            +matVy2Vz(2,1,i,j)*DyVx(n) &
            +matVy2Vz(2,2,i,j)*DyVy(n) &
            +matVy2Vz(2,3,i,j)*DyVz(n)
   DzVz(n) = matVx2Vz(3,1,i,j)*DxVx(n) &
            +matVx2Vz(3,2,i,j)*DxVy(n) &
            +matVx2Vz(3,3,i,j)*DxVz(n) &
            +matVy2Vz(3,1,i,j)*DyVx(n) &
            +matVy2Vz(3,2,i,j)*DyVy(n) &
            +matVy2Vz(3,3,i,j)*DyVz(n)
   DzVx(n)=DzVx(n)                    &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy(n)=DzVy(n)                    &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz(n)=DzVz(n)                    &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)

   !-- k=nk2-LenFD+1:nk2-1 --
do n=2,LenFD
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )/steph

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vx,i,j,k) &
     m3d_HOCzB2_RHS(Vx,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vy,i,j,k) &
     m3d_HOCzB2_RHS(Vy,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vz,i,j,k) &
     m3d_HOCzB2_RHS(Vz,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz
end do

do n=2,LenFD+1
   k=nk2-LenFD +n-1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   hTxx(i,j,k)= lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTyy(i,j,k)= lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTzz(i,j,k)= lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n)
   hTxy(i,j,k)= miu *(                    &
                xiy*DxVx(n) + xix*DxVy(n) &
               +ety*DyVx(n) + etx*DyVy(n) &
               +zty*DzVx(n) + ztx*DzVy(n) &
               )
   hTxz(i,j,k)= miu *(                    &
                xiz*DxVx(n) + xix*DxVz(n) &
               +etz*DyVx(n) + etx*DyVz(n) &
               +ztz*DzVx(n) + ztx*DzVz(n) &
               )
   hTyz(i,j,k)= miu *(                    &
                xiz*DxVy(n) + xiy*DxVz(n) &
               +etz*DyVy(n) + ety*DyVz(n) &
               +ztz*DzVy(n) + zty*DzVz(n) &
               )
end do
end do loop_xi
end do loop_eta
end subroutine LxB_LyB_LzB_VHOC

subroutine LxF_LyF_LzB_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz

#ifdef WithoutVHOC
   return
#endif
!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz-LenFD --
   n=1; k=nk2-LenFD
   DzVx(n) =  (          &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )/steph
   DzVy(n) = (           &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )/steph
   DzVz(n) = (           &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )/steph

   !-- k=nk2 --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )/steph

   DzVx(n) = matVx2Vz(1,1,i,j)*DxVx(n) &
            +matVx2Vz(1,2,i,j)*DxVy(n) &
            +matVx2Vz(1,3,i,j)*DxVz(n) &
            +matVy2Vz(1,1,i,j)*DyVx(n) &
            +matVy2Vz(1,2,i,j)*DyVy(n) &
            +matVy2Vz(1,3,i,j)*DyVz(n)
   DzVy(n) = matVx2Vz(2,1,i,j)*DxVx(n) &
            +matVx2Vz(2,2,i,j)*DxVy(n) &
            +matVx2Vz(2,3,i,j)*DxVz(n) &
            +matVy2Vz(2,1,i,j)*DyVx(n) &
            +matVy2Vz(2,2,i,j)*DyVy(n) &
            +matVy2Vz(2,3,i,j)*DyVz(n)
   DzVz(n) = matVx2Vz(3,1,i,j)*DxVx(n) &
            +matVx2Vz(3,2,i,j)*DxVy(n) &
            +matVx2Vz(3,3,i,j)*DxVz(n) &
            +matVy2Vz(3,1,i,j)*DyVx(n) &
            +matVy2Vz(3,2,i,j)*DyVy(n) &
            +matVy2Vz(3,3,i,j)*DyVz(n)
   DzVx(n)=DzVx(n)                    &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy(n)=DzVy(n)                    &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz(n)=DzVz(n)                    &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)

   !-- k=nk2-LenFD+1:nk2-1 --
do n=2,LenFD
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )/steph

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vx,i,j,k) &
     m3d_HOCzB2_RHS(Vx,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vy,i,j,k) &
     m3d_HOCzB2_RHS(Vy,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vz,i,j,k) &
     m3d_HOCzB2_RHS(Vz,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz
end do

do n=2,LenFD+1
   k=nk2-LenFD +n-1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   hTxx(i,j,k)= lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTyy(i,j,k)= lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTzz(i,j,k)= lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n)
   hTxy(i,j,k)= miu *(                    &
                xiy*DxVx(n) + xix*DxVy(n) &
               +ety*DyVx(n) + etx*DyVy(n) &
               +zty*DzVx(n) + ztx*DzVy(n) &
               )
   hTxz(i,j,k)= miu *(                    &
                xiz*DxVx(n) + xix*DxVz(n) &
               +etz*DyVx(n) + etx*DyVz(n) &
               +ztz*DzVx(n) + ztx*DzVz(n) &
               )
   hTyz(i,j,k)= miu *(                    &
                xiz*DxVy(n) + xiy*DxVz(n) &
               +etz*DyVy(n) + ety*DyVz(n) &
               +ztz*DzVy(n) + zty*DzVz(n) &
               )
end do
end do loop_xi
end do loop_eta
end subroutine LxF_LyF_LzB_VHOC
subroutine LxB_LyB_LzF_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
#ifdef WithoutVHOC
   return
#endif
!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )/steph

   DzVx(n) = matVx2Vz(1,1,i,j)*DxVx(n) &
            +matVx2Vz(1,2,i,j)*DxVy(n) &
            +matVx2Vz(1,3,i,j)*DxVz(n) &
            +matVy2Vz(1,1,i,j)*DyVx(n) &
            +matVy2Vz(1,2,i,j)*DyVy(n) &
            +matVy2Vz(1,3,i,j)*DyVz(n)
   DzVy(n) = matVx2Vz(2,1,i,j)*DxVx(n) &
            +matVx2Vz(2,2,i,j)*DxVy(n) &
            +matVx2Vz(2,3,i,j)*DxVz(n) &
            +matVy2Vz(2,1,i,j)*DyVx(n) &
            +matVy2Vz(2,2,i,j)*DyVy(n) &
            +matVy2Vz(2,3,i,j)*DyVz(n)
   DzVz(n) = matVx2Vz(3,1,i,j)*DxVx(n) &
            +matVx2Vz(3,2,i,j)*DxVy(n) &
            +matVx2Vz(3,3,i,j)*DxVz(n) &
            +matVy2Vz(3,1,i,j)*DyVx(n) &
            +matVy2Vz(3,2,i,j)*DyVy(n) &
            +matVy2Vz(3,3,i,j)*DyVz(n)
   DzVx(n)=DzVx(n)                    &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy(n)=DzVy(n)                    &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz(n)=DzVz(n)                    &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)

   !-- j=ny-LenFD+1:ny-1 --
do n=LenFD,2,-1
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )/steph

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vx,i,j,k) &
     m3d_HOCzF2_RHS(Vx,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vy,i,j,k) &
     m3d_HOCzF2_RHS(Vy,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vz,i,j,k) &
     m3d_HOCzF2_RHS(Vz,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz

end do
do n=2,LenFD+1
   k=nk2-LenFD +n-1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   hTxx(i,j,k)= lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTyy(i,j,k)= lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTzz(i,j,k)= lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n)
   hTxy(i,j,k)= miu *(                    &
                xiy*DxVx(n) + xix*DxVy(n) &
               +ety*DyVx(n) + etx*DyVy(n) &
               +zty*DzVx(n) + ztx*DzVy(n) &
               )
   hTxz(i,j,k)= miu *(                    &
                xiz*DxVx(n) + xix*DxVz(n) &
               +etz*DyVx(n) + etx*DyVz(n) &
               +ztz*DzVx(n) + ztx*DzVz(n) &
               )
   hTyz(i,j,k)= miu *(                    &
                xiz*DxVy(n) + xiy*DxVz(n) &
               +etz*DyVy(n) + ety*DyVz(n) &
               +ztz*DzVy(n) + zty*DzVz(n) &
               )
end do
end do loop_xi
end do loop_eta
end subroutine LxB_LyB_LzF_VHOC

subroutine LxB_LyF_LzF_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz

#ifdef WithoutVHOC
   return
#endif
!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )/steph

   DzVx(n) = matVx2Vz(1,1,i,j)*DxVx(n) &
            +matVx2Vz(1,2,i,j)*DxVy(n) &
            +matVx2Vz(1,3,i,j)*DxVz(n) &
            +matVy2Vz(1,1,i,j)*DyVx(n) &
            +matVy2Vz(1,2,i,j)*DyVy(n) &
            +matVy2Vz(1,3,i,j)*DyVz(n)
   DzVy(n) = matVx2Vz(2,1,i,j)*DxVx(n) &
            +matVx2Vz(2,2,i,j)*DxVy(n) &
            +matVx2Vz(2,3,i,j)*DxVz(n) &
            +matVy2Vz(2,1,i,j)*DyVx(n) &
            +matVy2Vz(2,2,i,j)*DyVy(n) &
            +matVy2Vz(2,3,i,j)*DyVz(n)
   DzVz(n) = matVx2Vz(3,1,i,j)*DxVx(n) &
            +matVx2Vz(3,2,i,j)*DxVy(n) &
            +matVx2Vz(3,3,i,j)*DxVz(n) &
            +matVy2Vz(3,1,i,j)*DyVx(n) &
            +matVy2Vz(3,2,i,j)*DyVy(n) &
            +matVy2Vz(3,3,i,j)*DyVz(n)
   DzVx(n)=DzVx(n)                    &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy(n)=DzVy(n)                    &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz(n)=DzVz(n)                    &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)

   !-- j=ny-LenFD+1:ny-1 --
do n=LenFD,2,-1
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )/steph

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vx,i,j,k) &
     m3d_HOCzF2_RHS(Vx,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vy,i,j,k) &
     m3d_HOCzF2_RHS(Vy,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vz,i,j,k) &
     m3d_HOCzF2_RHS(Vz,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz

end do
do n=2,LenFD+1
   k=nk2-LenFD +n-1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   hTxx(i,j,k)= lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTyy(i,j,k)= lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTzz(i,j,k)= lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n)
   hTxy(i,j,k)= miu *(                    &
                xiy*DxVx(n) + xix*DxVy(n) &
               +ety*DyVx(n) + etx*DyVy(n) &
               +zty*DzVx(n) + ztx*DzVy(n) &
               )
   hTxz(i,j,k)= miu *(                    &
                xiz*DxVx(n) + xix*DxVz(n) &
               +etz*DyVx(n) + etx*DyVz(n) &
               +ztz*DzVx(n) + ztx*DzVz(n) &
               )
   hTyz(i,j,k)= miu *(                    &
                xiz*DxVy(n) + xiy*DxVz(n) &
               +etz*DyVy(n) + ety*DyVz(n) &
               +ztz*DzVy(n) + zty*DzVz(n) &
               )
end do
end do loop_xi
end do loop_eta
end subroutine LxB_LyF_LzF_VHOC
subroutine LxF_LyB_LzB_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz

#ifdef WithoutVHOC
   return
#endif
!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz-LenFD --
   n=1; k=nk2-LenFD
   DzVx(n) =  (          &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )/steph
   DzVy(n) = (           &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )/steph
   DzVz(n) = (           &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )/steph

   !-- k=nk2 --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )/steph

   DzVx(n) = matVx2Vz(1,1,i,j)*DxVx(n) &
            +matVx2Vz(1,2,i,j)*DxVy(n) &
            +matVx2Vz(1,3,i,j)*DxVz(n) &
            +matVy2Vz(1,1,i,j)*DyVx(n) &
            +matVy2Vz(1,2,i,j)*DyVy(n) &
            +matVy2Vz(1,3,i,j)*DyVz(n)
   DzVy(n) = matVx2Vz(2,1,i,j)*DxVx(n) &
            +matVx2Vz(2,2,i,j)*DxVy(n) &
            +matVx2Vz(2,3,i,j)*DxVz(n) &
            +matVy2Vz(2,1,i,j)*DyVx(n) &
            +matVy2Vz(2,2,i,j)*DyVy(n) &
            +matVy2Vz(2,3,i,j)*DyVz(n)
   DzVz(n) = matVx2Vz(3,1,i,j)*DxVx(n) &
            +matVx2Vz(3,2,i,j)*DxVy(n) &
            +matVx2Vz(3,3,i,j)*DxVz(n) &
            +matVy2Vz(3,1,i,j)*DyVx(n) &
            +matVy2Vz(3,2,i,j)*DyVy(n) &
            +matVy2Vz(3,3,i,j)*DyVz(n)
   DzVx(n)=DzVx(n)                    &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy(n)=DzVy(n)                    &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz(n)=DzVz(n)                    &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)

   !-- k=nk2-LenFD+1:nk2-1 --
do n=2,LenFD
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )/steph

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vx,i,j,k) &
     m3d_HOCzB2_RHS(Vx,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vy,i,j,k) &
     m3d_HOCzB2_RHS(Vy,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vz,i,j,k) &
     m3d_HOCzB2_RHS(Vz,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz
end do

do n=2,LenFD+1
   k=nk2-LenFD +n-1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   hTxx(i,j,k)= lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTyy(i,j,k)= lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTzz(i,j,k)= lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n)
   hTxy(i,j,k)= miu *(                    &
                xiy*DxVx(n) + xix*DxVy(n) &
               +ety*DyVx(n) + etx*DyVy(n) &
               +zty*DzVx(n) + ztx*DzVy(n) &
               )
   hTxz(i,j,k)= miu *(                    &
                xiz*DxVx(n) + xix*DxVz(n) &
               +etz*DyVx(n) + etx*DyVz(n) &
               +ztz*DzVx(n) + ztx*DzVz(n) &
               )
   hTyz(i,j,k)= miu *(                    &
                xiz*DxVy(n) + xiy*DxVz(n) &
               +etz*DyVy(n) + ety*DyVz(n) &
               +ztz*DzVy(n) + zty*DzVz(n) &
               )
end do
end do loop_xi
end do loop_eta
end subroutine LxF_LyB_LzB_VHOC

subroutine LxF_LyB_LzF_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz

#ifdef WithoutVHOC
   return
#endif
!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )/steph

   DzVx(n) = matVx2Vz(1,1,i,j)*DxVx(n) &
            +matVx2Vz(1,2,i,j)*DxVy(n) &
            +matVx2Vz(1,3,i,j)*DxVz(n) &
            +matVy2Vz(1,1,i,j)*DyVx(n) &
            +matVy2Vz(1,2,i,j)*DyVy(n) &
            +matVy2Vz(1,3,i,j)*DyVz(n)
   DzVy(n) = matVx2Vz(2,1,i,j)*DxVx(n) &
            +matVx2Vz(2,2,i,j)*DxVy(n) &
            +matVx2Vz(2,3,i,j)*DxVz(n) &
            +matVy2Vz(2,1,i,j)*DyVx(n) &
            +matVy2Vz(2,2,i,j)*DyVy(n) &
            +matVy2Vz(2,3,i,j)*DyVz(n)
   DzVz(n) = matVx2Vz(3,1,i,j)*DxVx(n) &
            +matVx2Vz(3,2,i,j)*DxVy(n) &
            +matVx2Vz(3,3,i,j)*DxVz(n) &
            +matVy2Vz(3,1,i,j)*DyVx(n) &
            +matVy2Vz(3,2,i,j)*DyVy(n) &
            +matVy2Vz(3,3,i,j)*DyVz(n)
   DzVx(n)=DzVx(n)                    &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy(n)=DzVy(n)                    &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz(n)=DzVz(n)                    &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)

   !-- j=ny-LenFD+1:ny-1 --
do n=LenFD,2,-1
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxF1(Vx,i,j,k) &
     m3d_FDxF2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxF1(Vy,i,j,k) &
     m3d_FDxF2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxF1(Vz,i,j,k) &
     m3d_FDxF2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyB1(Vx,i,j,k) &
     m3d_FDyB2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyB1(Vy,i,j,k) &
     m3d_FDyB2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyB1(Vz,i,j,k) &
     m3d_FDyB2(Vz,i,j,k) &
     )/steph

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vx,i,j,k) &
     m3d_HOCzF2_RHS(Vx,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vy,i,j,k) &
     m3d_HOCzF2_RHS(Vy,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzF1_RHS(Vz,i,j,k) &
     m3d_HOCzF2_RHS(Vz,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_F_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz

end do
do n=2,LenFD+1
   k=nk2-LenFD +n-1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   hTxx(i,j,k)= lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTyy(i,j,k)= lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTzz(i,j,k)= lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n)
   hTxy(i,j,k)= miu *(                    &
                xiy*DxVx(n) + xix*DxVy(n) &
               +ety*DyVx(n) + etx*DyVy(n) &
               +zty*DzVx(n) + ztx*DzVy(n) &
               )
   hTxz(i,j,k)= miu *(                    &
                xiz*DxVx(n) + xix*DxVz(n) &
               +etz*DyVx(n) + etx*DyVz(n) &
               +ztz*DzVx(n) + ztx*DzVz(n) &
               )
   hTyz(i,j,k)= miu *(                    &
                xiz*DxVy(n) + xiy*DxVz(n) &
               +etz*DyVy(n) + ety*DyVz(n) &
               +ztz*DzVy(n) + zty*DzVz(n) &
               )
end do
end do loop_xi
end do loop_eta
end subroutine LxF_LyB_LzF_VHOC
subroutine LxB_LyF_LzB_VHOC
integer :: i,j,k,n
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz

#ifdef WithoutVHOC
   return
#endif
!-- velocity fd --
loop_eta: do j=nj1,nj2
loop_xi:  do i=ni1,ni2
   !-- k=nz-LenFD --
   n=1; k=nk2-LenFD
   DzVx(n) =  (          &
     m3d_FDzB1(Vx,i,j,k) &
     m3d_FDzB2(Vx,i,j,k) &
     )/steph
   DzVy(n) = (           &
     m3d_FDzB1(Vy,i,j,k) &
     m3d_FDzB2(Vy,i,j,k) &
     )/steph
   DzVz(n) = (           &
     m3d_FDzB1(Vz,i,j,k) &
     m3d_FDzB2(Vz,i,j,k) &
     )/steph

   !-- k=nk2 --
   n=LenFD+1; k=nk2
   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )/steph

   DzVx(n) = matVx2Vz(1,1,i,j)*DxVx(n) &
            +matVx2Vz(1,2,i,j)*DxVy(n) &
            +matVx2Vz(1,3,i,j)*DxVz(n) &
            +matVy2Vz(1,1,i,j)*DyVx(n) &
            +matVy2Vz(1,2,i,j)*DyVy(n) &
            +matVy2Vz(1,3,i,j)*DyVz(n)
   DzVy(n) = matVx2Vz(2,1,i,j)*DxVx(n) &
            +matVx2Vz(2,2,i,j)*DxVy(n) &
            +matVx2Vz(2,3,i,j)*DxVz(n) &
            +matVy2Vz(2,1,i,j)*DyVx(n) &
            +matVy2Vz(2,2,i,j)*DyVy(n) &
            +matVy2Vz(2,3,i,j)*DyVz(n)
   DzVz(n) = matVx2Vz(3,1,i,j)*DxVx(n) &
            +matVx2Vz(3,2,i,j)*DxVy(n) &
            +matVx2Vz(3,3,i,j)*DxVz(n) &
            +matVy2Vz(3,1,i,j)*DyVx(n) &
            +matVy2Vz(3,2,i,j)*DyVy(n) &
            +matVy2Vz(3,3,i,j)*DyVz(n)
   DzVx(n)=DzVx(n)                    &
        + matF2Vz(1,1,i,j)*VxSrc(i,j) &
        + matF2Vz(1,2,i,j)*VySrc(i,j) &
        + matF2Vz(1,3,i,j)*VzSrc(i,j)
   DzVy(n)=DzVy(n)                    &
        + matF2Vz(2,1,i,j)*VxSrc(i,j) &
        + matF2Vz(2,2,i,j)*VySrc(i,j) &
        + matF2Vz(2,3,i,j)*VzSrc(i,j)
   DzVz(n)=DzVz(n)                    &
        + matF2Vz(3,1,i,j)*VxSrc(i,j) &
        + matF2Vz(3,2,i,j)*VySrc(i,j) &
        + matF2Vz(3,3,i,j)*VzSrc(i,j)

   !-- k=nk2-LenFD+1:nk2-1 --
do n=2,LenFD
   k=nk2-LenFD +n-1

   DxVx(n) =  (          &
     m3d_FDxB1(Vx,i,j,k) &
     m3d_FDxB2(Vx,i,j,k) &
     )/steph
   DxVy(n) = (           &
     m3d_FDxB1(Vy,i,j,k) &
     m3d_FDxB2(Vy,i,j,k) &
     )/steph
   DxVz(n) = (           &
     m3d_FDxB1(Vz,i,j,k) &
     m3d_FDxB2(Vz,i,j,k) &
     )/steph

   DyVx(n) =  (          &
     m3d_FDyF1(Vx,i,j,k) &
     m3d_FDyF2(Vx,i,j,k) &
     )/steph
   DyVy(n) = (           &
     m3d_FDyF1(Vy,i,j,k) &
     m3d_FDyF2(Vy,i,j,k) &
     )/steph
   DyVz(n) = (           &
     m3d_FDyF1(Vz,i,j,k) &
     m3d_FDyF2(Vz,i,j,k) &
     )/steph

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vx,i,j,k) &
     m3d_HOCzB2_RHS(Vx,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVx,n)    &
     )
   DzVx(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vy,i,j,k) &
     m3d_HOCzB2_RHS(Vy,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVy,n)    &
     )
   DzVy(n)=rhs_Dz - lhs_Dz

   rhs_Dz= (                  &
     m3d_HOCzB1_RHS(Vz,i,j,k) &
     m3d_HOCzB2_RHS(Vz,i,j,k) &
     )/steph
   lhs_Dz= (                  &
     vec_HOC_B_LHS(DzVz,n)    &
     )
   DzVz(n)=rhs_Dz - lhs_Dz
end do

do n=2,LenFD+1
   k=nk2-LenFD +n-1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   hTxx(i,j,k)= lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTyy(i,j,k)= lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n)
   hTzz(i,j,k)= lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) &
               +lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) &
               +lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n)
   hTxy(i,j,k)= miu *(                    &
                xiy*DxVx(n) + xix*DxVy(n) &
               +ety*DyVx(n) + etx*DyVy(n) &
               +zty*DzVx(n) + ztx*DzVy(n) &
               )
   hTxz(i,j,k)= miu *(                    &
                xiz*DxVx(n) + xix*DxVz(n) &
               +etz*DyVx(n) + etx*DyVz(n) &
               +ztz*DzVx(n) + ztx*DzVz(n) &
               )
   hTyz(i,j,k)= miu *(                    &
                xiz*DxVy(n) + xiy*DxVz(n) &
               +etz*DyVy(n) + ety*DyVz(n) &
               +ztz*DzVy(n) + zty*DzVz(n) &
               )
end do
end do loop_xi
end do loop_eta
end subroutine LxB_LyF_LzB_VHOC

subroutine free_charac_flat
integer i,j,k
real(SP) v1,v2,v3,t11,t22,t33,t12,t13,t23
real(SP) lam,miu,lam2mu,rrho,f1,f2,fct
real(SP) :: Tx,Ty,Tz

if (.not. freenode) return

   k=nk2
do j=nj1,nj2
do i=ni1,ni2

   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=rho(i,j,k)
   f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu
   v1 = Vx (i,j,k)
   v2 = Vy (i,j,k)
   v3 = Vz (i,j,k)
   t11= Txx(i,j,k)
   t22= Tyy(i,j,k)
   t33= Tzz(i,j,k)
   t12= Txy(i,j,k)
   t13= Txz(i,j,k)
   t23= Tyz(i,j,k)
   Tx=-TxSrc(i,j)
   Ty=-TySrc(i,j)
   Tz= TzSrc(i,j)
   Vx (i,j,k)=v1-(t13+Tx)/f2
   Vy (i,j,k)=v2-(t23+Ty)/f2
   Vz (i,j,k)=v3-(t33-Tz)/f1
   Txx(i,j,k)=t11-fct*(t33-Tz)
   Tyy(i,j,k)=t22-fct*(t33-Tz)
   Tzz(i,j,k)=Tz
   Txz(i,j,k)=Tx
   Tyz(i,j,k)=Ty
end do
end do
end subroutine free_charac_flat
subroutine free_charac
integer i,j,k

real(SP),dimension(SEIS_GEO) :: vecV1,vecV2
real(SP),dimension(SEIS_GEO,SEIS_GEO) :: matT1,matT2
real(SP),dimension(SEIS_GEO,SEIS_GEO) :: convD

real(SP),dimension(SEIS_GEO) :: e1,e2,e3
real(SP) :: lam,miu,lam2mu,rrho,f1,f2,fct
real(SP) :: Tx,Ty,Tz

if (.not. freenode) return

   k=nk2
do j=nj1,nj2
do i=ni1,ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=rho(i,j,k)
   f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu
   Tx=TxSrc(i,j); Ty=TySrc(i,j); Tz=TzSrc(i,j)

   e3=(/ zeta_x(i,j,k),zeta_y(i,j,k),zeta_z(i,j,k) /)
   e3=e3/sqrt(dot_product(e3,e3))

   call grid_covariant(i,j,k,2,e2)
   e2=e2/sqrt(dot_product(e2,e2))

   call times_product(e2,e3,e1)
   e1=e1/sqrt(dot_product(e1,e1))

   convD(1,:)=e1; convD(2,:)=e2; convD(3,:)=e3

   vecV1 =matmul(convD,(/Tx,Ty,Tz/))
   Tx=vecV1(1); Ty=vecV1(2); Tz=vecV1(3)

   vecV1=(/ Vx(i,j,k),Vy(i,j,k),Vz(i,j,k) /)
   vecV2=matmul(convD,vecV1)

   matT1(1,1)=Txx(i,j,k); matT1(2,2)=Tyy(i,j,k); matT1(3,3)=Tzz(i,j,k)
   matT1(1,2)=Txy(i,j,k); matT1(2,1)=Txy(i,j,k)
   matT1(1,3)=Txz(i,j,k); matT1(3,1)=Txz(i,j,k)
   matT1(2,3)=Tyz(i,j,k); matT1(3,2)=Tyz(i,j,k)
   matT2=matmul(matmul(convD,matT1),transpose(convD))

   vecV2(1)=vecV2(1)-(matT2(1,3)+Tx)/f2
   vecV2(2)=vecV2(2)-(matT2(2,3)+Ty)/f2
   vecV2(3)=vecV2(3)-(matT2(3,3)-Tz)/f1
   matT2(1,1)=matT2(1,1)-fct*(matT2(3,3)-Tz)
   matT2(2,2)=matT2(2,2)-fct*(matT2(3,3)-Tz)
   matT2(3,3)=Tz
   matT2(1,3)=Tx; matT2(3,1)=matT2(1,3)
   matT2(2,3)=Ty; matT2(3,2)=matT2(2,3)

   vecV1=matmul(transpose(convD),vecV2)
   matT1=matmul( matmul(transpose(convD),matT2),convD)

   Vx (i,j,k)=vecV1(1)
   Vy (i,j,k)=vecV1(2)
   Vz (i,j,k)=vecV1(3)
   Txx(i,j,k)=matT1(1,1)
   Tyy(i,j,k)=matT1(2,2)
   Tzz(i,j,k)=matT1(3,3)
   Txy(i,j,k)=matT1(1,2)
   Txz(i,j,k)=matT1(1,3)
   Tyz(i,j,k)=matT1(2,3)
end do
end do
end subroutine free_charac

subroutine free_extrap
integer i,j,k

if (.not. freenode) return

do k=nk2+1,nz2
do j=nj1,nj2
do i=ni1,ni2
   Vx (i,j,k)=4.0*Vx (i,j,k-1)-6.0*Vx (i,j,k-2)+4.0*Vx (i,j,k-3)-Vx (i,j,k-4)
   Vy (i,j,k)=4.0*Vy (i,j,k-1)-6.0*Vy (i,j,k-2)+4.0*Vy (i,j,k-3)-Vy (i,j,k-4)
   Vz (i,j,k)=4.0*Vz (i,j,k-1)-6.0*Vz (i,j,k-2)+4.0*Vz (i,j,k-3)-Vz (i,j,k-4)
   Txx(i,j,k)=4.0*Txx(i,j,k-1)-6.0*Txx(i,j,k-2)+4.0*Txx(i,j,k-3)-Txx(i,j,k-4)
   Tzz(i,j,k)=4.0*Tzz(i,j,k-1)-6.0*Tzz(i,j,k-2)+4.0*Tzz(i,j,k-3)-Tzz(i,j,k-4)
   Tyy(i,j,k)=4.0*Tyy(i,j,k-1)-6.0*Tyy(i,j,k-2)+4.0*Tyy(i,j,k-3)-Tyy(i,j,k-4)
   Txy(i,j,k)=4.0*Txy(i,j,k-1)-6.0*Txy(i,j,k-2)+4.0*Txy(i,j,k-3)-Txy(i,j,k-4)
   Txz(i,j,k)=4.0*Txz(i,j,k-1)-6.0*Txz(i,j,k-2)+4.0*Txz(i,j,k-3)-Txz(i,j,k-4)
   Tyz(i,j,k)=4.0*Tyz(i,j,k-1)-6.0*Tyz(i,j,k-2)+4.0*Tyz(i,j,k-3)-Tyz(i,j,k-4)
end do
end do
end do
end subroutine free_extrap

!--------------------------------------------------------------------}
subroutine macdrp_mesg_init
  call mesg_init_LxF
  call mesg_init_LxB
  call mesg_init_LyF
  call mesg_init_LyB
  call mesg_init_LzF
  call mesg_init_LzB
end subroutine macdrp_mesg_init

subroutine mesg_init_LxF
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
! --- LxF ------------------------------------------------------------------
! to X1
call MPI_SEND_INIT(Txx(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1131,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1132,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1133,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1134,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1135,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1136,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1137,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1138,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1139,SWMPI_COMM,s9,ierr)
! from X2
call MPI_RECV_INIT(Txx(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1131,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1132,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1133,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1134,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1135,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1136,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1137,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1138,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1139,SWMPI_COMM,r9,ierr)
! put into array
reqXF(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

!to X2
call MPI_SEND_INIT(Txx(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1211,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1212,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1213,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1214,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1215,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1219,SWMPI_COMM,s9,ierr)
!from X1
call MPI_RECV_INIT(Txx(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1211,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1212,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1213,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1214,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1215,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1219,SWMPI_COMM,r9,ierr)
! put into array
reqXF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
end subroutine mesg_init_LxF

subroutine mesg_init_LxB
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
! --- LxB ------------------------------------------------------------------
! to X1
call MPI_SEND_INIT(Txx(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1111,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1112,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1113,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1114,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1115,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1116,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1117,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1118,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1119,SWMPI_COMM,s9,ierr)
! from X2
call MPI_RECV_INIT(Txx(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1111,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1112,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1113,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1114,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1115,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1116,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1117,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1118,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1119,SWMPI_COMM,r9,ierr)
! put into array
reqXB(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

! to X2
call MPI_SEND_INIT(Txx(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1231,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1232,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1233,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1234,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1235,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1239,SWMPI_COMM,s9,ierr)
! from  X1
call MPI_RECV_INIT(Txx(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1231,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1232,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1233,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1234,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1235,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1239,SWMPI_COMM,r9,ierr)
! put into array
reqXB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
end subroutine mesg_init_LxB

subroutine mesg_init_LyF
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
! --- LyF ------------------------------------------------------------------
! to Y1
call MPI_SEND_INIT(Txx(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2131,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2132,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2133,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2134,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2135,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2136,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2137,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2138,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2139,SWMPI_COMM,s9,ierr)
! from Y2
call MPI_RECV_INIT(Txx(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2131,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2132,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2133,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2134,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2135,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2136,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2137,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2138,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2139,SWMPI_COMM,r9,ierr)
! put into array
reqYF(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
!-----

! to Y2
call MPI_SEND_INIT(Txx(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2211,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2212,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2213,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2214,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2215,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2219,SWMPI_COMM,s9,ierr)
! from Y1
call MPI_RECV_INIT(Txx(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2211,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2212,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2213,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2214,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2215,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2219,SWMPI_COMM,r9,ierr)
! put into array
reqYF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
end subroutine mesg_init_LyF

subroutine mesg_init_LyB
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
! --- LyB ------------------------------------------------------------------
! to Y1
call MPI_SEND_INIT(Txx(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2111,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2112,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2113,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2114,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2115,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2116,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2117,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2118,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2119,SWMPI_COMM,s9,ierr)
! from Y2
call MPI_RECV_INIT(Txx(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2111,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2112,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2113,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2114,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2115,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2116,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2117,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2118,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2119,SWMPI_COMM,r9,ierr)
! put into array
reqYB(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

! to Y2
call MPI_SEND_INIT(Txx(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2231,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2232,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2233,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2234,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2235,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2239,SWMPI_COMM,s9,ierr)
! from Y1
call MPI_RECV_INIT(Txx(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2231,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2232,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2233,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2234,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2235,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2239,SWMPI_COMM,r9,ierr)
! put into array
reqYB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
end subroutine mesg_init_LyB


subroutine mesg_init_LzF
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
! --- LzF ------------------------------------------------------------------
!to Z1
call MPI_SEND_INIT(Txx(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3131,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3132,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3133,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3134,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3135,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3136,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3137,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3138,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3139,SWMPI_COMM,s9,ierr)
!from Z2
call MPI_RECV_INIT(Txx(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3131,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3132,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3133,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3134,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3135,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3136,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3137,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3138,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3139,SWMPI_COMM,r9,ierr)
! put into array
reqZF(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

!to Z2
call MPI_SEND_INIT(Txx(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3211,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3212,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3213,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3214,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3215,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3219,SWMPI_COMM,s9,ierr)
!from Z1
call MPI_RECV_INIT(Txx(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3211,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3212,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3213,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3214,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3215,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3219,SWMPI_COMM,r9,ierr)
! put into array
reqZF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
!---------------------------------------------------------------------------
end subroutine mesg_init_LzF

subroutine mesg_init_LzB
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
! --- LzB ------------------------------------------------------------------
!to Z1
call MPI_SEND_INIT(Txx(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3111,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3112,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3113,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3114,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3115,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3116,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3117,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3118,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3119,SWMPI_COMM,s9,ierr)
!from Z2
call MPI_RECV_INIT(Txx(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3111,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3112,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3113,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3114,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3115,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3116,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3117,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3118,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3119,SWMPI_COMM,r9,ierr)
! put into array
reqZB(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

!to Z2
call MPI_SEND_INIT(Txx(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3231,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3232,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3233,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3234,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3235,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3239,SWMPI_COMM,s9,ierr)
!from Z2
call MPI_RECV_INIT(Txx(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3231,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3232,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3233,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3234,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3235,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3239,SWMPI_COMM,r9,ierr)
! put into array
reqZB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
end subroutine mesg_init_LzB

end module macdrp_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
