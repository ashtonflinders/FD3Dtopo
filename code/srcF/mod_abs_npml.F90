module abs_mod

! This module is used for absorbing outgoing waves based on
! unsplit-field ADE CFS-PML (auxiliary differential equation
! complex frequncy shifted PML)
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2008 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-14 16:30:51 -0500 (Wed, 14 Jan 2009) $
! $Revision: 510 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************
 
#include "mod_macdrp.h"
#define AbsVzero
!#define DEBUG
!#define CorrAbs

use constants_mod
use string_mod
use para_mod
use nfseis_mod
use macdrp_mod, only :                        &
    Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,         &
    hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,hVx,hVy,hVz,&
    matVx2Vz,matVy2Vz
use media_mod
use grid_mod
use mpi_mod
use mpi
implicit none

private
public ::            &
    abs_init,        &
    abs_destroy,     &
    abs_LxF_LyF_LzF, &
    abs_LxB_LyB_LzB, &
    abs_LxF_LyF_LzB, &
    abs_LxB_LyB_LzF, &
    abs_LxB_LyF_LzF, &
    abs_LxF_LyB_LzB, &
    abs_LxF_LyB_LzF, &
    abs_LxB_LyF_LzB, &
    abs_syn,         &
    abs_RK_beg,      &
    abs_RK_inn,      &
    abs_RK_fin

interface abs_LxF_LyF_LzF
  module procedure in_LxF_LyF_LzF
end interface
interface abs_LxB_LyB_LzB
  module procedure in_LxB_LyB_LzB
end interface
interface abs_LxF_LyF_LzB
  module procedure in_LxF_LyF_LzB
end interface
interface abs_LxB_LyB_LzF
  module procedure in_LxB_LyB_LzF
end interface
interface abs_LxB_LyF_LzF
  module procedure in_LxB_LyF_LzF
end interface
interface abs_LxF_LyB_LzB
  module procedure in_LxF_LyB_LzB
end interface
interface abs_LxF_LyB_LzF
  module procedure in_LxF_LyB_LzF
end interface
interface abs_LxB_LyF_LzB
  module procedure in_LxB_LyF_LzB
end interface

!-----------------------------------------------------------------------------
DEFFDWET 
DEFFDWET24
DEFFDWET22
DEFLDDRK4A
DEFLDDRK4B
HOCWETL
HOCWETR

type ELEM
     logical :: isabs,isx,isy,isz
     integer :: ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk
     real(SP),dimension(:,:,:),allocatable ::                &
          Txx1, Tyy1, Txy1, Vx1, Vy1, Tzz1, Txz1, Tyz1, Vz1, &
         mTxx1,mTyy1,mTxy1,mVx1,mVy1,mTzz1,mTxz1,mTyz1,mVz1, &
         tTxx1,tTyy1,tTxy1,tVx1,tVy1,tTzz1,tTxz1,tTyz1,tVz1, &
         hTxx1,hTyy1,hTxy1,hVx1,hVy1,hTzz1,hTxz1,hTyz1,hVz1
     real(SP),dimension(:,:,:),allocatable ::                &
          Txx2, Tyy2, Txy2, Vx2, Vy2, Tzz2, Txz2, Tyz2, Vz2, &
         mTxx2,mTyy2,mTxy2,mVx2,mVy2,mTzz2,mTxz2,mTyz2,mVz2, &
         tTxx2,tTyy2,tTxy2,tVx2,tVy2,tTzz2,tTxz2,tTyz2,tVz2, &
         hTxx2,hTyy2,hTxy2,hVx2,hVy2,hTzz2,hTxz2,hTyz2,hVz2
     real(SP),dimension(:,:,:),allocatable ::                &
          Txx3, Tyy3, Txy3, Vx3, Vy3, Tzz3, Txz3, Tyz3, Vz3, &
         mTxx3,mTyy3,mTxy3,mVx3,mVy3,mTzz3,mTxz3,mTyz3,mVz3, &
         tTxx3,tTyy3,tTxy3,tVx3,tVy3,tTzz3,tTxz3,tTyz3,tVz3, &
         hTxx3,hTyy3,hTxy3,hVx3,hVy3,hTzz3,hTxz3,hTyz3,hVz3
end type ELEM

integer :: num_blk=26

type (ELEM), dimension(26) :: W

integer,dimension(SEIS_GEO,2) :: abs_number
real(SP),dimension(:),allocatable :: Ax,Ay,Az,Bx,By,Bz,Dx,Dy,Dz
real(SP),dimension(:),allocatable :: Ex,Ey,Ez
!real(SP),dimension(SEIS_GEO,SEIS_GEO) :: MulD

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine abs_init(fnm_conf)
use mpi_mod, only : absnode
character (len=*),intent(in) :: fnm_conf
integer fid,n,m,i,j,k
real(SP),dimension(SEIS_GEO,2) :: Vs,fc,bmax
real(SP),dimension(SEIS_GEO) :: vec1
real(SP) :: L

fid=1001
abs_number=0
Vs=0.0
open(fid,file=trim(fnm_conf),status="old")
do n=1,SEIS_GEO
do m=1,2
if (absnode(n,m)) then
   call string_conf(fid,1,'abs_number',(n-1)*2+m+1,abs_number(n,m))
   call string_conf(fid,1,'abs_velocity',(n-1)*2+m+1,Vs(n,m))
   call string_conf(fid,1,'CFS_bmax',(n-1)*2+m+1,bmax(n,m))
   call string_conf(fid,1,'CFS_fc',(n-1)*2+m+1,fc(n,m))
   if (bmax(n,m)<1.0) then
      print *, "CFS_bmax should be large or equal to 1"
      print *, "n,m,CFS_bmax(n,m)=",n,m,bmax(n,m)
      stop 1
   end if
end if
end do
end do
close(fid)

num_blk=26

do n=1,num_blk
   W(n)%isabs=.false.
   W(n)%isx=.false.
   W(n)%isy=.false.
   W(n)%isz=.false.
end do

allocate(Ax(nx1:nx2)); Ax=0.0_SP
allocate(Bx(nx1:nx2)); Bx=1.0_SP
allocate(Dx(nx1:nx2)); Dx=0.0_SP
allocate(Ay(ny1:ny2)); Ay=0.0_SP
allocate(By(ny1:ny2)); By=1.0_SP
allocate(Dy(ny1:ny2)); Dy=0.0_SP
allocate(Az(nz1:nz2)); Az=0.0_SP
allocate(Bz(nz1:nz2)); Bz=1.0_SP
allocate(Dz(nz1:nz2)); Dz=0.0_SP
allocate(Ex(nx1:nx2)); Ex=1.0_SP
allocate(Ey(ny1:ny2)); Ey=1.0_SP
allocate(Ez(nz1:nz2)); Ez=1.0_SP

do i=ni1,abs_number(1,1)+ni1-1
   j=nj1; k=nk1
   call grid_covariant(i,j,k,1,vec1)
   L=sqrt(dot_product(vec1,vec1))*steph
   Dx(i)=cal_pml_d(abs_number(1,1)-(i-ni1),Vs(1,1),L,abs_number(1,1))
   Ax(i)=cal_pml_a(abs_number(1,1)-(i-ni1),fc(1,1),abs_number(1,1))
   Bx(i)=cal_pml_b(abs_number(1,1)-(i-ni1),bmax(1,1),abs_number(1,1))
   Ex(i)=cal_pml_e(abs_number(1,1)-(i-ni1),Vs(1,1),L,abs_number(1,1))
end do
do i=ni2-abs_number(1,2)+1,ni2
   j=nj1; k=nk1
   call grid_covariant(i,j,k,1,vec1)
   L=sqrt(dot_product(vec1,vec1))*steph
   Dx(i)=cal_pml_d(i-(ni2-abs_number(1,2)),Vs(1,2),L,abs_number(1,2))
   Ax(i)=cal_pml_a(i-(ni2-abs_number(1,2)),fc(1,2),abs_number(1,2))
   Bx(i)=cal_pml_b(i-(ni2-abs_number(1,2)),bmax(1,2),abs_number(1,2))
   Ex(i)=cal_pml_e(i-(ni2-abs_number(1,2)),Vs(1,2),L,abs_number(1,2))
end do
do j=nj1,abs_number(2,1)+nj1-1
   i=ni1; k=nk1
   call grid_covariant(i,j,k,2,vec1)
   L=sqrt(dot_product(vec1,vec1))*steph
   Dy(j)=cal_pml_d(abs_number(2,1)-(j-nj1),Vs(2,1),L,abs_number(2,1))
   Ay(j)=cal_pml_a(abs_number(2,1)-(j-nj1),fc(2,1),abs_number(2,1))
   By(j)=cal_pml_b(abs_number(2,1)-(j-nj1),bmax(2,1),abs_number(2,1))
   Ey(j)=cal_pml_e(abs_number(2,1)-(j-nj1),Vs(2,1),L,abs_number(2,1))
end do
do j=nj2-abs_number(2,2)+1,nj2
   i=ni1; k=nk1
   call grid_covariant(i,j,k,2,vec1)
   L=sqrt(dot_product(vec1,vec1))*steph
   Dy(j)=cal_pml_d(j-(nj2-abs_number(2,2)),Vs(2,2),L,abs_number(2,2))
   Ay(j)=cal_pml_a(j-(nj2-abs_number(2,2)),fc(2,2),abs_number(2,2))
   By(j)=cal_pml_b(j-(nj2-abs_number(2,2)),bmax(2,2),abs_number(2,2))
   Ey(j)=cal_pml_e(j-(nj2-abs_number(2,2)),Vs(2,2),L,abs_number(2,2))
end do
do k=nk1,nk1+abs_number(3,1)-1
   i=ni1; j=nj1
   call grid_covariant(i,j,k,3,vec1)
   L=sqrt(dot_product(vec1,vec1))*steph
   Dz(k)=cal_pml_d(abs_number(3,1)-(k-nk1),Vs(3,1),L,abs_number(3,1))
   Az(k)=cal_pml_a(abs_number(3,1)-(k-nk1),fc(3,1),abs_number(3,1))
   Bz(k)=cal_pml_b(abs_number(3,1)-(k-nk1),bmax(3,1),abs_number(3,1))
   Ez(k)=cal_pml_e(abs_number(3,1)-(k-nk1),Vs(3,1),L,abs_number(3,1))
end do
do k=nk2-abs_number(3,2)+1,nk2
   i=ni1; j=nj1
   call grid_covariant(i,j,k,3,vec1)
   L=sqrt(dot_product(vec1,vec1))*steph
   Dz(k)=cal_pml_d(k-(nk2-abs_number(3,2)),Vs(3,2),L,abs_number(3,2))
   Az(k)=cal_pml_a(k-(nk2-abs_number(3,2)),fc(3,2),abs_number(3,2))
   Bz(k)=cal_pml_b(k-(nk2-abs_number(3,2)),bmax(3,2),abs_number(3,2))
   Ez(k)=cal_pml_e(k-(nk2-abs_number(3,2)),Vs(3,2),L,abs_number(3,2))
end do

!x1y1z1
n=1
W(n)%ni=abs_number(1,1);W(n)%ni1=ni1; W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1);W(n)%nj1=nj1; W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1);W(n)%nk1=nk1; W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x2y1z1
n=n+1
W(n)%ni=abs_number(1,2);W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1);W(n)%nj1=nj1                  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1);W(n)%nk1=nk1                  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x1y2z1
n=n+1
W(n)%ni=abs_number(1,1);W(n)%ni1=ni1                  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2);W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1);W(n)%nk1=nk1                  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x2y2z1
n=n+1
W(n)%ni=abs_number(1,2);W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2);W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1);W(n)%nk1=nk1                  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x1y1z2
n=n+1
W(n)%ni=abs_number(1,1);W(n)%ni1=ni1                  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1);W(n)%nj1=nj1                  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2);W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x2y1z2
n=n+1
W(n)%ni=abs_number(1,2);W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1);W(n)%nj1=nj1                  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2);W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x1y2z2
n=n+1
W(n)%ni=abs_number(1,1);W(n)%ni1=ni1                  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2);W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2);W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x2y2z2
n=n+1
W(n)%ni=abs_number(1,2);W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2);W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2);W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if

!x1y1
n=n+1
W(n)%ni=abs_number(1,1)        ;W(n)%ni1=ni1                ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1)        ;W(n)%nj1=nj1                ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk-sum(abs_number(3,:));W(n)%nk1=nk1+abs_number(3,1);W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
else
   n=n-1
end if
!x2y1
n=n+1
W(n)%ni=abs_number(1,2)        ;W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1)        ;W(n)%nj1=nj1                  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk-sum(abs_number(3,:));W(n)%nk1=nk1+abs_number(3,1)  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
else
   n=n-1
end if
!x1y2
n=n+1
W(n)%ni=abs_number(1,1)        ;W(n)%ni1=ni1                  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2)        ;W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk-sum(abs_number(3,:));W(n)%nk1=nk1+abs_number(3,1)  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
else
   n=n-1
end if
!x2y2
n=n+1
W(n)%ni=abs_number(1,2)        ;W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2)        ;W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk-sum(abs_number(3,:));W(n)%nk1=nk1+abs_number(3,1)  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isy=.true.
else
   n=n-1
end if
!x1z1
n=n+1
W(n)%ni=abs_number(1,1)        ;W(n)%ni1=ni1                ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1);W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1)        ;W(n)%nk1=nk1                ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x2z1
n=n+1
W(n)%ni=abs_number(1,2)        ;W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1)  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1)        ;W(n)%nk1=nk1                  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x1z2
n=n+1
W(n)%ni=abs_number(1,1)        ;W(n)%ni1=ni1                  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1)  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2)        ;W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!x2z2
n=n+1
W(n)%ni=abs_number(1,2)        ;W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1)  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2)        ;W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!y1z1
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1);W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1)        ;W(n)%nj1=nj1                ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1)        ;W(n)%nk1=nk1                ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!y2z1
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1);W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2)        ;W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1)        ;W(n)%nk1=nk1                ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!y1z2
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1)  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1)        ;W(n)%nj1=nj1                  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2)        ;W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!y2z2
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1)  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2)        ;W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2)        ;W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isy=.true.
   W(n)%isz=.true.
else
   n=n-1
end if

!x1
n=n+1
W(n)%ni=abs_number(1,1)        ;W(n)%ni1=ni1                ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1);W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk-sum(abs_number(3,:));W(n)%nk1=nk1+abs_number(3,1);W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
else
   n=n-1
end if
!x2
n=n+1
W(n)%ni=abs_number(1,2)        ;W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1)  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk-sum(abs_number(3,:));W(n)%nk1=nk1+abs_number(3,1)  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isx=.true.
else
   n=n-1
end if
!y1
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1);W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,1)        ;W(n)%nj1=nj1                ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk-sum(abs_number(3,:));W(n)%nk1=nk1+abs_number(3,1);W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isy=.true.
else
   n=n-1
end if
!y2
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1)  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=abs_number(2,2)        ;W(n)%nj1=nj2-abs_number(2,2)+1;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=nk-sum(abs_number(3,:));W(n)%nk1=nk1+abs_number(3,1)  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isy=.true.
else
   n=n-1
end if
!z1
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1);W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1);W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,1)        ;W(n)%nk1=nk1                ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isz=.true.
else
   n=n-1
end if
!z2
n=n+1
W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1)  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
W(n)%nj=nj-sum(abs_number(2,:));W(n)%nj1=nj1+abs_number(2,1)  ;W(n)%nj2=W(n)%nj1+W(n)%nj-1; 
W(n)%nk=abs_number(3,2)        ;W(n)%nk1=nk2-abs_number(3,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
if (W(n)%ni>0 .and. W(n)%nj>0 .and. W(n)%nk>0) then
   W(n)%isabs=.true.
   W(n)%isz=.true.
else
   n=n-1
end if

num_blk=n
!if (n/=num_blk) then
!   print *, 'n/=num_blk'
!   stop 1
!end if

do n=1,num_blk
if (W(n)%isx) then
   allocate(W(n)%  Vx1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%  Vx1=0.0_SP
   allocate(W(n)%  Vy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%  Vy1=0.0_SP
   allocate(W(n)%  Vz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%  Vz1=0.0_SP
   allocate(W(n)% Txx1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Txx1=0.0_SP
   allocate(W(n)% Tyy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Tyy1=0.0_SP
   allocate(W(n)% Tzz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Tzz1=0.0_SP
   allocate(W(n)% Txy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Txy1=0.0_SP
   allocate(W(n)% Txz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Txz1=0.0_SP
   allocate(W(n)% Tyz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Tyz1=0.0_SP
   allocate(W(n)% hVx1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% hVx1=0.0_SP
   allocate(W(n)% hVy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% hVy1=0.0_SP
   allocate(W(n)% hVz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% hVz1=0.0_SP
   allocate(W(n)%hTxx1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTxx1=0.0_SP
   allocate(W(n)%hTyy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTyy1=0.0_SP
   allocate(W(n)%hTzz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTzz1=0.0_SP
   allocate(W(n)%hTxy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTxy1=0.0_SP
   allocate(W(n)%hTxz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTxz1=0.0_SP
   allocate(W(n)%hTyz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTyz1=0.0_SP
   allocate(W(n)% mVx1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% mVx1=0.0_SP
   allocate(W(n)% mVy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% mVy1=0.0_SP
   allocate(W(n)% mVz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% mVz1=0.0_SP
   allocate(W(n)%mTxx1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTxx1=0.0_SP
   allocate(W(n)%mTyy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTyy1=0.0_SP
   allocate(W(n)%mTzz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTzz1=0.0_SP
   allocate(W(n)%mTxy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTxy1=0.0_SP
   allocate(W(n)%mTxz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTxz1=0.0_SP
   allocate(W(n)%mTyz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTyz1=0.0_SP
   allocate(W(n)% tVx1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% tVx1=0.0_SP
   allocate(W(n)% tVy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% tVy1=0.0_SP
   allocate(W(n)% tVz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% tVz1=0.0_SP
   allocate(W(n)%tTxx1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTxx1=0.0_SP
   allocate(W(n)%tTyy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTyy1=0.0_SP
   allocate(W(n)%tTzz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTzz1=0.0_SP
   allocate(W(n)%tTxy1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTxy1=0.0_SP
   allocate(W(n)%tTxz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTxz1=0.0_SP
   allocate(W(n)%tTyz1(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTyz1=0.0_SP
end if
if (W(n)%isy) then
   allocate(W(n)%  Vx2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%  Vx2=0.0_SP
   allocate(W(n)%  Vy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%  Vy2=0.0_SP
   allocate(W(n)%  Vz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%  Vz2=0.0_SP
   allocate(W(n)% Txx2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Txx2=0.0_SP
   allocate(W(n)% Tyy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Tyy2=0.0_SP
   allocate(W(n)% Tzz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Tzz2=0.0_SP
   allocate(W(n)% Txy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Txy2=0.0_SP
   allocate(W(n)% Txz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Txz2=0.0_SP
   allocate(W(n)% Tyz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Tyz2=0.0_SP
   allocate(W(n)% hVx2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% hVx2=0.0_SP
   allocate(W(n)% hVy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% hVy2=0.0_SP
   allocate(W(n)% hVz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% hVz2=0.0_SP
   allocate(W(n)%hTxx2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTxx2=0.0_SP
   allocate(W(n)%hTyy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTyy2=0.0_SP
   allocate(W(n)%hTzz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTzz2=0.0_SP
   allocate(W(n)%hTxy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTxy2=0.0_SP
   allocate(W(n)%hTxz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTxz2=0.0_SP
   allocate(W(n)%hTyz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTyz2=0.0_SP
   allocate(W(n)% mVx2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% mVx2=0.0_SP
   allocate(W(n)% mVy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% mVy2=0.0_SP
   allocate(W(n)% mVz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% mVz2=0.0_SP
   allocate(W(n)%mTxx2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTxx2=0.0_SP
   allocate(W(n)%mTyy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTyy2=0.0_SP
   allocate(W(n)%mTzz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTzz2=0.0_SP
   allocate(W(n)%mTxy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTxy2=0.0_SP
   allocate(W(n)%mTxz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTxz2=0.0_SP
   allocate(W(n)%mTyz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTyz2=0.0_SP
   allocate(W(n)% tVx2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% tVx2=0.0_SP
   allocate(W(n)% tVy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% tVy2=0.0_SP
   allocate(W(n)% tVz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% tVz2=0.0_SP
   allocate(W(n)%tTxx2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTxx2=0.0_SP
   allocate(W(n)%tTyy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTyy2=0.0_SP
   allocate(W(n)%tTzz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTzz2=0.0_SP
   allocate(W(n)%tTxy2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTxy2=0.0_SP
   allocate(W(n)%tTxz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTxz2=0.0_SP
   allocate(W(n)%tTyz2(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTyz2=0.0_SP
end if
if (W(n)%isz) then
   allocate(W(n)%  Vx3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%  Vx3=0.0_SP
   allocate(W(n)%  Vy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%  Vy3=0.0_SP
   allocate(W(n)%  Vz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%  Vz3=0.0_SP
   allocate(W(n)% Txx3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Txx3=0.0_SP
   allocate(W(n)% Tyy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Tyy3=0.0_SP
   allocate(W(n)% Tzz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Tzz3=0.0_SP
   allocate(W(n)% Txy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Txy3=0.0_SP
   allocate(W(n)% Txz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Txz3=0.0_SP
   allocate(W(n)% Tyz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% Tyz3=0.0_SP
   allocate(W(n)% hVx3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% hVx3=0.0_SP
   allocate(W(n)% hVy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% hVy3=0.0_SP
   allocate(W(n)% hVz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% hVz3=0.0_SP
   allocate(W(n)%hTxx3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTxx3=0.0_SP
   allocate(W(n)%hTyy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTyy3=0.0_SP
   allocate(W(n)%hTzz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTzz3=0.0_SP
   allocate(W(n)%hTxy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTxy3=0.0_SP
   allocate(W(n)%hTxz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTxz3=0.0_SP
   allocate(W(n)%hTyz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%hTyz3=0.0_SP
   allocate(W(n)% mVx3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% mVx3=0.0_SP
   allocate(W(n)% mVy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% mVy3=0.0_SP
   allocate(W(n)% mVz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% mVz3=0.0_SP
   allocate(W(n)%mTxx3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTxx3=0.0_SP
   allocate(W(n)%mTyy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTyy3=0.0_SP
   allocate(W(n)%mTzz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTzz3=0.0_SP
   allocate(W(n)%mTxy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTxy3=0.0_SP
   allocate(W(n)%mTxz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTxz3=0.0_SP
   allocate(W(n)%mTyz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%mTyz3=0.0_SP
   allocate(W(n)% tVx3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% tVx3=0.0_SP
   allocate(W(n)% tVy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% tVy3=0.0_SP
   allocate(W(n)% tVz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)% tVz3=0.0_SP
   allocate(W(n)%tTxx3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTxx3=0.0_SP
   allocate(W(n)%tTyy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTyy3=0.0_SP
   allocate(W(n)%tTzz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTzz3=0.0_SP
   allocate(W(n)%tTxy3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTxy3=0.0_SP
   allocate(W(n)%tTxz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTxz3=0.0_SP
   allocate(W(n)%tTyz3(W(n)%ni1:W(n)%ni2,W(n)%nj1:W(n)%nj2,W(n)%nk1:W(n)%nk2));W(n)%tTyz3=0.0_SP
end if
end do

if (freenode) then
do n=1,num_blk

if (W(n)%isx) then
do k=W(n)%nk1,W(n)%nk2
if (k==nk2) then
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   call coef_fdxy2fdz_abs(i,j)
end do
end do
end if
end do
end if

if (W(n)%isy) then
do k=W(n)%nk1,W(n)%nk2
if (k==nk2) then
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   call coef_fdxy2fdz_abs(i,j)
end do
end do
end if
end do
end if

if (W(n)%isz) then
do k=W(n)%nk1,W(n)%nk2
if (k==nk2) then
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   call coef_fdxy2fdz_abs(i,j)
end do
end do
end if
end do
end if

end do
end if

#ifdef DEBUG
do i=nx1,nx2
write(200+myid,"(3es14.5)") Ax(i),Bx(i),Dx(i)
end do
do j=ny1,ny2
write(200+myid,"(3es14.5)") Ay(j),By(j),Dy(j)
end do
do k=nz1,nz2
write(200+myid,"(3es14.5)") Az(k),Bz(k),Dz(k)
end do
do n=1,num_blk
write(200+myid,*) "n=",n
write(200+myid,*) W(n)%isabs,W(n)%isx,W(n)%isy,W(n)%isz
write(200+myid,"(3i6)") W(n)%ni1,W(n)%ni2,W(n)%ni
write(200+myid,"(3i6)") W(n)%nj1,W(n)%nj2,W(n)%nj
write(200+myid,"(3i6)") W(n)%nk1,W(n)%nk2,W(n)%nk
end do
#endif
#ifdef AbsExport
call absnc_create
write(300+myid,*) W(3)%nk1,W(3)%nk2,W(3)%nk
i_nc=30
nt_nc=0
#endif
#ifndef CorrAbs
call abs_indx
#endif
end subroutine abs_init

subroutine abs_indx
use macdrp_mod, only : indx
if (absnode(1,1)) then
   indx(1,SEIS_GEO*2+1)=ni1+abs_number(1,1)
   indx(1,SEIS_GEO*2  )=ni1+abs_number(1,1)
   indx(1,SEIS_GEO*2-1)=ni1+abs_number(1,1)
   indx(1,           1)=ni1+abs_number(1,1)
   indx(1,           2)=ni1+abs_number(1,1)
   indx(1,           3)=ni1+abs_number(1,1)
end if
if (absnode(1,2)) then
   indx(2,SEIS_GEO*2+1)=ni2-abs_number(1,2)
   indx(2,SEIS_GEO*2  )=ni2-abs_number(1,2)
   indx(2,SEIS_GEO*2-1)=ni2-abs_number(1,2)
   indx(2,           1)=ni2-abs_number(1,2)
   indx(2,           2)=ni2-abs_number(1,2)
   indx(2,           4)=ni2-abs_number(1,2)
end if
if (absnode(2,1)) then
   indx(3,SEIS_GEO*2+1)=nj1+abs_number(2,1)
   indx(3,SEIS_GEO*2  )=nj1+abs_number(2,1)
   indx(3,SEIS_GEO*2-1)=nj1+abs_number(2,1)
   indx(3,           1)=nj1+abs_number(2,1)
   indx(3,           3)=nj1+abs_number(2,1)
   indx(3,           4)=nj1+abs_number(2,1)
end if
if (absnode(2,2)) then
   indx(4,SEIS_GEO*2+1)=nj2-abs_number(2,2)
   indx(4,SEIS_GEO*2  )=nj2-abs_number(2,2)
   indx(4,SEIS_GEO*2-1)=nj2-abs_number(2,2)
   indx(4,           2)=nj2-abs_number(2,2)
   indx(4,           3)=nj2-abs_number(2,2)
   indx(4,           4)=nj2-abs_number(2,2)
end if
if (absnode(3,1)) then
   indx(5,SEIS_GEO*2+1)=nk1+abs_number(3,1)
   indx(5,SEIS_GEO*2-1)=nk1+abs_number(3,1)
   indx(5,           1)=nk1+abs_number(3,1)
   indx(5,           2)=nk1+abs_number(3,1)
   indx(5,           3)=nk1+abs_number(3,1)
   indx(5,           4)=nk1+abs_number(3,1)
end if
if (absnode(3,2)) then
   indx(6,SEIS_GEO*2+1)=nk2-abs_number(3,2)
   indx(6,SEIS_GEO*2  )=nk2-abs_number(3,2)
   indx(6,           1)=nk2-abs_number(3,2)
   indx(6,           2)=nk2-abs_number(3,2)
   indx(6,           3)=nk2-abs_number(3,2)
   indx(6,           4)=nk2-abs_number(3,2)
end if
end subroutine abs_indx

subroutine abs_destroy
integer n
do n=1,num_blk
if (allocated(W(n)%Txx1)) deallocate(W(n)%Txx1)
end do
end subroutine abs_destroy

!-----------------------------------------------------------------------------

subroutine abs_syn
 integer n,i,j,k
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(n,i,j,k)
do n=1,num_blk
if (W(n)%isx) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)% mVx1(i,j,k)=W(n)% Vx1(i,j,k)
   W(n)% mVy1(i,j,k)=W(n)% Vy1(i,j,k)
   W(n)% mVz1(i,j,k)=W(n)% Vz1(i,j,k)
   W(n)%mTxx1(i,j,k)=W(n)%Txx1(i,j,k)
   W(n)%mTyy1(i,j,k)=W(n)%Tyy1(i,j,k)
   W(n)%mTzz1(i,j,k)=W(n)%Tzz1(i,j,k)
   W(n)%mTxy1(i,j,k)=W(n)%Txy1(i,j,k)
   W(n)%mTxz1(i,j,k)=W(n)%Txz1(i,j,k)
   W(n)%mTyz1(i,j,k)=W(n)%Tyz1(i,j,k)
end do
end do
end do
end if
if (W(n)%isy) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)% mVx2(i,j,k)=W(n)% Vx2(i,j,k)
   W(n)% mVy2(i,j,k)=W(n)% Vy2(i,j,k)
   W(n)% mVz2(i,j,k)=W(n)% Vz2(i,j,k)
   W(n)%mTxx2(i,j,k)=W(n)%Txx2(i,j,k)
   W(n)%mTyy2(i,j,k)=W(n)%Tyy2(i,j,k)
   W(n)%mTzz2(i,j,k)=W(n)%Tzz2(i,j,k)
   W(n)%mTxy2(i,j,k)=W(n)%Txy2(i,j,k)
   W(n)%mTxz2(i,j,k)=W(n)%Txz2(i,j,k)
   W(n)%mTyz2(i,j,k)=W(n)%Tyz2(i,j,k)
end do
end do
end do
end if
if (W(n)%isz) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)% mVx3(i,j,k)=W(n)% Vx3(i,j,k)
   W(n)% mVy3(i,j,k)=W(n)% Vy3(i,j,k)
   W(n)% mVz3(i,j,k)=W(n)% Vz3(i,j,k)
   W(n)%mTxx3(i,j,k)=W(n)%Txx3(i,j,k)
   W(n)%mTyy3(i,j,k)=W(n)%Tyy3(i,j,k)
   W(n)%mTzz3(i,j,k)=W(n)%Tzz3(i,j,k)
   W(n)%mTxy3(i,j,k)=W(n)%Txy3(i,j,k)
   W(n)%mTxz3(i,j,k)=W(n)%Txz3(i,j,k)
   W(n)%mTyz3(i,j,k)=W(n)%Tyz3(i,j,k)
end do
end do
end do
end if
end do
!$OMP END PARALLEL DO
end subroutine abs_syn

subroutine abs_RK_beg(rka,rkb)
 real(SP),intent(in) :: rka,rkb
 real(SP) a,b
 real(SP) e
 integer n,i,j,k
 a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
do n=1,num_blk
if (W(n)%isx) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)%  Vx1(i,j,k)=W(n)% mVx1(i,j,k)+a*W(n)% hVx1(i,j,k)
   W(n)%  Vy1(i,j,k)=W(n)% mVy1(i,j,k)+a*W(n)% hVy1(i,j,k)
   W(n)%  Vz1(i,j,k)=W(n)% mVz1(i,j,k)+a*W(n)% hVz1(i,j,k)
   W(n)% Txx1(i,j,k)=W(n)%mTxx1(i,j,k)+a*W(n)%hTxx1(i,j,k)
   W(n)% Tyy1(i,j,k)=W(n)%mTyy1(i,j,k)+a*W(n)%hTyy1(i,j,k)
   W(n)% Tzz1(i,j,k)=W(n)%mTzz1(i,j,k)+a*W(n)%hTzz1(i,j,k)
   W(n)% Txy1(i,j,k)=W(n)%mTxy1(i,j,k)+a*W(n)%hTxy1(i,j,k)
   W(n)% Txz1(i,j,k)=W(n)%mTxz1(i,j,k)+a*W(n)%hTxz1(i,j,k)
   W(n)% Tyz1(i,j,k)=W(n)%mTyz1(i,j,k)+a*W(n)%hTyz1(i,j,k)
   W(n)% tVx1(i,j,k)=W(n)% mVx1(i,j,k)+b*W(n)% hVx1(i,j,k)
   W(n)% tVy1(i,j,k)=W(n)% mVy1(i,j,k)+b*W(n)% hVy1(i,j,k)
   W(n)% tVz1(i,j,k)=W(n)% mVz1(i,j,k)+b*W(n)% hVz1(i,j,k)
   W(n)%tTxx1(i,j,k)=W(n)%mTxx1(i,j,k)+b*W(n)%hTxx1(i,j,k)
   W(n)%tTyy1(i,j,k)=W(n)%mTyy1(i,j,k)+b*W(n)%hTyy1(i,j,k)
   W(n)%tTzz1(i,j,k)=W(n)%mTzz1(i,j,k)+b*W(n)%hTzz1(i,j,k)
   W(n)%tTxy1(i,j,k)=W(n)%mTxy1(i,j,k)+b*W(n)%hTxy1(i,j,k)
   W(n)%tTxz1(i,j,k)=W(n)%mTxz1(i,j,k)+b*W(n)%hTxz1(i,j,k)
   W(n)%tTyz1(i,j,k)=W(n)%mTyz1(i,j,k)+b*W(n)%hTyz1(i,j,k)
end do
end do
end do
end if
if (W(n)%isy) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)%  Vx2(i,j,k)=W(n)% mVx2(i,j,k)+a*W(n)% hVx2(i,j,k)
   W(n)%  Vy2(i,j,k)=W(n)% mVy2(i,j,k)+a*W(n)% hVy2(i,j,k)
   W(n)%  Vz2(i,j,k)=W(n)% mVz2(i,j,k)+a*W(n)% hVz2(i,j,k)
   W(n)% Txx2(i,j,k)=W(n)%mTxx2(i,j,k)+a*W(n)%hTxx2(i,j,k)
   W(n)% Tyy2(i,j,k)=W(n)%mTyy2(i,j,k)+a*W(n)%hTyy2(i,j,k)
   W(n)% Tzz2(i,j,k)=W(n)%mTzz2(i,j,k)+a*W(n)%hTzz2(i,j,k)
   W(n)% Txy2(i,j,k)=W(n)%mTxy2(i,j,k)+a*W(n)%hTxy2(i,j,k)
   W(n)% Txz2(i,j,k)=W(n)%mTxz2(i,j,k)+a*W(n)%hTxz2(i,j,k)
   W(n)% Tyz2(i,j,k)=W(n)%mTyz2(i,j,k)+a*W(n)%hTyz2(i,j,k)
   W(n)% tVx2(i,j,k)=W(n)% mVx2(i,j,k)+b*W(n)% hVx2(i,j,k)
   W(n)% tVy2(i,j,k)=W(n)% mVy2(i,j,k)+b*W(n)% hVy2(i,j,k)
   W(n)% tVz2(i,j,k)=W(n)% mVz2(i,j,k)+b*W(n)% hVz2(i,j,k)
   W(n)%tTxx2(i,j,k)=W(n)%mTxx2(i,j,k)+b*W(n)%hTxx2(i,j,k)
   W(n)%tTyy2(i,j,k)=W(n)%mTyy2(i,j,k)+b*W(n)%hTyy2(i,j,k)
   W(n)%tTzz2(i,j,k)=W(n)%mTzz2(i,j,k)+b*W(n)%hTzz2(i,j,k)
   W(n)%tTxy2(i,j,k)=W(n)%mTxy2(i,j,k)+b*W(n)%hTxy2(i,j,k)
   W(n)%tTxz2(i,j,k)=W(n)%mTxz2(i,j,k)+b*W(n)%hTxz2(i,j,k)
   W(n)%tTyz2(i,j,k)=W(n)%mTyz2(i,j,k)+b*W(n)%hTyz2(i,j,k)
end do
end do
end do
end if
if (W(n)%isz) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)%  Vx3(i,j,k)=W(n)% mVx3(i,j,k)+a*W(n)% hVx3(i,j,k)
   W(n)%  Vy3(i,j,k)=W(n)% mVy3(i,j,k)+a*W(n)% hVy3(i,j,k)
   W(n)%  Vz3(i,j,k)=W(n)% mVz3(i,j,k)+a*W(n)% hVz3(i,j,k)
   W(n)% Txx3(i,j,k)=W(n)%mTxx3(i,j,k)+a*W(n)%hTxx3(i,j,k)
   W(n)% Tyy3(i,j,k)=W(n)%mTyy3(i,j,k)+a*W(n)%hTyy3(i,j,k)
   W(n)% Tzz3(i,j,k)=W(n)%mTzz3(i,j,k)+a*W(n)%hTzz3(i,j,k)
   W(n)% Txy3(i,j,k)=W(n)%mTxy3(i,j,k)+a*W(n)%hTxy3(i,j,k)
   W(n)% Txz3(i,j,k)=W(n)%mTxz3(i,j,k)+a*W(n)%hTxz3(i,j,k)
   W(n)% Tyz3(i,j,k)=W(n)%mTyz3(i,j,k)+a*W(n)%hTyz3(i,j,k)
   W(n)% tVx3(i,j,k)=W(n)% mVx3(i,j,k)+b*W(n)% hVx3(i,j,k)
   W(n)% tVy3(i,j,k)=W(n)% mVy3(i,j,k)+b*W(n)% hVy3(i,j,k)
   W(n)% tVz3(i,j,k)=W(n)% mVz3(i,j,k)+b*W(n)% hVz3(i,j,k)
   W(n)%tTxx3(i,j,k)=W(n)%mTxx3(i,j,k)+b*W(n)%hTxx3(i,j,k)
   W(n)%tTyy3(i,j,k)=W(n)%mTyy3(i,j,k)+b*W(n)%hTyy3(i,j,k)
   W(n)%tTzz3(i,j,k)=W(n)%mTzz3(i,j,k)+b*W(n)%hTzz3(i,j,k)
   W(n)%tTxy3(i,j,k)=W(n)%mTxy3(i,j,k)+b*W(n)%hTxy3(i,j,k)
   W(n)%tTxz3(i,j,k)=W(n)%mTxz3(i,j,k)+b*W(n)%hTxz3(i,j,k)
   W(n)%tTyz3(i,j,k)=W(n)%mTyz3(i,j,k)+b*W(n)%hTyz3(i,j,k)
end do
end do
end do
end if
end do
!$OMP END PARALLEL DO
!if (freenode) then
!do n=1,num_blk
!if ( W(n)%nk2==nk2) then
!   k=W(n)%nk2
!
!if ( W(n)%isx ) then
!do j=W(n)%nj1,W(n)%nj2
!do i=W(n)%ni1,W(n)%ni2
!   W(n)%  Vx1(i,j,k)=0.0_SP
!   W(n)%  Vy1(i,j,k)=0.0_SP
!   W(n)%  Vz1(i,j,k)=0.0_SP
!   W(n)% Txx1(i,j,k)=0.0_SP
!   W(n)% Tyy1(i,j,k)=0.0_SP
!   W(n)% Tzz1(i,j,k)=0.0_SP
!   W(n)% Txy1(i,j,k)=0.0_SP
!   W(n)% Txz1(i,j,k)=0.0_SP
!   W(n)% Tyz1(i,j,k)=0.0_SP
!end do
!end do
!end if !isx
!if ( W(n)%isy ) then
!do j=W(n)%nj1,W(n)%nj2
!do i=W(n)%ni1,W(n)%ni2
!   W(n)%  Vx2(i,j,k)=0.0_SP
!   W(n)%  Vy2(i,j,k)=0.0_SP
!   W(n)%  Vz2(i,j,k)=0.0_SP
!   W(n)% Txx2(i,j,k)=0.0_SP
!   W(n)% Tyy2(i,j,k)=0.0_SP
!   W(n)% Tzz2(i,j,k)=0.0_SP
!   W(n)% Txy2(i,j,k)=0.0_SP
!   W(n)% Txz2(i,j,k)=0.0_SP
!   W(n)% Tyz2(i,j,k)=0.0_SP
!end do
!end do
!end if !isy
!
!if ( (W(n)%isx .or. W(n)%isy) ) then
!do j=W(n)%nj1,W(n)%nj2
!do i=W(n)%ni1,W(n)%ni2
!   e=min(Ex(i),Ey(j),Ez(k))
!    Vx (i,j,k)= Vx (i,j,k)*e
!    Vy (i,j,k)= Vy (i,j,k)*e
!    Vz (i,j,k)= Vz (i,j,k)*e
!   Txx (i,j,k)=Txx (i,j,k)*e
!   Tyy (i,j,k)=Tyy (i,j,k)*e
!   Tzz (i,j,k)=Tzz (i,j,k)*e
!   Txy (i,j,k)=Txy (i,j,k)*e
!   Txz (i,j,k)=Txz (i,j,k)*e
!   Tyz (i,j,k)=Tyz (i,j,k)*e
!end do
!end do
!end if !damp
!
!end if !nk2
!end do !nblk
!end if !freenode
#ifdef CondFreeCharac
   call swmpi_except('PML conflicts with charac boundary')
#endif
end subroutine abs_RK_beg

subroutine abs_RK_inn(rka,rkb)
 real(SP),intent(in) :: rka,rkb
 real(SP) a,b
 real(SP) e
 integer n,i,j,k
 a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,k)
do n=1,num_blk
if (W(n)%isx) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)%  Vx1(i,j,k)=W(n)% mVx1(i,j,k)+a*W(n)% hVx1(i,j,k)
   W(n)%  Vy1(i,j,k)=W(n)% mVy1(i,j,k)+a*W(n)% hVy1(i,j,k)
   W(n)%  Vz1(i,j,k)=W(n)% mVz1(i,j,k)+a*W(n)% hVz1(i,j,k)
   W(n)% Txx1(i,j,k)=W(n)%mTxx1(i,j,k)+a*W(n)%hTxx1(i,j,k)
   W(n)% Tyy1(i,j,k)=W(n)%mTyy1(i,j,k)+a*W(n)%hTyy1(i,j,k)
   W(n)% Tzz1(i,j,k)=W(n)%mTzz1(i,j,k)+a*W(n)%hTzz1(i,j,k)
   W(n)% Txy1(i,j,k)=W(n)%mTxy1(i,j,k)+a*W(n)%hTxy1(i,j,k)
   W(n)% Txz1(i,j,k)=W(n)%mTxz1(i,j,k)+a*W(n)%hTxz1(i,j,k)
   W(n)% Tyz1(i,j,k)=W(n)%mTyz1(i,j,k)+a*W(n)%hTyz1(i,j,k)
   W(n)% tVx1(i,j,k)=W(n)% tVx1(i,j,k)+b*W(n)% hVx1(i,j,k)
   W(n)% tVy1(i,j,k)=W(n)% tVy1(i,j,k)+b*W(n)% hVy1(i,j,k)
   W(n)% tVz1(i,j,k)=W(n)% tVz1(i,j,k)+b*W(n)% hVz1(i,j,k)
   W(n)%tTxx1(i,j,k)=W(n)%tTxx1(i,j,k)+b*W(n)%hTxx1(i,j,k)
   W(n)%tTyy1(i,j,k)=W(n)%tTyy1(i,j,k)+b*W(n)%hTyy1(i,j,k)
   W(n)%tTzz1(i,j,k)=W(n)%tTzz1(i,j,k)+b*W(n)%hTzz1(i,j,k)
   W(n)%tTxy1(i,j,k)=W(n)%tTxy1(i,j,k)+b*W(n)%hTxy1(i,j,k)
   W(n)%tTxz1(i,j,k)=W(n)%tTxz1(i,j,k)+b*W(n)%hTxz1(i,j,k)
   W(n)%tTyz1(i,j,k)=W(n)%tTyz1(i,j,k)+b*W(n)%hTyz1(i,j,k)
end do
end do
end do
end if
if (W(n)%isy) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)%  Vx2(i,j,k)=W(n)% mVx2(i,j,k)+a*W(n)% hVx2(i,j,k)
   W(n)%  Vy2(i,j,k)=W(n)% mVy2(i,j,k)+a*W(n)% hVy2(i,j,k)
   W(n)%  Vz2(i,j,k)=W(n)% mVz2(i,j,k)+a*W(n)% hVz2(i,j,k)
   W(n)% Txx2(i,j,k)=W(n)%mTxx2(i,j,k)+a*W(n)%hTxx2(i,j,k)
   W(n)% Tyy2(i,j,k)=W(n)%mTyy2(i,j,k)+a*W(n)%hTyy2(i,j,k)
   W(n)% Tzz2(i,j,k)=W(n)%mTzz2(i,j,k)+a*W(n)%hTzz2(i,j,k)
   W(n)% Txy2(i,j,k)=W(n)%mTxy2(i,j,k)+a*W(n)%hTxy2(i,j,k)
   W(n)% Txz2(i,j,k)=W(n)%mTxz2(i,j,k)+a*W(n)%hTxz2(i,j,k)
   W(n)% Tyz2(i,j,k)=W(n)%mTyz2(i,j,k)+a*W(n)%hTyz2(i,j,k)
   W(n)% tVx2(i,j,k)=W(n)% tVx2(i,j,k)+b*W(n)% hVx2(i,j,k)
   W(n)% tVy2(i,j,k)=W(n)% tVy2(i,j,k)+b*W(n)% hVy2(i,j,k)
   W(n)% tVz2(i,j,k)=W(n)% tVz2(i,j,k)+b*W(n)% hVz2(i,j,k)
   W(n)%tTxx2(i,j,k)=W(n)%tTxx2(i,j,k)+b*W(n)%hTxx2(i,j,k)
   W(n)%tTyy2(i,j,k)=W(n)%tTyy2(i,j,k)+b*W(n)%hTyy2(i,j,k)
   W(n)%tTzz2(i,j,k)=W(n)%tTzz2(i,j,k)+b*W(n)%hTzz2(i,j,k)
   W(n)%tTxy2(i,j,k)=W(n)%tTxy2(i,j,k)+b*W(n)%hTxy2(i,j,k)
   W(n)%tTxz2(i,j,k)=W(n)%tTxz2(i,j,k)+b*W(n)%hTxz2(i,j,k)
   W(n)%tTyz2(i,j,k)=W(n)%tTyz2(i,j,k)+b*W(n)%hTyz2(i,j,k)
end do
end do
end do
end if
if (W(n)%isz) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)%  Vx3(i,j,k)=W(n)% mVx3(i,j,k)+a*W(n)% hVx3(i,j,k)
   W(n)%  Vy3(i,j,k)=W(n)% mVy3(i,j,k)+a*W(n)% hVy3(i,j,k)
   W(n)%  Vz3(i,j,k)=W(n)% mVz3(i,j,k)+a*W(n)% hVz3(i,j,k)
   W(n)% Txx3(i,j,k)=W(n)%mTxx3(i,j,k)+a*W(n)%hTxx3(i,j,k)
   W(n)% Tyy3(i,j,k)=W(n)%mTyy3(i,j,k)+a*W(n)%hTyy3(i,j,k)
   W(n)% Tzz3(i,j,k)=W(n)%mTzz3(i,j,k)+a*W(n)%hTzz3(i,j,k)
   W(n)% Txy3(i,j,k)=W(n)%mTxy3(i,j,k)+a*W(n)%hTxy3(i,j,k)
   W(n)% Txz3(i,j,k)=W(n)%mTxz3(i,j,k)+a*W(n)%hTxz3(i,j,k)
   W(n)% Tyz3(i,j,k)=W(n)%mTyz3(i,j,k)+a*W(n)%hTyz3(i,j,k)
   W(n)% tVx3(i,j,k)=W(n)% tVx3(i,j,k)+b*W(n)% hVx3(i,j,k)
   W(n)% tVy3(i,j,k)=W(n)% tVy3(i,j,k)+b*W(n)% hVy3(i,j,k)
   W(n)% tVz3(i,j,k)=W(n)% tVz3(i,j,k)+b*W(n)% hVz3(i,j,k)
   W(n)%tTxx3(i,j,k)=W(n)%tTxx3(i,j,k)+b*W(n)%hTxx3(i,j,k)
   W(n)%tTyy3(i,j,k)=W(n)%tTyy3(i,j,k)+b*W(n)%hTyy3(i,j,k)
   W(n)%tTzz3(i,j,k)=W(n)%tTzz3(i,j,k)+b*W(n)%hTzz3(i,j,k)
   W(n)%tTxy3(i,j,k)=W(n)%tTxy3(i,j,k)+b*W(n)%hTxy3(i,j,k)
   W(n)%tTxz3(i,j,k)=W(n)%tTxz3(i,j,k)+b*W(n)%hTxz3(i,j,k)
   W(n)%tTyz3(i,j,k)=W(n)%tTyz3(i,j,k)+b*W(n)%hTyz3(i,j,k)
end do
end do
end do
end if
end do
!$OMP END PARALLEL DO
!if (freenode) then
!do n=1,num_blk
!if ( W(n)%nk2==nk2) then
!   k=W(n)%nk2
!
!if ( W(n)%isx ) then
!do j=W(n)%nj1,W(n)%nj2
!do i=W(n)%ni1,W(n)%ni2
!   W(n)%  Vx1(i,j,k)=0.0_SP
!   W(n)%  Vy1(i,j,k)=0.0_SP
!   W(n)%  Vz1(i,j,k)=0.0_SP
!   W(n)% Txx1(i,j,k)=0.0_SP
!   W(n)% Tyy1(i,j,k)=0.0_SP
!   W(n)% Tzz1(i,j,k)=0.0_SP
!   W(n)% Txy1(i,j,k)=0.0_SP
!   W(n)% Txz1(i,j,k)=0.0_SP
!   W(n)% Tyz1(i,j,k)=0.0_SP
!end do
!end do
!end if !isx
!if ( W(n)%isy ) then
!do j=W(n)%nj1,W(n)%nj2
!do i=W(n)%ni1,W(n)%ni2
!   W(n)%  Vx2(i,j,k)=0.0_SP
!   W(n)%  Vy2(i,j,k)=0.0_SP
!   W(n)%  Vz2(i,j,k)=0.0_SP
!   W(n)% Txx2(i,j,k)=0.0_SP
!   W(n)% Tyy2(i,j,k)=0.0_SP
!   W(n)% Tzz2(i,j,k)=0.0_SP
!   W(n)% Txy2(i,j,k)=0.0_SP
!   W(n)% Txz2(i,j,k)=0.0_SP
!   W(n)% Tyz2(i,j,k)=0.0_SP
!end do
!end do
!end if !isy
!
!if ( (W(n)%isx .or. W(n)%isy) ) then
!do j=W(n)%nj1,W(n)%nj2
!do i=W(n)%ni1,W(n)%ni2
!   e=min(Ex(i),Ey(j),Ez(k))
!    Vx (i,j,k)= Vx (i,j,k)*e
!    Vy (i,j,k)= Vy (i,j,k)*e
!    Vz (i,j,k)= Vz (i,j,k)*e
!   Txx (i,j,k)=Txx (i,j,k)*e
!   Tyy (i,j,k)=Tyy (i,j,k)*e
!   Tzz (i,j,k)=Tzz (i,j,k)*e
!   Txy (i,j,k)=Txy (i,j,k)*e
!   Txz (i,j,k)=Txz (i,j,k)*e
!   Tyz (i,j,k)=Tyz (i,j,k)*e
!end do
!end do
!end if !damp
!
!end if !nk2
!end do !nblk
!end if !freenode
end subroutine abs_RK_inn

subroutine abs_RK_fin(rkb)
 real(SP),intent(in) :: rkb
 real(SP) b
 real(SP) e
 integer n,i,j,k
 b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
do n=1,num_blk
if (W(n)%isx) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)% Vx1(i,j,k)=W(n)% tVx1(i,j,k)+b*W(n)% hVx1(i,j,k)
   W(n)% Vy1(i,j,k)=W(n)% tVy1(i,j,k)+b*W(n)% hVy1(i,j,k)
   W(n)% Vz1(i,j,k)=W(n)% tVz1(i,j,k)+b*W(n)% hVz1(i,j,k)
   W(n)%Txx1(i,j,k)=W(n)%tTxx1(i,j,k)+b*W(n)%hTxx1(i,j,k)
   W(n)%Tyy1(i,j,k)=W(n)%tTyy1(i,j,k)+b*W(n)%hTyy1(i,j,k)
   W(n)%Tzz1(i,j,k)=W(n)%tTzz1(i,j,k)+b*W(n)%hTzz1(i,j,k)
   W(n)%Txy1(i,j,k)=W(n)%tTxy1(i,j,k)+b*W(n)%hTxy1(i,j,k)
   W(n)%Txz1(i,j,k)=W(n)%tTxz1(i,j,k)+b*W(n)%hTxz1(i,j,k)
   W(n)%Tyz1(i,j,k)=W(n)%tTyz1(i,j,k)+b*W(n)%hTyz1(i,j,k)
end do
end do
end do
end if
if (W(n)%isy) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)% Vx2(i,j,k)=W(n)% tVx2(i,j,k)+b*W(n)% hVx2(i,j,k)
   W(n)% Vy2(i,j,k)=W(n)% tVy2(i,j,k)+b*W(n)% hVy2(i,j,k)
   W(n)% Vz2(i,j,k)=W(n)% tVz2(i,j,k)+b*W(n)% hVz2(i,j,k)
   W(n)%Txx2(i,j,k)=W(n)%tTxx2(i,j,k)+b*W(n)%hTxx2(i,j,k)
   W(n)%Tyy2(i,j,k)=W(n)%tTyy2(i,j,k)+b*W(n)%hTyy2(i,j,k)
   W(n)%Tzz2(i,j,k)=W(n)%tTzz2(i,j,k)+b*W(n)%hTzz2(i,j,k)
   W(n)%Txy2(i,j,k)=W(n)%tTxy2(i,j,k)+b*W(n)%hTxy2(i,j,k)
   W(n)%Txz2(i,j,k)=W(n)%tTxz2(i,j,k)+b*W(n)%hTxz2(i,j,k)
   W(n)%Tyz2(i,j,k)=W(n)%tTyz2(i,j,k)+b*W(n)%hTyz2(i,j,k)
end do
end do
end do
end if
if (W(n)%isz) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   W(n)% Vx3(i,j,k)=W(n)% tVx3(i,j,k)+b*W(n)% hVx3(i,j,k)
   W(n)% Vy3(i,j,k)=W(n)% tVy3(i,j,k)+b*W(n)% hVy3(i,j,k)
   W(n)% Vz3(i,j,k)=W(n)% tVz3(i,j,k)+b*W(n)% hVz3(i,j,k)
   W(n)%Txx3(i,j,k)=W(n)%tTxx3(i,j,k)+b*W(n)%hTxx3(i,j,k)
   W(n)%Tyy3(i,j,k)=W(n)%tTyy3(i,j,k)+b*W(n)%hTyy3(i,j,k)
   W(n)%Tzz3(i,j,k)=W(n)%tTzz3(i,j,k)+b*W(n)%hTzz3(i,j,k)
   W(n)%Txy3(i,j,k)=W(n)%tTxy3(i,j,k)+b*W(n)%hTxy3(i,j,k)
   W(n)%Txz3(i,j,k)=W(n)%tTxz3(i,j,k)+b*W(n)%hTxz3(i,j,k)
   W(n)%Tyz3(i,j,k)=W(n)%tTyz3(i,j,k)+b*W(n)%hTyz3(i,j,k)
end do
end do
end do
end if
end do
!$OMP END PARALLEL DO
!if (freenode) then
!do n=1,num_blk
!if ( W(n)%nk2==nk2) then
!   k=W(n)%nk2
!
!if ( W(n)%isx ) then
!do j=W(n)%nj1,W(n)%nj2
!do i=W(n)%ni1,W(n)%ni2
!   W(n)%  Vx1(i,j,k)=0.0_SP
!   W(n)%  Vy1(i,j,k)=0.0_SP
!   W(n)%  Vz1(i,j,k)=0.0_SP
!   W(n)% Txx1(i,j,k)=0.0_SP
!   W(n)% Tyy1(i,j,k)=0.0_SP
!   W(n)% Tzz1(i,j,k)=0.0_SP
!   W(n)% Txy1(i,j,k)=0.0_SP
!   W(n)% Txz1(i,j,k)=0.0_SP
!   W(n)% Tyz1(i,j,k)=0.0_SP
!end do
!end do
!end if !isx
!if ( W(n)%isy ) then
!do j=W(n)%nj1,W(n)%nj2
!do i=W(n)%ni1,W(n)%ni2
!   W(n)%  Vx2(i,j,k)=0.0_SP
!   W(n)%  Vy2(i,j,k)=0.0_SP
!   W(n)%  Vz2(i,j,k)=0.0_SP
!   W(n)% Txx2(i,j,k)=0.0_SP
!   W(n)% Tyy2(i,j,k)=0.0_SP
!   W(n)% Tzz2(i,j,k)=0.0_SP
!   W(n)% Txy2(i,j,k)=0.0_SP
!   W(n)% Txz2(i,j,k)=0.0_SP
!   W(n)% Tyz2(i,j,k)=0.0_SP
!end do
!end do
!end if !isy
!
!if ( (W(n)%isx .or. W(n)%isy) ) then
!do j=W(n)%nj1,W(n)%nj2
!do i=W(n)%ni1,W(n)%ni2
!   e=min(Ex(i),Ey(j),Ez(k))
!    Vx (i,j,k)= Vx (i,j,k)*e
!    Vy (i,j,k)= Vy (i,j,k)*e
!    Vz (i,j,k)= Vz (i,j,k)*e
!   Txx (i,j,k)=Txx (i,j,k)*e
!   Tyy (i,j,k)=Tyy (i,j,k)*e
!   Tzz (i,j,k)=Tzz (i,j,k)*e
!   Txy (i,j,k)=Txy (i,j,k)*e
!   Txz (i,j,k)=Txz (i,j,k)*e
!   Tyz (i,j,k)=Tyz (i,j,k)*e
!end do
!end do
!end if !damp
!
!end if !nk2
!end do !nblk
!end if !freenode
end subroutine abs_RK_fin

!-----------------------------------------------------------------------------
subroutine in_LxF_LyF_LzF
  call LxF_LyF_LzF
#ifndef CorrAbs
  call LxF_LyF_LzF_TIMG
  call LxF_LyF_LzF_VHOC
#endif
end subroutine in_LxF_LyF_LzF
subroutine in_LxB_LyB_LzB
  call LxB_LyB_LzB
#ifndef CorrAbs
  call LxB_LyB_LzB_TIMG
  call LxB_LyB_LzB_VHOC
#endif
end subroutine in_LxB_LyB_LzB
subroutine in_LxF_LyF_LzB
  call LxF_LyF_LzB
#ifndef CorrAbs
  call LxF_LyF_LzB_TIMG
  call LxF_LyF_LzB_VHOC
#endif
end subroutine in_LxF_LyF_LzB
subroutine in_LxB_LyB_LzF
  call LxB_LyB_LzF
#ifndef CorrAbs
  call LxB_LyB_LzF_TIMG
  call LxB_LyB_LzF_VHOC
#endif
end subroutine in_LxB_LyB_LzF
subroutine in_LxB_LyF_LzF
  call LxB_LyF_LzF
#ifndef CorrAbs
  call LxB_LyF_LzF_TIMG
  call LxB_LyF_LzF_VHOC
#endif
end subroutine in_LxB_LyF_LzF
subroutine in_LxF_LyB_LzB
  call LxF_LyB_LzB
#ifndef CorrAbs
  call LxF_LyB_LzB_TIMG
  call LxF_LyB_LzB_VHOC
#endif
end subroutine in_LxF_LyB_LzB
subroutine in_LxF_LyB_LzF
  call LxF_LyB_LzF
#ifndef CorrAbs
  call LxF_LyB_LzF_TIMG
  call LxF_LyB_LzF_VHOC
#endif
end subroutine in_LxF_LyB_LzF
subroutine in_LxB_LyF_LzB
  call LxB_LyF_LzB
#ifndef CorrAbs
  call LxB_LyF_LzB_TIMG
  call LxB_LyF_LzB_VHOC
#endif
end subroutine in_LxB_LyF_LzB

!-----------------------------------------------------------------------------

subroutine LxF_LyF_LzF
integer :: n,i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: DzVx1,DzVx2,DzVy1,DzVy2,DzVz1,DzVz2
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef CorrAbs
call abs_tsymm
call abs_vsymm
#endif

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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

#ifdef AbsVzero
if (freenode .and. k==nk2) then
   !call coef_fdxy2fdz_abs(i,j)
   DzVx1= matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz
   DzVx2= matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy1= matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz
   DzVy2= matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz1= matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz 
   DzVz2= matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx = DzVx1 + DzVx2
   DzVy = DzVy1 + DzVy2
   DzVz = DzVz1 + DzVz2
   !if (W(n)%isx) then
   !   DzVx=DzVx+MulD(1,1)*W(n)%Txz1(i,j,k)/b1
   !   DzVy=DzVy+MulD(1,2)*W(n)%Tyz1(i,j,k)/b1
   !   DzVz=DzVz+MulD(1,3)*W(n)%Tzz1(i,j,k)/b1
   !end if
   !if (W(n)%isy) then
   !   DzVx=DzVx+MulD(2,1)*W(n)%Txz2(i,j,k) /b2
   !   DzVy=DzVy+MulD(2,2)*W(n)%Tyz2(i,j,k) /b2
   !   DzVz=DzVz+MulD(2,3)*W(n)%Tzz2(i,j,k) /b2
   !end if
   !if (W(n)%isz) then
   !   DzVx=DzVx+MulD(3,1)*W(n)%Txz3(i,j,k)/b3
   !   DzVy=DzVy+MulD(3,2)*W(n)%Tyz3(i,j,k)/b3
   !   DzVz=DzVz+MulD(3,3)*W(n)%Tzz3(i,j,k)/b3
   !end if
end if
#endif

#ifndef CorrAbs
    hVx(i,j,k)= rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )/b1 &
              + rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )/b2 &
              + rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )/b3
    hVy(i,j,k)= rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )/b1 &
              + rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )/b2 &
              + rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )/b3
    hVz(i,j,k)= rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )/b1 &
              + rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )/b2 &
              + rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )/b3
   hTxx(i,j,k)= ( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz )/b3
   hTxy(i,j,k)= miu *(                                          &
                ( xiy*DxVx + xix*DxVy )/b1                      &
               +( ety*DyVx + etx*DyVy )/b2                      &
               +( zty*DzVx + ztx*DzVy )/b3                      &
               )
   hTxz(i,j,k)= miu *(                                          &
                ( xiz*DxVx + xix*DxVz )/b1                      &
               +( etz*DyVx + etx*DyVz )/b2                      &
               +( ztz*DzVx + ztx*DzVz )/b3                      &
               )
   hTyz(i,j,k)= miu *(                                          &
                ( xiz*DxVy + xiy*DxVz )/b1                      &
               +( etz*DyVy + ety*DyVz )/b2                      &
               +( ztz*DzVy + zty*DzVz )/b3                      &
               )
#endif
if (W(n)%isx) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx1(i,j,k)/b1
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy1(i,j,k)/b1
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz1(i,j,k)/b1
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz1(i,j,k)/b1
   W(n)% hVx1(i,j,k)= d1*rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )  &
               - (d1+a1)*W(n)% Vx1(i,j,k)
   W(n)% hVy1(i,j,k)= d1*rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )  &
               - (d1+a1)*W(n)% Vy1(i,j,k)  
   W(n)% hVz1(i,j,k)= d1*rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )  &
               - (d1+a1)*W(n)% Vz1(i,j,k)  
   W(n)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Txx1(i,j,k)  
   W(n)%hTyy1(i,j,k)= d1*( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tyy1(i,j,k)  
   W(n)%hTzz1(i,j,k)= d1*( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tzz1(i,j,k) 
   W(n)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx + xix*DxVy ) &
               - (d1+a1)*W(n)%Txy1(i,j,k)  
   W(n)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx + xix*DxVz ) &
               - (d1+a1)*W(n)%Txz1(i,j,k)  
   W(n)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy + xiy*DxVz ) &
               - (d1+a1)*W(n)%Tyz1(i,j,k)  
end if
if (W(n)%isy) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx2(i,j,k)/b2
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy2(i,j,k)/b2
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz2(i,j,k)/b2
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )  &
               - (d2+a2)*W(n)%Vx2(i,j,k) 
   W(n)%hVy2(i,j,k)= d2*rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )  &
               - (d2+a2)*W(n)%Vy2(i,j,k)  
   W(n)%hVz2(i,j,k)= d2*rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )  &
               - (d2+a2)*W(n)%Vz2(i,j,k)  
   W(n)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Txx2(i,j,k)  
   W(n)%hTyy2(i,j,k)= d2*( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Tyy2(i,j,k)  
   W(n)%hTzz2(i,j,k)= d2*( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz ) &
               - (d2+a2)*W(n)%Tzz2(i,j,k) 
   W(n)%hTxy2(i,j,k)= d2*miu*( ety*DyVx + etx*DyVy ) &
               - (d2+a2)*W(n)%Txy2(i,j,k)  
   W(n)%hTxz2(i,j,k)= d2*miu*( etz*DyVx + etx*DyVz ) &
               - (d2+a2)*W(n)%Txz2(i,j,k)  
   W(n)%hTyz2(i,j,k)= d2*miu*( etz*DyVy + ety*DyVz ) &
               - (d2+a2)*W(n)%Tyz2(i,j,k)  
end if
if (W(n)%isz) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx3(i,j,k)/b3
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy3(i,j,k)/b3
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz3(i,j,k)/b3
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )  &
               - (d3+a3)*W(n)%Vx3(i,j,k) 
   W(n)%hVy3(i,j,k)= d3*rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )  &
               - (d3+a3)*W(n)%Vy3(i,j,k)  
   W(n)%hVz3(i,j,k)= d3*rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )  &
               - (d3+a3)*W(n)%Vz3(i,j,k)  
   W(n)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Txx3(i,j,k)  
   W(n)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tyy3(i,j,k)  
   W(n)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tzz3(i,j,k) 
   W(n)%hTxy3(i,j,k)= d3*miu*( zty*DzVx + ztx*DzVy ) &
               - (d3+a3)*W(n)%Txy3(i,j,k)  
   W(n)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx + ztx*DzVz ) &
               - (d3+a3)*W(n)%Txz3(i,j,k) 
   W(n)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy + zty*DzVz ) &
               - (d3+a3)*W(n)%Tyz3(i,j,k) 
end if

!if (freenode .and. k==nk2) then
!if (W(n)%isx) then
!   W(n)%hTxz1(i,j,k)=0.0_SP
!   W(n)%hTyz1(i,j,k)=0.0_SP
!   W(n)%hTzz1(i,j,k)=0.0_SP
!   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx 
!   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!if (W(n)%isy) then
!   W(n)%hTxz2(i,j,k)=0.0_SP
!   W(n)%hTyz2(i,j,k)=0.0_SP
!   W(n)%hTzz2(i,j,k)=0.0_SP
!   W(n)%hTxx2(i,j,k)= W(n)%hTxx2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   W(n)%hTyy2(i,j,k)= W(n)%hTyy2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!end if
#ifdef AbsVzero
if (freenode .and. k==nk2) then
if (W(n)%isx) then
   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) &
               + d1 * ( lam2mu*ztx*DzVx1 + lam*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam2mu*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTzz1(i,j,k)=W(n)%hTzz1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam*zty*DzVy1 + lam2mu*ztz*DzVz1 ) *b1
   W(n)%hTyz1(i,j,k)=W(n)%hTyz1(i,j,k) &
               + d1 * miu * ( ztz*DzVy1 + zty*DzVz1 ) *b1
   W(n)%hTxz1(i,j,k)=W(n)%hTxz1(i,j,k) &
               + d1 * miu * ( ztz*DzVx1 + ztx*DzVz1 ) *b1
   W(n)%hTxy1(i,j,k)=W(n)%hTxy1(i,j,k) &
               + d1 * miu * ( zty*DzVx1 + ztx*DzVy1 ) *b1
end if
if (W(n)%isy) then
   W(n)%hTxx2(i,j,k)=W(n)%hTxx2(i,j,k) &
               + d2 * ( lam2mu*ztx*DzVx2 + lam*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTyy2(i,j,k)=W(n)%hTyy2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam2mu*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTzz2(i,j,k)=W(n)%hTzz2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam*zty*DzVy2 + lam2mu*ztz*DzVz2 ) *b2
   W(n)%hTyz2(i,j,k)=W(n)%hTyz2(i,j,k) &
               + d2 * miu * ( ztz*DzVy2 + zty*DzVz2 ) *b2
   W(n)%hTxz2(i,j,k)=W(n)%hTxz2(i,j,k) &
               + d2 * miu * ( ztz*DzVx2 + ztx*DzVz2 ) *b2
   W(n)%hTxy2(i,j,k)=W(n)%hTxy2(i,j,k) &
               + d2 * miu * ( zty*DzVx2 + ztx*DzVy2 ) *b2
end if
end if
#endif

end do
end do
end do
end do
end subroutine LxF_LyF_LzF

subroutine LxB_LyB_LzB
integer :: n,i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: DzVx1,DzVx2,DzVy1,DzVy2,DzVz1,DzVz2
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef CorrAbs
call abs_tsymm
call abs_vsymm
#endif

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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

#ifdef AbsVzero
if (freenode .and. k==nk2) then
   !call coef_fdxy2fdz_abs(i,j)
   DzVx1= matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz
   DzVx2= matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy1= matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz
   DzVy2= matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz1= matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz 
   DzVz2= matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx = DzVx1 + DzVx2
   DzVy = DzVy1 + DzVy2
   DzVz = DzVz1 + DzVz2
   !if (W(n)%isx) then
   !   DzVx=DzVx+MulD(1,1)*W(n)%Txz1(i,j,k)/b1
   !   DzVy=DzVy+MulD(1,2)*W(n)%Tyz1(i,j,k)/b1
   !   DzVz=DzVz+MulD(1,3)*W(n)%Tzz1(i,j,k)/b1
   !end if
   !if (W(n)%isy) then
   !   DzVx=DzVx+MulD(2,1)*W(n)%Txz2(i,j,k) /b2
   !   DzVy=DzVy+MulD(2,2)*W(n)%Tyz2(i,j,k) /b2
   !   DzVz=DzVz+MulD(2,3)*W(n)%Tzz2(i,j,k) /b2
   !end if
   !if (W(n)%isz) then
   !   DzVx=DzVx+MulD(3,1)*W(n)%Txz3(i,j,k)/b3
   !   DzVy=DzVy+MulD(3,2)*W(n)%Tyz3(i,j,k)/b3
   !   DzVz=DzVz+MulD(3,3)*W(n)%Tzz3(i,j,k)/b3
   !end if
end if
#endif

#ifndef CorrAbs
    hVx(i,j,k)= rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )/b1 &
              + rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )/b2 &
              + rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )/b3
    hVy(i,j,k)= rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )/b1 &
              + rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )/b2 &
              + rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )/b3
    hVz(i,j,k)= rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )/b1 &
              + rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )/b2 &
              + rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )/b3
   hTxx(i,j,k)= ( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz )/b3
   hTxy(i,j,k)= miu *(                                          &
                ( xiy*DxVx + xix*DxVy )/b1                      &
               +( ety*DyVx + etx*DyVy )/b2                      &
               +( zty*DzVx + ztx*DzVy )/b3                      &
               )
   hTxz(i,j,k)= miu *(                                          &
                ( xiz*DxVx + xix*DxVz )/b1                      &
               +( etz*DyVx + etx*DyVz )/b2                      &
               +( ztz*DzVx + ztx*DzVz )/b3                      &
               )
   hTyz(i,j,k)= miu *(                                          &
                ( xiz*DxVy + xiy*DxVz )/b1                      &
               +( etz*DyVy + ety*DyVz )/b2                      &
               +( ztz*DzVy + zty*DzVz )/b3                      &
               )
#endif
if (W(n)%isx) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx1(i,j,k)/b1
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy1(i,j,k)/b1
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz1(i,j,k)/b1
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz1(i,j,k)/b1
   W(n)% hVx1(i,j,k)= d1*rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )  &
               - (d1+a1)*W(n)% Vx1(i,j,k)
   W(n)% hVy1(i,j,k)= d1*rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )  &
               - (d1+a1)*W(n)% Vy1(i,j,k)  
   W(n)% hVz1(i,j,k)= d1*rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )  &
               - (d1+a1)*W(n)% Vz1(i,j,k)  
   W(n)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Txx1(i,j,k)  
   W(n)%hTyy1(i,j,k)= d1*( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tyy1(i,j,k)  
   W(n)%hTzz1(i,j,k)= d1*( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tzz1(i,j,k) 
   W(n)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx + xix*DxVy ) &
               - (d1+a1)*W(n)%Txy1(i,j,k)  
   W(n)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx + xix*DxVz ) &
               - (d1+a1)*W(n)%Txz1(i,j,k)  
   W(n)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy + xiy*DxVz ) &
               - (d1+a1)*W(n)%Tyz1(i,j,k)  
end if
if (W(n)%isy) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx2(i,j,k)/b2
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy2(i,j,k)/b2
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz2(i,j,k)/b2
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )  &
               - (d2+a2)*W(n)%Vx2(i,j,k) 
   W(n)%hVy2(i,j,k)= d2*rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )  &
               - (d2+a2)*W(n)%Vy2(i,j,k)  
   W(n)%hVz2(i,j,k)= d2*rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )  &
               - (d2+a2)*W(n)%Vz2(i,j,k)  
   W(n)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Txx2(i,j,k)  
   W(n)%hTyy2(i,j,k)= d2*( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Tyy2(i,j,k)  
   W(n)%hTzz2(i,j,k)= d2*( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz ) &
               - (d2+a2)*W(n)%Tzz2(i,j,k) 
   W(n)%hTxy2(i,j,k)= d2*miu*( ety*DyVx + etx*DyVy ) &
               - (d2+a2)*W(n)%Txy2(i,j,k)  
   W(n)%hTxz2(i,j,k)= d2*miu*( etz*DyVx + etx*DyVz ) &
               - (d2+a2)*W(n)%Txz2(i,j,k)  
   W(n)%hTyz2(i,j,k)= d2*miu*( etz*DyVy + ety*DyVz ) &
               - (d2+a2)*W(n)%Tyz2(i,j,k)  
end if
if (W(n)%isz) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx3(i,j,k)/b3
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy3(i,j,k)/b3
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz3(i,j,k)/b3
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )  &
               - (d3+a3)*W(n)%Vx3(i,j,k) 
   W(n)%hVy3(i,j,k)= d3*rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )  &
               - (d3+a3)*W(n)%Vy3(i,j,k)  
   W(n)%hVz3(i,j,k)= d3*rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )  &
               - (d3+a3)*W(n)%Vz3(i,j,k)  
   W(n)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Txx3(i,j,k)  
   W(n)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tyy3(i,j,k)  
   W(n)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tzz3(i,j,k) 
   W(n)%hTxy3(i,j,k)= d3*miu*( zty*DzVx + ztx*DzVy ) &
               - (d3+a3)*W(n)%Txy3(i,j,k)  
   W(n)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx + ztx*DzVz ) &
               - (d3+a3)*W(n)%Txz3(i,j,k) 
   W(n)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy + zty*DzVz ) &
               - (d3+a3)*W(n)%Tyz3(i,j,k) 
end if

!if (freenode .and. k==nk2) then
!if (W(n)%isx) then
!   W(n)%hTxz1(i,j,k)=0.0_SP
!   W(n)%hTyz1(i,j,k)=0.0_SP
!   W(n)%hTzz1(i,j,k)=0.0_SP
!   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx 
!   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!if (W(n)%isy) then
!   W(n)%hTxz2(i,j,k)=0.0_SP
!   W(n)%hTyz2(i,j,k)=0.0_SP
!   W(n)%hTzz2(i,j,k)=0.0_SP
!   W(n)%hTxx2(i,j,k)= W(n)%hTxx2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   W(n)%hTyy2(i,j,k)= W(n)%hTyy2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!end if
#ifdef AbsVzero
if (freenode .and. k==nk2) then
if (W(n)%isx) then
   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) &
               + d1 * ( lam2mu*ztx*DzVx1 + lam*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam2mu*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTzz1(i,j,k)=W(n)%hTzz1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam*zty*DzVy1 + lam2mu*ztz*DzVz1 ) *b1
   W(n)%hTyz1(i,j,k)=W(n)%hTyz1(i,j,k) &
               + d1 * miu * ( ztz*DzVy1 + zty*DzVz1 ) *b1
   W(n)%hTxz1(i,j,k)=W(n)%hTxz1(i,j,k) &
               + d1 * miu * ( ztz*DzVx1 + ztx*DzVz1 ) *b1
   W(n)%hTxy1(i,j,k)=W(n)%hTxy1(i,j,k) &
               + d1 * miu * ( zty*DzVx1 + ztx*DzVy1 ) *b1
end if
if (W(n)%isy) then
   W(n)%hTxx2(i,j,k)=W(n)%hTxx2(i,j,k) &
               + d2 * ( lam2mu*ztx*DzVx2 + lam*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTyy2(i,j,k)=W(n)%hTyy2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam2mu*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTzz2(i,j,k)=W(n)%hTzz2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam*zty*DzVy2 + lam2mu*ztz*DzVz2 ) *b2
   W(n)%hTyz2(i,j,k)=W(n)%hTyz2(i,j,k) &
               + d2 * miu * ( ztz*DzVy2 + zty*DzVz2 ) *b2
   W(n)%hTxz2(i,j,k)=W(n)%hTxz2(i,j,k) &
               + d2 * miu * ( ztz*DzVx2 + ztx*DzVz2 ) *b2
   W(n)%hTxy2(i,j,k)=W(n)%hTxy2(i,j,k) &
               + d2 * miu * ( zty*DzVx2 + ztx*DzVy2 ) *b2
end if
end if
#endif

end do
end do
end do
end do
end subroutine LxB_LyB_LzB

subroutine LxF_LyF_LzB
integer :: n,i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: DzVx1,DzVx2,DzVy1,DzVy2,DzVz1,DzVz2
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef CorrAbs
call abs_tsymm
call abs_vsymm
#endif

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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

#ifdef AbsVzero
if (freenode .and. k==nk2) then
   !call coef_fdxy2fdz_abs(i,j)
   DzVx1= matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz
   DzVx2= matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy1= matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz
   DzVy2= matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz1= matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz 
   DzVz2= matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx = DzVx1 + DzVx2
   DzVy = DzVy1 + DzVy2
   DzVz = DzVz1 + DzVz2
   !if (W(n)%isx) then
   !   DzVx=DzVx+MulD(1,1)*W(n)%Txz1(i,j,k)/b1
   !   DzVy=DzVy+MulD(1,2)*W(n)%Tyz1(i,j,k)/b1
   !   DzVz=DzVz+MulD(1,3)*W(n)%Tzz1(i,j,k)/b1
   !end if
   !if (W(n)%isy) then
   !   DzVx=DzVx+MulD(2,1)*W(n)%Txz2(i,j,k) /b2
   !   DzVy=DzVy+MulD(2,2)*W(n)%Tyz2(i,j,k) /b2
   !   DzVz=DzVz+MulD(2,3)*W(n)%Tzz2(i,j,k) /b2
   !end if
   !if (W(n)%isz) then
   !   DzVx=DzVx+MulD(3,1)*W(n)%Txz3(i,j,k)/b3
   !   DzVy=DzVy+MulD(3,2)*W(n)%Tyz3(i,j,k)/b3
   !   DzVz=DzVz+MulD(3,3)*W(n)%Tzz3(i,j,k)/b3
   !end if
end if
#endif

#ifndef CorrAbs
    hVx(i,j,k)= rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )/b1 &
              + rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )/b2 &
              + rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )/b3
    hVy(i,j,k)= rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )/b1 &
              + rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )/b2 &
              + rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )/b3
    hVz(i,j,k)= rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )/b1 &
              + rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )/b2 &
              + rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )/b3
   hTxx(i,j,k)= ( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz )/b3
   hTxy(i,j,k)= miu *(                                          &
                ( xiy*DxVx + xix*DxVy )/b1                      &
               +( ety*DyVx + etx*DyVy )/b2                      &
               +( zty*DzVx + ztx*DzVy )/b3                      &
               )
   hTxz(i,j,k)= miu *(                                          &
                ( xiz*DxVx + xix*DxVz )/b1                      &
               +( etz*DyVx + etx*DyVz )/b2                      &
               +( ztz*DzVx + ztx*DzVz )/b3                      &
               )
   hTyz(i,j,k)= miu *(                                          &
                ( xiz*DxVy + xiy*DxVz )/b1                      &
               +( etz*DyVy + ety*DyVz )/b2                      &
               +( ztz*DzVy + zty*DzVz )/b3                      &
               )
#endif
if (W(n)%isx) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx1(i,j,k)/b1
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy1(i,j,k)/b1
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz1(i,j,k)/b1
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz1(i,j,k)/b1
   W(n)% hVx1(i,j,k)= d1*rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )  &
               - (d1+a1)*W(n)% Vx1(i,j,k)
   W(n)% hVy1(i,j,k)= d1*rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )  &
               - (d1+a1)*W(n)% Vy1(i,j,k)  
   W(n)% hVz1(i,j,k)= d1*rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )  &
               - (d1+a1)*W(n)% Vz1(i,j,k)  
   W(n)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Txx1(i,j,k)  
   W(n)%hTyy1(i,j,k)= d1*( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tyy1(i,j,k)  
   W(n)%hTzz1(i,j,k)= d1*( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tzz1(i,j,k) 
   W(n)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx + xix*DxVy ) &
               - (d1+a1)*W(n)%Txy1(i,j,k)  
   W(n)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx + xix*DxVz ) &
               - (d1+a1)*W(n)%Txz1(i,j,k)  
   W(n)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy + xiy*DxVz ) &
               - (d1+a1)*W(n)%Tyz1(i,j,k)  
end if
if (W(n)%isy) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx2(i,j,k)/b2
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy2(i,j,k)/b2
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz2(i,j,k)/b2
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )  &
               - (d2+a2)*W(n)%Vx2(i,j,k) 
   W(n)%hVy2(i,j,k)= d2*rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )  &
               - (d2+a2)*W(n)%Vy2(i,j,k)  
   W(n)%hVz2(i,j,k)= d2*rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )  &
               - (d2+a2)*W(n)%Vz2(i,j,k)  
   W(n)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Txx2(i,j,k)  
   W(n)%hTyy2(i,j,k)= d2*( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Tyy2(i,j,k)  
   W(n)%hTzz2(i,j,k)= d2*( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz ) &
               - (d2+a2)*W(n)%Tzz2(i,j,k) 
   W(n)%hTxy2(i,j,k)= d2*miu*( ety*DyVx + etx*DyVy ) &
               - (d2+a2)*W(n)%Txy2(i,j,k)  
   W(n)%hTxz2(i,j,k)= d2*miu*( etz*DyVx + etx*DyVz ) &
               - (d2+a2)*W(n)%Txz2(i,j,k)  
   W(n)%hTyz2(i,j,k)= d2*miu*( etz*DyVy + ety*DyVz ) &
               - (d2+a2)*W(n)%Tyz2(i,j,k)  
end if
if (W(n)%isz) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx3(i,j,k)/b3
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy3(i,j,k)/b3
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz3(i,j,k)/b3
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )  &
               - (d3+a3)*W(n)%Vx3(i,j,k) 
   W(n)%hVy3(i,j,k)= d3*rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )  &
               - (d3+a3)*W(n)%Vy3(i,j,k)  
   W(n)%hVz3(i,j,k)= d3*rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )  &
               - (d3+a3)*W(n)%Vz3(i,j,k)  
   W(n)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Txx3(i,j,k)  
   W(n)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tyy3(i,j,k)  
   W(n)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tzz3(i,j,k) 
   W(n)%hTxy3(i,j,k)= d3*miu*( zty*DzVx + ztx*DzVy ) &
               - (d3+a3)*W(n)%Txy3(i,j,k)  
   W(n)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx + ztx*DzVz ) &
               - (d3+a3)*W(n)%Txz3(i,j,k) 
   W(n)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy + zty*DzVz ) &
               - (d3+a3)*W(n)%Tyz3(i,j,k) 
end if

!if (freenode .and. k==nk2) then
!if (W(n)%isx) then
!   W(n)%hTxz1(i,j,k)=0.0_SP
!   W(n)%hTyz1(i,j,k)=0.0_SP
!   W(n)%hTzz1(i,j,k)=0.0_SP
!   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx 
!   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!if (W(n)%isy) then
!   W(n)%hTxz2(i,j,k)=0.0_SP
!   W(n)%hTyz2(i,j,k)=0.0_SP
!   W(n)%hTzz2(i,j,k)=0.0_SP
!   W(n)%hTxx2(i,j,k)= W(n)%hTxx2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   W(n)%hTyy2(i,j,k)= W(n)%hTyy2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!end if
#ifdef AbsVzero
if (freenode .and. k==nk2) then
if (W(n)%isx) then
   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) &
               + d1 * ( lam2mu*ztx*DzVx1 + lam*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam2mu*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTzz1(i,j,k)=W(n)%hTzz1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam*zty*DzVy1 + lam2mu*ztz*DzVz1 ) *b1
   W(n)%hTyz1(i,j,k)=W(n)%hTyz1(i,j,k) &
               + d1 * miu * ( ztz*DzVy1 + zty*DzVz1 ) *b1
   W(n)%hTxz1(i,j,k)=W(n)%hTxz1(i,j,k) &
               + d1 * miu * ( ztz*DzVx1 + ztx*DzVz1 ) *b1
   W(n)%hTxy1(i,j,k)=W(n)%hTxy1(i,j,k) &
               + d1 * miu * ( zty*DzVx1 + ztx*DzVy1 ) *b1
end if
if (W(n)%isy) then
   W(n)%hTxx2(i,j,k)=W(n)%hTxx2(i,j,k) &
               + d2 * ( lam2mu*ztx*DzVx2 + lam*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTyy2(i,j,k)=W(n)%hTyy2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam2mu*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTzz2(i,j,k)=W(n)%hTzz2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam*zty*DzVy2 + lam2mu*ztz*DzVz2 ) *b2
   W(n)%hTyz2(i,j,k)=W(n)%hTyz2(i,j,k) &
               + d2 * miu * ( ztz*DzVy2 + zty*DzVz2 ) *b2
   W(n)%hTxz2(i,j,k)=W(n)%hTxz2(i,j,k) &
               + d2 * miu * ( ztz*DzVx2 + ztx*DzVz2 ) *b2
   W(n)%hTxy2(i,j,k)=W(n)%hTxy2(i,j,k) &
               + d2 * miu * ( zty*DzVx2 + ztx*DzVy2 ) *b2
end if
end if
#endif

end do
end do
end do
end do
end subroutine LxF_LyF_LzB

subroutine LxB_LyB_LzF
integer :: n,i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: DzVx1,DzVx2,DzVy1,DzVy2,DzVz1,DzVz2
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef CorrAbs
call abs_tsymm
call abs_vsymm
#endif

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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

#ifdef AbsVzero
if (freenode .and. k==nk2) then
   !call coef_fdxy2fdz_abs(i,j)
   DzVx1= matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz
   DzVx2= matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy1= matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz
   DzVy2= matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz1= matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz 
   DzVz2= matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx = DzVx1 + DzVx2
   DzVy = DzVy1 + DzVy2
   DzVz = DzVz1 + DzVz2
   !if (W(n)%isx) then
   !   DzVx=DzVx+MulD(1,1)*W(n)%Txz1(i,j,k)/b1
   !   DzVy=DzVy+MulD(1,2)*W(n)%Tyz1(i,j,k)/b1
   !   DzVz=DzVz+MulD(1,3)*W(n)%Tzz1(i,j,k)/b1
   !end if
   !if (W(n)%isy) then
   !   DzVx=DzVx+MulD(2,1)*W(n)%Txz2(i,j,k) /b2
   !   DzVy=DzVy+MulD(2,2)*W(n)%Tyz2(i,j,k) /b2
   !   DzVz=DzVz+MulD(2,3)*W(n)%Tzz2(i,j,k) /b2
   !end if
   !if (W(n)%isz) then
   !   DzVx=DzVx+MulD(3,1)*W(n)%Txz3(i,j,k)/b3
   !   DzVy=DzVy+MulD(3,2)*W(n)%Tyz3(i,j,k)/b3
   !   DzVz=DzVz+MulD(3,3)*W(n)%Tzz3(i,j,k)/b3
   !end if
end if
#endif

#ifndef CorrAbs
    hVx(i,j,k)= rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )/b1 &
              + rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )/b2 &
              + rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )/b3
    hVy(i,j,k)= rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )/b1 &
              + rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )/b2 &
              + rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )/b3
    hVz(i,j,k)= rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )/b1 &
              + rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )/b2 &
              + rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )/b3
   hTxx(i,j,k)= ( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz )/b3
   hTxy(i,j,k)= miu *(                                          &
                ( xiy*DxVx + xix*DxVy )/b1                      &
               +( ety*DyVx + etx*DyVy )/b2                      &
               +( zty*DzVx + ztx*DzVy )/b3                      &
               )
   hTxz(i,j,k)= miu *(                                          &
                ( xiz*DxVx + xix*DxVz )/b1                      &
               +( etz*DyVx + etx*DyVz )/b2                      &
               +( ztz*DzVx + ztx*DzVz )/b3                      &
               )
   hTyz(i,j,k)= miu *(                                          &
                ( xiz*DxVy + xiy*DxVz )/b1                      &
               +( etz*DyVy + ety*DyVz )/b2                      &
               +( ztz*DzVy + zty*DzVz )/b3                      &
               )
#endif
if (W(n)%isx) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx1(i,j,k)/b1
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy1(i,j,k)/b1
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz1(i,j,k)/b1
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz1(i,j,k)/b1
   W(n)% hVx1(i,j,k)= d1*rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )  &
               - (d1+a1)*W(n)% Vx1(i,j,k)
   W(n)% hVy1(i,j,k)= d1*rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )  &
               - (d1+a1)*W(n)% Vy1(i,j,k)  
   W(n)% hVz1(i,j,k)= d1*rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )  &
               - (d1+a1)*W(n)% Vz1(i,j,k)  
   W(n)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Txx1(i,j,k)  
   W(n)%hTyy1(i,j,k)= d1*( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tyy1(i,j,k)  
   W(n)%hTzz1(i,j,k)= d1*( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tzz1(i,j,k) 
   W(n)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx + xix*DxVy ) &
               - (d1+a1)*W(n)%Txy1(i,j,k)  
   W(n)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx + xix*DxVz ) &
               - (d1+a1)*W(n)%Txz1(i,j,k)  
   W(n)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy + xiy*DxVz ) &
               - (d1+a1)*W(n)%Tyz1(i,j,k)  
end if
if (W(n)%isy) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx2(i,j,k)/b2
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy2(i,j,k)/b2
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz2(i,j,k)/b2
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )  &
               - (d2+a2)*W(n)%Vx2(i,j,k) 
   W(n)%hVy2(i,j,k)= d2*rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )  &
               - (d2+a2)*W(n)%Vy2(i,j,k)  
   W(n)%hVz2(i,j,k)= d2*rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )  &
               - (d2+a2)*W(n)%Vz2(i,j,k)  
   W(n)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Txx2(i,j,k)  
   W(n)%hTyy2(i,j,k)= d2*( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Tyy2(i,j,k)  
   W(n)%hTzz2(i,j,k)= d2*( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz ) &
               - (d2+a2)*W(n)%Tzz2(i,j,k) 
   W(n)%hTxy2(i,j,k)= d2*miu*( ety*DyVx + etx*DyVy ) &
               - (d2+a2)*W(n)%Txy2(i,j,k)  
   W(n)%hTxz2(i,j,k)= d2*miu*( etz*DyVx + etx*DyVz ) &
               - (d2+a2)*W(n)%Txz2(i,j,k)  
   W(n)%hTyz2(i,j,k)= d2*miu*( etz*DyVy + ety*DyVz ) &
               - (d2+a2)*W(n)%Tyz2(i,j,k)  
end if
if (W(n)%isz) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx3(i,j,k)/b3
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy3(i,j,k)/b3
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz3(i,j,k)/b3
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )  &
               - (d3+a3)*W(n)%Vx3(i,j,k) 
   W(n)%hVy3(i,j,k)= d3*rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )  &
               - (d3+a3)*W(n)%Vy3(i,j,k)  
   W(n)%hVz3(i,j,k)= d3*rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )  &
               - (d3+a3)*W(n)%Vz3(i,j,k)  
   W(n)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Txx3(i,j,k)  
   W(n)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tyy3(i,j,k)  
   W(n)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tzz3(i,j,k) 
   W(n)%hTxy3(i,j,k)= d3*miu*( zty*DzVx + ztx*DzVy ) &
               - (d3+a3)*W(n)%Txy3(i,j,k)  
   W(n)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx + ztx*DzVz ) &
               - (d3+a3)*W(n)%Txz3(i,j,k) 
   W(n)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy + zty*DzVz ) &
               - (d3+a3)*W(n)%Tyz3(i,j,k) 
end if

!if (freenode .and. k==nk2) then
!if (W(n)%isx) then
!   W(n)%hTxz1(i,j,k)=0.0_SP
!   W(n)%hTyz1(i,j,k)=0.0_SP
!   W(n)%hTzz1(i,j,k)=0.0_SP
!   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx 
!   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!if (W(n)%isy) then
!   W(n)%hTxz2(i,j,k)=0.0_SP
!   W(n)%hTyz2(i,j,k)=0.0_SP
!   W(n)%hTzz2(i,j,k)=0.0_SP
!   W(n)%hTxx2(i,j,k)= W(n)%hTxx2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   W(n)%hTyy2(i,j,k)= W(n)%hTyy2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!end if
#ifdef AbsVzero
if (freenode .and. k==nk2) then
if (W(n)%isx) then
   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) &
               + d1 * ( lam2mu*ztx*DzVx1 + lam*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam2mu*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTzz1(i,j,k)=W(n)%hTzz1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam*zty*DzVy1 + lam2mu*ztz*DzVz1 ) *b1
   W(n)%hTyz1(i,j,k)=W(n)%hTyz1(i,j,k) &
               + d1 * miu * ( ztz*DzVy1 + zty*DzVz1 ) *b1
   W(n)%hTxz1(i,j,k)=W(n)%hTxz1(i,j,k) &
               + d1 * miu * ( ztz*DzVx1 + ztx*DzVz1 ) *b1
   W(n)%hTxy1(i,j,k)=W(n)%hTxy1(i,j,k) &
               + d1 * miu * ( zty*DzVx1 + ztx*DzVy1 ) *b1
end if
if (W(n)%isy) then
   W(n)%hTxx2(i,j,k)=W(n)%hTxx2(i,j,k) &
               + d2 * ( lam2mu*ztx*DzVx2 + lam*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTyy2(i,j,k)=W(n)%hTyy2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam2mu*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTzz2(i,j,k)=W(n)%hTzz2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam*zty*DzVy2 + lam2mu*ztz*DzVz2 ) *b2
   W(n)%hTyz2(i,j,k)=W(n)%hTyz2(i,j,k) &
               + d2 * miu * ( ztz*DzVy2 + zty*DzVz2 ) *b2
   W(n)%hTxz2(i,j,k)=W(n)%hTxz2(i,j,k) &
               + d2 * miu * ( ztz*DzVx2 + ztx*DzVz2 ) *b2
   W(n)%hTxy2(i,j,k)=W(n)%hTxy2(i,j,k) &
               + d2 * miu * ( zty*DzVx2 + ztx*DzVy2 ) *b2
end if
end if
#endif

end do
end do
end do
end do
end subroutine LxB_LyB_LzF

subroutine LxB_LyF_LzF
integer :: n,i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: DzVx1,DzVx2,DzVy1,DzVy2,DzVz1,DzVz2
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef CorrAbs
call abs_tsymm
call abs_vsymm
#endif

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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

#ifdef AbsVzero
if (freenode .and. k==nk2) then
   !call coef_fdxy2fdz_abs(i,j)
   DzVx1= matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz
   DzVx2= matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy1= matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz
   DzVy2= matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz1= matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz 
   DzVz2= matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx = DzVx1 + DzVx2
   DzVy = DzVy1 + DzVy2
   DzVz = DzVz1 + DzVz2
   !if (W(n)%isx) then
   !   DzVx=DzVx+MulD(1,1)*W(n)%Txz1(i,j,k)/b1
   !   DzVy=DzVy+MulD(1,2)*W(n)%Tyz1(i,j,k)/b1
   !   DzVz=DzVz+MulD(1,3)*W(n)%Tzz1(i,j,k)/b1
   !end if
   !if (W(n)%isy) then
   !   DzVx=DzVx+MulD(2,1)*W(n)%Txz2(i,j,k) /b2
   !   DzVy=DzVy+MulD(2,2)*W(n)%Tyz2(i,j,k) /b2
   !   DzVz=DzVz+MulD(2,3)*W(n)%Tzz2(i,j,k) /b2
   !end if
   !if (W(n)%isz) then
   !   DzVx=DzVx+MulD(3,1)*W(n)%Txz3(i,j,k)/b3
   !   DzVy=DzVy+MulD(3,2)*W(n)%Tyz3(i,j,k)/b3
   !   DzVz=DzVz+MulD(3,3)*W(n)%Tzz3(i,j,k)/b3
   !end if
end if
#endif

#ifndef CorrAbs
    hVx(i,j,k)= rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )/b1 &
              + rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )/b2 &
              + rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )/b3
    hVy(i,j,k)= rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )/b1 &
              + rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )/b2 &
              + rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )/b3
    hVz(i,j,k)= rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )/b1 &
              + rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )/b2 &
              + rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )/b3
   hTxx(i,j,k)= ( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz )/b3
   hTxy(i,j,k)= miu *(                                          &
                ( xiy*DxVx + xix*DxVy )/b1                      &
               +( ety*DyVx + etx*DyVy )/b2                      &
               +( zty*DzVx + ztx*DzVy )/b3                      &
               )
   hTxz(i,j,k)= miu *(                                          &
                ( xiz*DxVx + xix*DxVz )/b1                      &
               +( etz*DyVx + etx*DyVz )/b2                      &
               +( ztz*DzVx + ztx*DzVz )/b3                      &
               )
   hTyz(i,j,k)= miu *(                                          &
                ( xiz*DxVy + xiy*DxVz )/b1                      &
               +( etz*DyVy + ety*DyVz )/b2                      &
               +( ztz*DzVy + zty*DzVz )/b3                      &
               )
#endif
if (W(n)%isx) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx1(i,j,k)/b1
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy1(i,j,k)/b1
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz1(i,j,k)/b1
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz1(i,j,k)/b1
   W(n)% hVx1(i,j,k)= d1*rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )  &
               - (d1+a1)*W(n)% Vx1(i,j,k)
   W(n)% hVy1(i,j,k)= d1*rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )  &
               - (d1+a1)*W(n)% Vy1(i,j,k)  
   W(n)% hVz1(i,j,k)= d1*rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )  &
               - (d1+a1)*W(n)% Vz1(i,j,k)  
   W(n)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Txx1(i,j,k)  
   W(n)%hTyy1(i,j,k)= d1*( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tyy1(i,j,k)  
   W(n)%hTzz1(i,j,k)= d1*( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tzz1(i,j,k) 
   W(n)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx + xix*DxVy ) &
               - (d1+a1)*W(n)%Txy1(i,j,k)  
   W(n)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx + xix*DxVz ) &
               - (d1+a1)*W(n)%Txz1(i,j,k)  
   W(n)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy + xiy*DxVz ) &
               - (d1+a1)*W(n)%Tyz1(i,j,k)  
end if
if (W(n)%isy) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx2(i,j,k)/b2
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy2(i,j,k)/b2
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz2(i,j,k)/b2
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )  &
               - (d2+a2)*W(n)%Vx2(i,j,k) 
   W(n)%hVy2(i,j,k)= d2*rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )  &
               - (d2+a2)*W(n)%Vy2(i,j,k)  
   W(n)%hVz2(i,j,k)= d2*rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )  &
               - (d2+a2)*W(n)%Vz2(i,j,k)  
   W(n)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Txx2(i,j,k)  
   W(n)%hTyy2(i,j,k)= d2*( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Tyy2(i,j,k)  
   W(n)%hTzz2(i,j,k)= d2*( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz ) &
               - (d2+a2)*W(n)%Tzz2(i,j,k) 
   W(n)%hTxy2(i,j,k)= d2*miu*( ety*DyVx + etx*DyVy ) &
               - (d2+a2)*W(n)%Txy2(i,j,k)  
   W(n)%hTxz2(i,j,k)= d2*miu*( etz*DyVx + etx*DyVz ) &
               - (d2+a2)*W(n)%Txz2(i,j,k)  
   W(n)%hTyz2(i,j,k)= d2*miu*( etz*DyVy + ety*DyVz ) &
               - (d2+a2)*W(n)%Tyz2(i,j,k)  
end if
if (W(n)%isz) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx3(i,j,k)/b3
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy3(i,j,k)/b3
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz3(i,j,k)/b3
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )  &
               - (d3+a3)*W(n)%Vx3(i,j,k) 
   W(n)%hVy3(i,j,k)= d3*rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )  &
               - (d3+a3)*W(n)%Vy3(i,j,k)  
   W(n)%hVz3(i,j,k)= d3*rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )  &
               - (d3+a3)*W(n)%Vz3(i,j,k)  
   W(n)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Txx3(i,j,k)  
   W(n)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tyy3(i,j,k)  
   W(n)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tzz3(i,j,k) 
   W(n)%hTxy3(i,j,k)= d3*miu*( zty*DzVx + ztx*DzVy ) &
               - (d3+a3)*W(n)%Txy3(i,j,k)  
   W(n)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx + ztx*DzVz ) &
               - (d3+a3)*W(n)%Txz3(i,j,k) 
   W(n)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy + zty*DzVz ) &
               - (d3+a3)*W(n)%Tyz3(i,j,k) 
end if

!if (freenode .and. k==nk2) then
!if (W(n)%isx) then
!   W(n)%hTxz1(i,j,k)=0.0_SP
!   W(n)%hTyz1(i,j,k)=0.0_SP
!   W(n)%hTzz1(i,j,k)=0.0_SP
!   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx 
!   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!if (W(n)%isy) then
!   W(n)%hTxz2(i,j,k)=0.0_SP
!   W(n)%hTyz2(i,j,k)=0.0_SP
!   W(n)%hTzz2(i,j,k)=0.0_SP
!   W(n)%hTxx2(i,j,k)= W(n)%hTxx2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   W(n)%hTyy2(i,j,k)= W(n)%hTyy2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!end if
#ifdef AbsVzero
if (freenode .and. k==nk2) then
if (W(n)%isx) then
   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) &
               + d1 * ( lam2mu*ztx*DzVx1 + lam*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam2mu*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTzz1(i,j,k)=W(n)%hTzz1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam*zty*DzVy1 + lam2mu*ztz*DzVz1 ) *b1
   W(n)%hTyz1(i,j,k)=W(n)%hTyz1(i,j,k) &
               + d1 * miu * ( ztz*DzVy1 + zty*DzVz1 ) *b1
   W(n)%hTxz1(i,j,k)=W(n)%hTxz1(i,j,k) &
               + d1 * miu * ( ztz*DzVx1 + ztx*DzVz1 ) *b1
   W(n)%hTxy1(i,j,k)=W(n)%hTxy1(i,j,k) &
               + d1 * miu * ( zty*DzVx1 + ztx*DzVy1 ) *b1
end if
if (W(n)%isy) then
   W(n)%hTxx2(i,j,k)=W(n)%hTxx2(i,j,k) &
               + d2 * ( lam2mu*ztx*DzVx2 + lam*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTyy2(i,j,k)=W(n)%hTyy2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam2mu*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTzz2(i,j,k)=W(n)%hTzz2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam*zty*DzVy2 + lam2mu*ztz*DzVz2 ) *b2
   W(n)%hTyz2(i,j,k)=W(n)%hTyz2(i,j,k) &
               + d2 * miu * ( ztz*DzVy2 + zty*DzVz2 ) *b2
   W(n)%hTxz2(i,j,k)=W(n)%hTxz2(i,j,k) &
               + d2 * miu * ( ztz*DzVx2 + ztx*DzVz2 ) *b2
   W(n)%hTxy2(i,j,k)=W(n)%hTxy2(i,j,k) &
               + d2 * miu * ( zty*DzVx2 + ztx*DzVy2 ) *b2
end if
end if
#endif

end do
end do
end do
end do
end subroutine LxB_LyF_LzF

subroutine LxF_LyB_LzB
integer :: n,i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: DzVx1,DzVx2,DzVy1,DzVy2,DzVz1,DzVz2
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef CorrAbs
call abs_tsymm
call abs_vsymm
#endif

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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

#ifdef AbsVzero
if (freenode .and. k==nk2) then
   !call coef_fdxy2fdz_abs(i,j)
   DzVx1= matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz
   DzVx2= matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy1= matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz
   DzVy2= matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz1= matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz 
   DzVz2= matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx = DzVx1 + DzVx2
   DzVy = DzVy1 + DzVy2
   DzVz = DzVz1 + DzVz2
   !if (W(n)%isx) then
   !   DzVx=DzVx+MulD(1,1)*W(n)%Txz1(i,j,k)/b1
   !   DzVy=DzVy+MulD(1,2)*W(n)%Tyz1(i,j,k)/b1
   !   DzVz=DzVz+MulD(1,3)*W(n)%Tzz1(i,j,k)/b1
   !end if
   !if (W(n)%isy) then
   !   DzVx=DzVx+MulD(2,1)*W(n)%Txz2(i,j,k) /b2
   !   DzVy=DzVy+MulD(2,2)*W(n)%Tyz2(i,j,k) /b2
   !   DzVz=DzVz+MulD(2,3)*W(n)%Tzz2(i,j,k) /b2
   !end if
   !if (W(n)%isz) then
   !   DzVx=DzVx+MulD(3,1)*W(n)%Txz3(i,j,k)/b3
   !   DzVy=DzVy+MulD(3,2)*W(n)%Tyz3(i,j,k)/b3
   !   DzVz=DzVz+MulD(3,3)*W(n)%Tzz3(i,j,k)/b3
   !end if
end if
#endif

#ifndef CorrAbs
    hVx(i,j,k)= rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )/b1 &
              + rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )/b2 &
              + rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )/b3
    hVy(i,j,k)= rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )/b1 &
              + rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )/b2 &
              + rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )/b3
    hVz(i,j,k)= rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )/b1 &
              + rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )/b2 &
              + rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )/b3
   hTxx(i,j,k)= ( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz )/b3
   hTxy(i,j,k)= miu *(                                          &
                ( xiy*DxVx + xix*DxVy )/b1                      &
               +( ety*DyVx + etx*DyVy )/b2                      &
               +( zty*DzVx + ztx*DzVy )/b3                      &
               )
   hTxz(i,j,k)= miu *(                                          &
                ( xiz*DxVx + xix*DxVz )/b1                      &
               +( etz*DyVx + etx*DyVz )/b2                      &
               +( ztz*DzVx + ztx*DzVz )/b3                      &
               )
   hTyz(i,j,k)= miu *(                                          &
                ( xiz*DxVy + xiy*DxVz )/b1                      &
               +( etz*DyVy + ety*DyVz )/b2                      &
               +( ztz*DzVy + zty*DzVz )/b3                      &
               )
#endif
if (W(n)%isx) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx1(i,j,k)/b1
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy1(i,j,k)/b1
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz1(i,j,k)/b1
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz1(i,j,k)/b1
   W(n)% hVx1(i,j,k)= d1*rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )  &
               - (d1+a1)*W(n)% Vx1(i,j,k)
   W(n)% hVy1(i,j,k)= d1*rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )  &
               - (d1+a1)*W(n)% Vy1(i,j,k)  
   W(n)% hVz1(i,j,k)= d1*rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )  &
               - (d1+a1)*W(n)% Vz1(i,j,k)  
   W(n)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Txx1(i,j,k)  
   W(n)%hTyy1(i,j,k)= d1*( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tyy1(i,j,k)  
   W(n)%hTzz1(i,j,k)= d1*( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tzz1(i,j,k) 
   W(n)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx + xix*DxVy ) &
               - (d1+a1)*W(n)%Txy1(i,j,k)  
   W(n)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx + xix*DxVz ) &
               - (d1+a1)*W(n)%Txz1(i,j,k)  
   W(n)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy + xiy*DxVz ) &
               - (d1+a1)*W(n)%Tyz1(i,j,k)  
end if
if (W(n)%isy) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx2(i,j,k)/b2
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy2(i,j,k)/b2
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz2(i,j,k)/b2
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )  &
               - (d2+a2)*W(n)%Vx2(i,j,k) 
   W(n)%hVy2(i,j,k)= d2*rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )  &
               - (d2+a2)*W(n)%Vy2(i,j,k)  
   W(n)%hVz2(i,j,k)= d2*rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )  &
               - (d2+a2)*W(n)%Vz2(i,j,k)  
   W(n)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Txx2(i,j,k)  
   W(n)%hTyy2(i,j,k)= d2*( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Tyy2(i,j,k)  
   W(n)%hTzz2(i,j,k)= d2*( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz ) &
               - (d2+a2)*W(n)%Tzz2(i,j,k) 
   W(n)%hTxy2(i,j,k)= d2*miu*( ety*DyVx + etx*DyVy ) &
               - (d2+a2)*W(n)%Txy2(i,j,k)  
   W(n)%hTxz2(i,j,k)= d2*miu*( etz*DyVx + etx*DyVz ) &
               - (d2+a2)*W(n)%Txz2(i,j,k)  
   W(n)%hTyz2(i,j,k)= d2*miu*( etz*DyVy + ety*DyVz ) &
               - (d2+a2)*W(n)%Tyz2(i,j,k)  
end if
if (W(n)%isz) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx3(i,j,k)/b3
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy3(i,j,k)/b3
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz3(i,j,k)/b3
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )  &
               - (d3+a3)*W(n)%Vx3(i,j,k) 
   W(n)%hVy3(i,j,k)= d3*rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )  &
               - (d3+a3)*W(n)%Vy3(i,j,k)  
   W(n)%hVz3(i,j,k)= d3*rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )  &
               - (d3+a3)*W(n)%Vz3(i,j,k)  
   W(n)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Txx3(i,j,k)  
   W(n)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tyy3(i,j,k)  
   W(n)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tzz3(i,j,k) 
   W(n)%hTxy3(i,j,k)= d3*miu*( zty*DzVx + ztx*DzVy ) &
               - (d3+a3)*W(n)%Txy3(i,j,k)  
   W(n)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx + ztx*DzVz ) &
               - (d3+a3)*W(n)%Txz3(i,j,k) 
   W(n)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy + zty*DzVz ) &
               - (d3+a3)*W(n)%Tyz3(i,j,k) 
end if

!if (freenode .and. k==nk2) then
!if (W(n)%isx) then
!   W(n)%hTxz1(i,j,k)=0.0_SP
!   W(n)%hTyz1(i,j,k)=0.0_SP
!   W(n)%hTzz1(i,j,k)=0.0_SP
!   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx 
!   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!if (W(n)%isy) then
!   W(n)%hTxz2(i,j,k)=0.0_SP
!   W(n)%hTyz2(i,j,k)=0.0_SP
!   W(n)%hTzz2(i,j,k)=0.0_SP
!   W(n)%hTxx2(i,j,k)= W(n)%hTxx2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   W(n)%hTyy2(i,j,k)= W(n)%hTyy2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!end if
#ifdef AbsVzero
if (freenode .and. k==nk2) then
if (W(n)%isx) then
   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) &
               + d1 * ( lam2mu*ztx*DzVx1 + lam*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam2mu*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTzz1(i,j,k)=W(n)%hTzz1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam*zty*DzVy1 + lam2mu*ztz*DzVz1 ) *b1
   W(n)%hTyz1(i,j,k)=W(n)%hTyz1(i,j,k) &
               + d1 * miu * ( ztz*DzVy1 + zty*DzVz1 ) *b1
   W(n)%hTxz1(i,j,k)=W(n)%hTxz1(i,j,k) &
               + d1 * miu * ( ztz*DzVx1 + ztx*DzVz1 ) *b1
   W(n)%hTxy1(i,j,k)=W(n)%hTxy1(i,j,k) &
               + d1 * miu * ( zty*DzVx1 + ztx*DzVy1 ) *b1
end if
if (W(n)%isy) then
   W(n)%hTxx2(i,j,k)=W(n)%hTxx2(i,j,k) &
               + d2 * ( lam2mu*ztx*DzVx2 + lam*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTyy2(i,j,k)=W(n)%hTyy2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam2mu*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTzz2(i,j,k)=W(n)%hTzz2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam*zty*DzVy2 + lam2mu*ztz*DzVz2 ) *b2
   W(n)%hTyz2(i,j,k)=W(n)%hTyz2(i,j,k) &
               + d2 * miu * ( ztz*DzVy2 + zty*DzVz2 ) *b2
   W(n)%hTxz2(i,j,k)=W(n)%hTxz2(i,j,k) &
               + d2 * miu * ( ztz*DzVx2 + ztx*DzVz2 ) *b2
   W(n)%hTxy2(i,j,k)=W(n)%hTxy2(i,j,k) &
               + d2 * miu * ( zty*DzVx2 + ztx*DzVy2 ) *b2
end if
end if
#endif

end do
end do
end do
end do
end subroutine LxF_LyB_LzB

subroutine LxF_LyB_LzF
integer :: n,i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: DzVx1,DzVx2,DzVy1,DzVy2,DzVz1,DzVz2
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef CorrAbs
call abs_tsymm
call abs_vsymm
#endif

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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

#ifdef AbsVzero
if (freenode .and. k==nk2) then
   !call coef_fdxy2fdz_abs(i,j)
   DzVx1= matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz
   DzVx2= matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy1= matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz
   DzVy2= matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz1= matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz 
   DzVz2= matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx = DzVx1 + DzVx2
   DzVy = DzVy1 + DzVy2
   DzVz = DzVz1 + DzVz2
   !if (W(n)%isx) then
   !   DzVx=DzVx+MulD(1,1)*W(n)%Txz1(i,j,k)/b1
   !   DzVy=DzVy+MulD(1,2)*W(n)%Tyz1(i,j,k)/b1
   !   DzVz=DzVz+MulD(1,3)*W(n)%Tzz1(i,j,k)/b1
   !end if
   !if (W(n)%isy) then
   !   DzVx=DzVx+MulD(2,1)*W(n)%Txz2(i,j,k) /b2
   !   DzVy=DzVy+MulD(2,2)*W(n)%Tyz2(i,j,k) /b2
   !   DzVz=DzVz+MulD(2,3)*W(n)%Tzz2(i,j,k) /b2
   !end if
   !if (W(n)%isz) then
   !   DzVx=DzVx+MulD(3,1)*W(n)%Txz3(i,j,k)/b3
   !   DzVy=DzVy+MulD(3,2)*W(n)%Tyz3(i,j,k)/b3
   !   DzVz=DzVz+MulD(3,3)*W(n)%Tzz3(i,j,k)/b3
   !end if
end if
#endif

#ifndef CorrAbs
    hVx(i,j,k)= rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )/b1 &
              + rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )/b2 &
              + rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )/b3
    hVy(i,j,k)= rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )/b1 &
              + rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )/b2 &
              + rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )/b3
    hVz(i,j,k)= rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )/b1 &
              + rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )/b2 &
              + rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )/b3
   hTxx(i,j,k)= ( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz )/b3
   hTxy(i,j,k)= miu *(                                          &
                ( xiy*DxVx + xix*DxVy )/b1                      &
               +( ety*DyVx + etx*DyVy )/b2                      &
               +( zty*DzVx + ztx*DzVy )/b3                      &
               )
   hTxz(i,j,k)= miu *(                                          &
                ( xiz*DxVx + xix*DxVz )/b1                      &
               +( etz*DyVx + etx*DyVz )/b2                      &
               +( ztz*DzVx + ztx*DzVz )/b3                      &
               )
   hTyz(i,j,k)= miu *(                                          &
                ( xiz*DxVy + xiy*DxVz )/b1                      &
               +( etz*DyVy + ety*DyVz )/b2                      &
               +( ztz*DzVy + zty*DzVz )/b3                      &
               )
#endif
if (W(n)%isx) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx1(i,j,k)/b1
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy1(i,j,k)/b1
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz1(i,j,k)/b1
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz1(i,j,k)/b1
   W(n)% hVx1(i,j,k)= d1*rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )  &
               - (d1+a1)*W(n)% Vx1(i,j,k)
   W(n)% hVy1(i,j,k)= d1*rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )  &
               - (d1+a1)*W(n)% Vy1(i,j,k)  
   W(n)% hVz1(i,j,k)= d1*rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )  &
               - (d1+a1)*W(n)% Vz1(i,j,k)  
   W(n)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Txx1(i,j,k)  
   W(n)%hTyy1(i,j,k)= d1*( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tyy1(i,j,k)  
   W(n)%hTzz1(i,j,k)= d1*( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tzz1(i,j,k) 
   W(n)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx + xix*DxVy ) &
               - (d1+a1)*W(n)%Txy1(i,j,k)  
   W(n)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx + xix*DxVz ) &
               - (d1+a1)*W(n)%Txz1(i,j,k)  
   W(n)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy + xiy*DxVz ) &
               - (d1+a1)*W(n)%Tyz1(i,j,k)  
end if
if (W(n)%isy) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx2(i,j,k)/b2
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy2(i,j,k)/b2
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz2(i,j,k)/b2
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )  &
               - (d2+a2)*W(n)%Vx2(i,j,k) 
   W(n)%hVy2(i,j,k)= d2*rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )  &
               - (d2+a2)*W(n)%Vy2(i,j,k)  
   W(n)%hVz2(i,j,k)= d2*rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )  &
               - (d2+a2)*W(n)%Vz2(i,j,k)  
   W(n)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Txx2(i,j,k)  
   W(n)%hTyy2(i,j,k)= d2*( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Tyy2(i,j,k)  
   W(n)%hTzz2(i,j,k)= d2*( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz ) &
               - (d2+a2)*W(n)%Tzz2(i,j,k) 
   W(n)%hTxy2(i,j,k)= d2*miu*( ety*DyVx + etx*DyVy ) &
               - (d2+a2)*W(n)%Txy2(i,j,k)  
   W(n)%hTxz2(i,j,k)= d2*miu*( etz*DyVx + etx*DyVz ) &
               - (d2+a2)*W(n)%Txz2(i,j,k)  
   W(n)%hTyz2(i,j,k)= d2*miu*( etz*DyVy + ety*DyVz ) &
               - (d2+a2)*W(n)%Tyz2(i,j,k)  
end if
if (W(n)%isz) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx3(i,j,k)/b3
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy3(i,j,k)/b3
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz3(i,j,k)/b3
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )  &
               - (d3+a3)*W(n)%Vx3(i,j,k) 
   W(n)%hVy3(i,j,k)= d3*rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )  &
               - (d3+a3)*W(n)%Vy3(i,j,k)  
   W(n)%hVz3(i,j,k)= d3*rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )  &
               - (d3+a3)*W(n)%Vz3(i,j,k)  
   W(n)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Txx3(i,j,k)  
   W(n)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tyy3(i,j,k)  
   W(n)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tzz3(i,j,k) 
   W(n)%hTxy3(i,j,k)= d3*miu*( zty*DzVx + ztx*DzVy ) &
               - (d3+a3)*W(n)%Txy3(i,j,k)  
   W(n)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx + ztx*DzVz ) &
               - (d3+a3)*W(n)%Txz3(i,j,k) 
   W(n)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy + zty*DzVz ) &
               - (d3+a3)*W(n)%Tyz3(i,j,k) 
end if

!if (freenode .and. k==nk2) then
!if (W(n)%isx) then
!   W(n)%hTxz1(i,j,k)=0.0_SP
!   W(n)%hTyz1(i,j,k)=0.0_SP
!   W(n)%hTzz1(i,j,k)=0.0_SP
!   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx 
!   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!if (W(n)%isy) then
!   W(n)%hTxz2(i,j,k)=0.0_SP
!   W(n)%hTyz2(i,j,k)=0.0_SP
!   W(n)%hTzz2(i,j,k)=0.0_SP
!   W(n)%hTxx2(i,j,k)= W(n)%hTxx2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   W(n)%hTyy2(i,j,k)= W(n)%hTyy2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!end if
#ifdef AbsVzero
if (freenode .and. k==nk2) then
if (W(n)%isx) then
   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) &
               + d1 * ( lam2mu*ztx*DzVx1 + lam*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam2mu*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTzz1(i,j,k)=W(n)%hTzz1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam*zty*DzVy1 + lam2mu*ztz*DzVz1 ) *b1
   W(n)%hTyz1(i,j,k)=W(n)%hTyz1(i,j,k) &
               + d1 * miu * ( ztz*DzVy1 + zty*DzVz1 ) *b1
   W(n)%hTxz1(i,j,k)=W(n)%hTxz1(i,j,k) &
               + d1 * miu * ( ztz*DzVx1 + ztx*DzVz1 ) *b1
   W(n)%hTxy1(i,j,k)=W(n)%hTxy1(i,j,k) &
               + d1 * miu * ( zty*DzVx1 + ztx*DzVy1 ) *b1
end if
if (W(n)%isy) then
   W(n)%hTxx2(i,j,k)=W(n)%hTxx2(i,j,k) &
               + d2 * ( lam2mu*ztx*DzVx2 + lam*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTyy2(i,j,k)=W(n)%hTyy2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam2mu*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTzz2(i,j,k)=W(n)%hTzz2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam*zty*DzVy2 + lam2mu*ztz*DzVz2 ) *b2
   W(n)%hTyz2(i,j,k)=W(n)%hTyz2(i,j,k) &
               + d2 * miu * ( ztz*DzVy2 + zty*DzVz2 ) *b2
   W(n)%hTxz2(i,j,k)=W(n)%hTxz2(i,j,k) &
               + d2 * miu * ( ztz*DzVx2 + ztx*DzVz2 ) *b2
   W(n)%hTxy2(i,j,k)=W(n)%hTxy2(i,j,k) &
               + d2 * miu * ( zty*DzVx2 + ztx*DzVy2 ) *b2
end if
end if
#endif

end do
end do
end do
end do
end subroutine LxF_LyB_LzF

subroutine LxB_LyF_LzB
integer :: n,i,j,k
real(SP) :: DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz
real(SP) :: DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz
real(SP) :: DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz
real(SP) :: DzVx1,DzVx2,DzVy1,DzVy2,DzVz1,DzVz2
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef CorrAbs
call abs_tsymm
call abs_vsymm
#endif

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=1.0/rho(i,j,k)
   xix=  xi_x(i,j,k); xiy=  xi_y(i,j,k); xiz=  xi_z(i,j,k)
   etx= eta_x(i,j,k); ety= eta_y(i,j,k); etz= eta_z(i,j,k)
   ztx=zeta_x(i,j,k); zty=zeta_y(i,j,k); ztz=zeta_z(i,j,k)
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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

#ifdef AbsVzero
if (freenode .and. k==nk2) then
   !call coef_fdxy2fdz_abs(i,j)
   DzVx1= matVx2Vz(1,1,i,j)*DxVx &
         +matVx2Vz(1,2,i,j)*DxVy &
         +matVx2Vz(1,3,i,j)*DxVz
   DzVx2= matVy2Vz(1,1,i,j)*DyVx &
         +matVy2Vz(1,2,i,j)*DyVy &
         +matVy2Vz(1,3,i,j)*DyVz
   DzVy1= matVx2Vz(2,1,i,j)*DxVx &
         +matVx2Vz(2,2,i,j)*DxVy &
         +matVx2Vz(2,3,i,j)*DxVz
   DzVy2= matVy2Vz(2,1,i,j)*DyVx &
         +matVy2Vz(2,2,i,j)*DyVy &
         +matVy2Vz(2,3,i,j)*DyVz
   DzVz1= matVx2Vz(3,1,i,j)*DxVx &
         +matVx2Vz(3,2,i,j)*DxVy &
         +matVx2Vz(3,3,i,j)*DxVz 
   DzVz2= matVy2Vz(3,1,i,j)*DyVx &
         +matVy2Vz(3,2,i,j)*DyVy &
         +matVy2Vz(3,3,i,j)*DyVz
   DzVx = DzVx1 + DzVx2
   DzVy = DzVy1 + DzVy2
   DzVz = DzVz1 + DzVz2
   !if (W(n)%isx) then
   !   DzVx=DzVx+MulD(1,1)*W(n)%Txz1(i,j,k)/b1
   !   DzVy=DzVy+MulD(1,2)*W(n)%Tyz1(i,j,k)/b1
   !   DzVz=DzVz+MulD(1,3)*W(n)%Tzz1(i,j,k)/b1
   !end if
   !if (W(n)%isy) then
   !   DzVx=DzVx+MulD(2,1)*W(n)%Txz2(i,j,k) /b2
   !   DzVy=DzVy+MulD(2,2)*W(n)%Tyz2(i,j,k) /b2
   !   DzVz=DzVz+MulD(2,3)*W(n)%Tzz2(i,j,k) /b2
   !end if
   !if (W(n)%isz) then
   !   DzVx=DzVx+MulD(3,1)*W(n)%Txz3(i,j,k)/b3
   !   DzVy=DzVy+MulD(3,2)*W(n)%Tyz3(i,j,k)/b3
   !   DzVz=DzVz+MulD(3,3)*W(n)%Tzz3(i,j,k)/b3
   !end if
end if
#endif

#ifndef CorrAbs
    hVx(i,j,k)= rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )/b1 &
              + rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )/b2 &
              + rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )/b3
    hVy(i,j,k)= rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )/b1 &
              + rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )/b2 &
              + rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )/b3
    hVz(i,j,k)= rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )/b1 &
              + rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )/b2 &
              + rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )/b3
   hTxx(i,j,k)= ( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz )/b1 &
               +( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz )/b2 &
               +( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz )/b3
   hTxy(i,j,k)= miu *(                                          &
                ( xiy*DxVx + xix*DxVy )/b1                      &
               +( ety*DyVx + etx*DyVy )/b2                      &
               +( zty*DzVx + ztx*DzVy )/b3                      &
               )
   hTxz(i,j,k)= miu *(                                          &
                ( xiz*DxVx + xix*DxVz )/b1                      &
               +( etz*DyVx + etx*DyVz )/b2                      &
               +( ztz*DzVx + ztx*DzVz )/b3                      &
               )
   hTyz(i,j,k)= miu *(                                          &
                ( xiz*DxVy + xiy*DxVz )/b1                      &
               +( etz*DyVy + ety*DyVz )/b2                      &
               +( ztz*DzVy + zty*DzVz )/b3                      &
               )
#endif
if (W(n)%isx) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx1(i,j,k)/b1
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy1(i,j,k)/b1
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz1(i,j,k)/b1
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz1(i,j,k)/b1
   W(n)% hVx1(i,j,k)= d1*rrho*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz )  &
               - (d1+a1)*W(n)% Vx1(i,j,k)
   W(n)% hVy1(i,j,k)= d1*rrho*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz )  &
               - (d1+a1)*W(n)% Vy1(i,j,k)  
   W(n)% hVz1(i,j,k)= d1*rrho*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz )  &
               - (d1+a1)*W(n)% Vz1(i,j,k)  
   W(n)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Txx1(i,j,k)  
   W(n)%hTyy1(i,j,k)= d1*( lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tyy1(i,j,k)  
   W(n)%hTzz1(i,j,k)= d1*( lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz ) &
               - (d1+a1)*W(n)%Tzz1(i,j,k) 
   W(n)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx + xix*DxVy ) &
               - (d1+a1)*W(n)%Txy1(i,j,k)  
   W(n)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx + xix*DxVz ) &
               - (d1+a1)*W(n)%Txz1(i,j,k)  
   W(n)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy + xiy*DxVz ) &
               - (d1+a1)*W(n)%Tyz1(i,j,k)  
end if
if (W(n)%isy) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx2(i,j,k)/b2
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy2(i,j,k)/b2
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz2(i,j,k)/b2
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrho*( etx*DyTxx + ety*DyTxy + etz*DyTxz )  &
               - (d2+a2)*W(n)%Vx2(i,j,k) 
   W(n)%hVy2(i,j,k)= d2*rrho*( etx*DyTxy + ety*DyTyy + etz*DyTyz )  &
               - (d2+a2)*W(n)%Vy2(i,j,k)  
   W(n)%hVz2(i,j,k)= d2*rrho*( etx*DyTxz + ety*DyTyz + etz*DyTzz )  &
               - (d2+a2)*W(n)%Vz2(i,j,k)  
   W(n)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Txx2(i,j,k)  
   W(n)%hTyy2(i,j,k)= d2*( lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz ) &
               - (d2+a2)*W(n)%Tyy2(i,j,k)  
   W(n)%hTzz2(i,j,k)= d2*( lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz ) &
               - (d2+a2)*W(n)%Tzz2(i,j,k) 
   W(n)%hTxy2(i,j,k)= d2*miu*( ety*DyVx + etx*DyVy ) &
               - (d2+a2)*W(n)%Txy2(i,j,k)  
   W(n)%hTxz2(i,j,k)= d2*miu*( etz*DyVx + etx*DyVz ) &
               - (d2+a2)*W(n)%Txz2(i,j,k)  
   W(n)%hTyz2(i,j,k)= d2*miu*( etz*DyVy + ety*DyVz ) &
               - (d2+a2)*W(n)%Tyz2(i,j,k)  
end if
if (W(n)%isz) then
    hVx(i,j,k)= hVx(i,j,k)-W(n)% Vx3(i,j,k)/b3
    hVy(i,j,k)= hVy(i,j,k)-W(n)% Vy3(i,j,k)/b3
    hVz(i,j,k)= hVz(i,j,k)-W(n)% Vz3(i,j,k)/b3
   hTxx(i,j,k)=hTxx(i,j,k)-W(n)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(n)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(n)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(n)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(n)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(n)%Tyz3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrho*( ztx*DzTxx + zty*DzTxy + ztz*DzTxz )  &
               - (d3+a3)*W(n)%Vx3(i,j,k) 
   W(n)%hVy3(i,j,k)= d3*rrho*( ztx*DzTxy + zty*DzTyy + ztz*DzTyz )  &
               - (d3+a3)*W(n)%Vy3(i,j,k)  
   W(n)%hVz3(i,j,k)= d3*rrho*( ztx*DzTxz + zty*DzTyz + ztz*DzTzz )  &
               - (d3+a3)*W(n)%Vz3(i,j,k)  
   W(n)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Txx3(i,j,k)  
   W(n)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tyy3(i,j,k)  
   W(n)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz ) &
               - (d3+a3)*W(n)%Tzz3(i,j,k) 
   W(n)%hTxy3(i,j,k)= d3*miu*( zty*DzVx + ztx*DzVy ) &
               - (d3+a3)*W(n)%Txy3(i,j,k)  
   W(n)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx + ztx*DzVz ) &
               - (d3+a3)*W(n)%Txz3(i,j,k) 
   W(n)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy + zty*DzVz ) &
               - (d3+a3)*W(n)%Tyz3(i,j,k) 
end if

!if (freenode .and. k==nk2) then
!if (W(n)%isx) then
!   W(n)%hTxz1(i,j,k)=0.0_SP
!   W(n)%hTyz1(i,j,k)=0.0_SP
!   W(n)%hTzz1(i,j,k)=0.0_SP
!   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx 
!   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) - d1* lam**2/lam2mu *xix *DxVx
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!if (W(n)%isy) then
!   W(n)%hTxz2(i,j,k)=0.0_SP
!   W(n)%hTyz2(i,j,k)=0.0_SP
!   W(n)%hTzz2(i,j,k)=0.0_SP
!   W(n)%hTxx2(i,j,k)= W(n)%hTxx2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   W(n)%hTyy2(i,j,k)= W(n)%hTyy2(i,j,k) - d2* lam**2/lam2mu *ety *DyVy
!   hTxz(i,j,k)=0.0_SP
!   hTyz(i,j,k)=0.0_SP
!   hTzz(i,j,k)=0.0_SP
!end if
!end if
#ifdef AbsVzero
if (freenode .and. k==nk2) then
if (W(n)%isx) then
   W(n)%hTxx1(i,j,k)=W(n)%hTxx1(i,j,k) &
               + d1 * ( lam2mu*ztx*DzVx1 + lam*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTyy1(i,j,k)=W(n)%hTyy1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam2mu*zty*DzVy1 + lam*ztz*DzVz1 ) *b1
   W(n)%hTzz1(i,j,k)=W(n)%hTzz1(i,j,k) &
               + d1 * ( lam*ztx*DzVx1 + lam*zty*DzVy1 + lam2mu*ztz*DzVz1 ) *b1
   W(n)%hTyz1(i,j,k)=W(n)%hTyz1(i,j,k) &
               + d1 * miu * ( ztz*DzVy1 + zty*DzVz1 ) *b1
   W(n)%hTxz1(i,j,k)=W(n)%hTxz1(i,j,k) &
               + d1 * miu * ( ztz*DzVx1 + ztx*DzVz1 ) *b1
   W(n)%hTxy1(i,j,k)=W(n)%hTxy1(i,j,k) &
               + d1 * miu * ( zty*DzVx1 + ztx*DzVy1 ) *b1
end if
if (W(n)%isy) then
   W(n)%hTxx2(i,j,k)=W(n)%hTxx2(i,j,k) &
               + d2 * ( lam2mu*ztx*DzVx2 + lam*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTyy2(i,j,k)=W(n)%hTyy2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam2mu*zty*DzVy2 + lam*ztz*DzVz2 ) *b2
   W(n)%hTzz2(i,j,k)=W(n)%hTzz2(i,j,k) &
               + d2 * ( lam*ztx*DzVx2 + lam*zty*DzVy2 + lam2mu*ztz*DzVz2 ) *b2
   W(n)%hTyz2(i,j,k)=W(n)%hTyz2(i,j,k) &
               + d2 * miu * ( ztz*DzVy2 + zty*DzVz2 ) *b2
   W(n)%hTxz2(i,j,k)=W(n)%hTxz2(i,j,k) &
               + d2 * miu * ( ztz*DzVx2 + ztx*DzVz2 ) *b2
   W(n)%hTxy2(i,j,k)=W(n)%hTxy2(i,j,k) &
               + d2 * miu * ( zty*DzVx2 + ztx*DzVy2 ) *b2
end if
end if
#endif

end do
end do
end do
end do
end subroutine LxB_LyF_LzB

!*************************************************************************
!* use Traction Image method to calculate fd of stress compoents         *
!*************************************************************************

subroutine LxF_LyF_LzF_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n,n1
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

if (.not. freenode) return

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
if (k>=nk2-LenFD+1) then
   n1=nk2-k+1
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVx(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx1(i,j,k)/b1
   W(n)%hVx1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vx1(i,j,k)
end if
if (W(n)%isy) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vx2(i,j,k)
end if
if (W(n)%isz) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vx3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVy(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy1(i,j,k)/b1
   W(n)%hVy1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vy1(i,j,k)
end if
if (W(n)%isy) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy2(i,j,k)/b2
   W(n)%hVy2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vy2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVy3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vy3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVz(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz1(i,j,k)/b1
   W(n)%hVz1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vz1(i,j,k)
end if
if (W(n)%isy) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz2(i,j,k)/b2
   W(n)%hVz2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vz2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVz3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vz3(i,j,k)
end if

end do
end do
end if
end do
end do
end subroutine LxF_LyF_LzF_TIMG
subroutine LxB_LyB_LzB_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n,n1
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

if (.not. freenode) return

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
if (k>=nk2-LenFD+1) then
   n1=nk2-k+1
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVx(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx1(i,j,k)/b1
   W(n)%hVx1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vx1(i,j,k)
end if
if (W(n)%isy) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vx2(i,j,k)
end if
if (W(n)%isz) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vx3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVy(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy1(i,j,k)/b1
   W(n)%hVy1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vy1(i,j,k)
end if
if (W(n)%isy) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy2(i,j,k)/b2
   W(n)%hVy2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vy2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVy3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vy3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVz(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz1(i,j,k)/b1
   W(n)%hVz1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vz1(i,j,k)
end if
if (W(n)%isy) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz2(i,j,k)/b2
   W(n)%hVz2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vz2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVz3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vz3(i,j,k)
end if

end do
end do
end if
end do
end do
end subroutine LxB_LyB_LzB_TIMG
subroutine LxF_LyF_LzB_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n,n1
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

if (.not. freenode) return

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
if (k>=nk2-LenFD+1) then
   n1=nk2-k+1
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVx(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx1(i,j,k)/b1
   W(n)%hVx1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vx1(i,j,k)
end if
if (W(n)%isy) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vx2(i,j,k)
end if
if (W(n)%isz) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vx3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVy(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy1(i,j,k)/b1
   W(n)%hVy1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vy1(i,j,k)
end if
if (W(n)%isy) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy2(i,j,k)/b2
   W(n)%hVy2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vy2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVy3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vy3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVz(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz1(i,j,k)/b1
   W(n)%hVz1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vz1(i,j,k)
end if
if (W(n)%isy) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz2(i,j,k)/b2
   W(n)%hVz2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vz2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVz3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vz3(i,j,k)
end if

end do
end do
end if
end do
end do
end subroutine LxF_LyF_LzB_TIMG
subroutine LxB_LyB_LzF_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n,n1
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

if (.not. freenode) return

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
if (k>=nk2-LenFD+1) then
   n1=nk2-k+1
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVx(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx1(i,j,k)/b1
   W(n)%hVx1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vx1(i,j,k)
end if
if (W(n)%isy) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vx2(i,j,k)
end if
if (W(n)%isz) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vx3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVy(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy1(i,j,k)/b1
   W(n)%hVy1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vy1(i,j,k)
end if
if (W(n)%isy) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy2(i,j,k)/b2
   W(n)%hVy2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vy2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVy3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vy3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVz(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz1(i,j,k)/b1
   W(n)%hVz1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vz1(i,j,k)
end if
if (W(n)%isy) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz2(i,j,k)/b2
   W(n)%hVz2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vz2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVz3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vz3(i,j,k)
end if

end do
end do
end if
end do
end do
end subroutine LxB_LyB_LzF_TIMG
subroutine LxB_LyF_LzF_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n,n1
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

if (.not. freenode) return

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
if (k>=nk2-LenFD+1) then
   n1=nk2-k+1
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVx(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx1(i,j,k)/b1
   W(n)%hVx1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vx1(i,j,k)
end if
if (W(n)%isy) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vx2(i,j,k)
end if
if (W(n)%isz) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vx3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVy(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy1(i,j,k)/b1
   W(n)%hVy1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vy1(i,j,k)
end if
if (W(n)%isy) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy2(i,j,k)/b2
   W(n)%hVy2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vy2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVy3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vy3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVz(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz1(i,j,k)/b1
   W(n)%hVz1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vz1(i,j,k)
end if
if (W(n)%isy) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz2(i,j,k)/b2
   W(n)%hVz2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vz2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVz3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vz3(i,j,k)
end if

end do
end do
end if
end do
end do
end subroutine LxB_LyF_LzF_TIMG
subroutine LxF_LyB_LzB_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n,n1
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

if (.not. freenode) return

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
if (k>=nk2-LenFD+1) then
   n1=nk2-k+1
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVx(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx1(i,j,k)/b1
   W(n)%hVx1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vx1(i,j,k)
end if
if (W(n)%isy) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vx2(i,j,k)
end if
if (W(n)%isz) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vx3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVy(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy1(i,j,k)/b1
   W(n)%hVy1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vy1(i,j,k)
end if
if (W(n)%isy) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy2(i,j,k)/b2
   W(n)%hVy2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vy2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVy3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vy3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVz(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz1(i,j,k)/b1
   W(n)%hVz1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vz1(i,j,k)
end if
if (W(n)%isy) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz2(i,j,k)/b2
   W(n)%hVz2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vz2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVz3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vz3(i,j,k)
end if

end do
end do
end if
end do
end do
end subroutine LxF_LyB_LzB_TIMG
subroutine LxF_LyB_LzF_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n,n1
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

if (.not. freenode) return

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
if (k>=nk2-LenFD+1) then
   n1=nk2-k+1
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVx(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx1(i,j,k)/b1
   W(n)%hVx1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vx1(i,j,k)
end if
if (W(n)%isy) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vx2(i,j,k)
end if
if (W(n)%isz) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vx3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVy(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy1(i,j,k)/b1
   W(n)%hVy1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vy1(i,j,k)
end if
if (W(n)%isy) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy2(i,j,k)/b2
   W(n)%hVy2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vy2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVy3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vy3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVz(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz1(i,j,k)/b1
   W(n)%hVz1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vz1(i,j,k)
end if
if (W(n)%isy) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz2(i,j,k)/b2
   W(n)%hVz2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vz2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVz3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vz3(i,j,k)
end if

end do
end do
end if
end do
end do
end subroutine LxF_LyB_LzF_TIMG
subroutine LxB_LyF_LzB_TIMG
real(SP),dimension(-LenFD:LenFD) :: vecTx,vecTy,vecTz
real(SP) :: DxTx,DyTy,DzTz
real(SP) :: rrhojac
integer :: i,j,k,n,n1
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

if (.not. freenode) return

do n=1,num_blk
do k=W(n)%nk1,W(n)%nk2
if (k>=nk2-LenFD+1) then
   n1=nk2-k+1
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)

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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVx(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx1(i,j,k)/b1
   W(n)%hVx1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vx1(i,j,k)
end if
if (W(n)%isy) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx2(i,j,k)/b2
   W(n)%hVx2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vx2(i,j,k)
end if
if (W(n)%isz) then
   hVx(i,j,k)= hVx(i,j,k)-W(n)%Vx3(i,j,k)/b3
   W(n)%hVx3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vx3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVy(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy1(i,j,k)/b1
   W(n)%hVy1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vy1(i,j,k)
end if
if (W(n)%isy) then
   hVy(i,j,k)= hVy(i,j,k)-W(n)%Vy2(i,j,k)/b2
   W(n)%hVy2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vy2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVy3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vy3(i,j,k)
end if
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
   vecTz(n1:LenFD)=-vecTz(n1-2:n1-2 -(LenFD-n1):-1 )
   vecTz(n1-1)=0.0
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

   hVz(i,j,k)=( (DxTx )/b1 &
               +(DyTy )/b2 &
               +(DzTz )/b3 &
               )*rrhojac
if (W(n)%isx) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz1(i,j,k)/b1
   W(n)%hVz1(i,j,k)= d1*rrhojac*( DxTx ) - (d1+a1)*W(n)%Vz1(i,j,k)
end if
if (W(n)%isy) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz2(i,j,k)/b2
   W(n)%hVz2(i,j,k)= d2*rrhojac*( DyTy ) - (d2+a2)*W(n)%Vz2(i,j,k)
end if
if (W(n)%isz) then
   hVz(i,j,k)= hVz(i,j,k)-W(n)%Vz3(i,j,k)/b3
   W(n)%hVz3(i,j,k)= d3*rrhojac*( DzTz ) - (d3+a3)*W(n)%Vz3(i,j,k)
end if

end do
end do
end if
end do
end do
end subroutine LxB_LyF_LzB_TIMG

subroutine LxF_LyF_LzF_VHOC
integer :: i,j,k,n,nblk
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef AbsVzero
   return
#endif

if (.not. freenode) return

do nblk=1,num_blk
   if (W(nblk)%nk2/=nk2) cycle
loop_eta: do j=W(nblk)%nj1,W(nblk)%nj2
loop_xi:  do i=W(nblk)%ni1,W(nblk)%ni2
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
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)
   hTxx(i,j,k)= ( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) )/b3
   hTxy(i,j,k)= miu *(                           &
                ( xiy*DxVx(n) + xix*DxVy(n) )/b1 &
               +( ety*DyVx(n) + etx*DyVy(n) )/b2 &
               +( zty*DzVx(n) + ztx*DzVy(n) )/b3 &
               )
   hTxz(i,j,k)= miu *(                           &
                ( xiz*DxVx(n) + xix*DxVz(n) )/b1 &
               +( etz*DyVx(n) + etx*DyVz(n) )/b2 &
               +( ztz*DzVx(n) + ztx*DzVz(n) )/b3 &
               )
   hTyz(i,j,k)= miu *(                           &
                ( xiz*DxVy(n) + xiy*DxVz(n) )/b1 &
               +( etz*DyVy(n) + ety*DyVz(n) )/b2 &
               +( ztz*DzVy(n) + zty*DzVz(n) )/b3 &
               )
if (W(nblk)%isx) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz1(i,j,k)/b1
   W(nblk)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txx1(i,j,k)  
   W(nblk)%hTyy1(i,j,k)= d1*( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyy1(i,j,k)  
   W(nblk)%hTzz1(i,j,k)= d1*( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tzz1(i,j,k) 
   W(nblk)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx(n) + xix*DxVy(n) ) &
               - (d1+a1)*W(nblk)%Txy1(i,j,k)  
   W(nblk)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx(n) + xix*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txz1(i,j,k)  
   W(nblk)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy(n) + xiy*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyz1(i,j,k)  
end if
if (W(nblk)%isy) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz2(i,j,k)/b2
   W(nblk)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txx2(i,j,k)  
   W(nblk)%hTyy2(i,j,k)= d2*( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyy2(i,j,k)  
   W(nblk)%hTzz2(i,j,k)= d2*( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tzz2(i,j,k) 
   W(nblk)%hTxy2(i,j,k)= d2*miu*( ety*DyVx(n) + etx*DyVy(n) ) &
               - (d2+a2)*W(nblk)%Txy2(i,j,k)  
   W(nblk)%hTxz2(i,j,k)= d2*miu*( etz*DyVx(n) + etx*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txz2(i,j,k)  
   W(nblk)%hTyz2(i,j,k)= d2*miu*( etz*DyVy(n) + ety*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyz2(i,j,k)  
end if
if (W(nblk)%isz) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz3(i,j,k)/b3
   W(nblk)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txx3(i,j,k)  
   W(nblk)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyy3(i,j,k)  
   W(nblk)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tzz3(i,j,k) 
   W(nblk)%hTxy3(i,j,k)= d3*miu*( zty*DzVx(n) + ztx*DzVy(n) ) &
               - (d3+a3)*W(nblk)%Txy3(i,j,k)  
   W(nblk)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx(n) + ztx*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txz3(i,j,k) 
   W(nblk)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy(n) + zty*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyz3(i,j,k) 
end if
end do
end do loop_xi
end do loop_eta
end do
end subroutine LxF_LyF_LzF_VHOC
subroutine LxB_LyB_LzB_VHOC
integer :: i,j,k,n,nblk
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef AbsVzero
   return
#endif

if (.not. freenode) return

do nblk=1,num_blk
   if (W(nblk)%nk2/=nk2) cycle
loop_eta: do j=W(nblk)%nj1,W(nblk)%nj2
loop_xi:  do i=W(nblk)%ni1,W(nblk)%ni2
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
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)
   hTxx(i,j,k)= ( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) )/b3
   hTxy(i,j,k)= miu *(                           &
                ( xiy*DxVx(n) + xix*DxVy(n) )/b1 &
               +( ety*DyVx(n) + etx*DyVy(n) )/b2 &
               +( zty*DzVx(n) + ztx*DzVy(n) )/b3 &
               )
   hTxz(i,j,k)= miu *(                           &
                ( xiz*DxVx(n) + xix*DxVz(n) )/b1 &
               +( etz*DyVx(n) + etx*DyVz(n) )/b2 &
               +( ztz*DzVx(n) + ztx*DzVz(n) )/b3 &
               )
   hTyz(i,j,k)= miu *(                           &
                ( xiz*DxVy(n) + xiy*DxVz(n) )/b1 &
               +( etz*DyVy(n) + ety*DyVz(n) )/b2 &
               +( ztz*DzVy(n) + zty*DzVz(n) )/b3 &
               )
if (W(nblk)%isx) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz1(i,j,k)/b1
   W(nblk)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txx1(i,j,k)  
   W(nblk)%hTyy1(i,j,k)= d1*( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyy1(i,j,k)  
   W(nblk)%hTzz1(i,j,k)= d1*( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tzz1(i,j,k) 
   W(nblk)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx(n) + xix*DxVy(n) ) &
               - (d1+a1)*W(nblk)%Txy1(i,j,k)  
   W(nblk)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx(n) + xix*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txz1(i,j,k)  
   W(nblk)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy(n) + xiy*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyz1(i,j,k)  
end if
if (W(nblk)%isy) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz2(i,j,k)/b2
   W(nblk)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txx2(i,j,k)  
   W(nblk)%hTyy2(i,j,k)= d2*( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyy2(i,j,k)  
   W(nblk)%hTzz2(i,j,k)= d2*( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tzz2(i,j,k) 
   W(nblk)%hTxy2(i,j,k)= d2*miu*( ety*DyVx(n) + etx*DyVy(n) ) &
               - (d2+a2)*W(nblk)%Txy2(i,j,k)  
   W(nblk)%hTxz2(i,j,k)= d2*miu*( etz*DyVx(n) + etx*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txz2(i,j,k)  
   W(nblk)%hTyz2(i,j,k)= d2*miu*( etz*DyVy(n) + ety*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyz2(i,j,k)  
end if
if (W(nblk)%isz) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz3(i,j,k)/b3
   W(nblk)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txx3(i,j,k)  
   W(nblk)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyy3(i,j,k)  
   W(nblk)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tzz3(i,j,k) 
   W(nblk)%hTxy3(i,j,k)= d3*miu*( zty*DzVx(n) + ztx*DzVy(n) ) &
               - (d3+a3)*W(nblk)%Txy3(i,j,k)  
   W(nblk)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx(n) + ztx*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txz3(i,j,k) 
   W(nblk)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy(n) + zty*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyz3(i,j,k) 
end if
end do
end do loop_xi
end do loop_eta
end do
end subroutine LxB_LyB_LzB_VHOC

subroutine LxF_LyF_LzB_VHOC
integer :: i,j,k,n,nblk
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef AbsVzero
   return
#endif

if (.not. freenode) return

do nblk=1,num_blk
   if (W(nblk)%nk2/=nk2) cycle
loop_eta: do j=W(nblk)%nj1,W(nblk)%nj2
loop_xi:  do i=W(nblk)%ni1,W(nblk)%ni2
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
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)
   hTxx(i,j,k)= ( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) )/b3
   hTxy(i,j,k)= miu *(                           &
                ( xiy*DxVx(n) + xix*DxVy(n) )/b1 &
               +( ety*DyVx(n) + etx*DyVy(n) )/b2 &
               +( zty*DzVx(n) + ztx*DzVy(n) )/b3 &
               )
   hTxz(i,j,k)= miu *(                           &
                ( xiz*DxVx(n) + xix*DxVz(n) )/b1 &
               +( etz*DyVx(n) + etx*DyVz(n) )/b2 &
               +( ztz*DzVx(n) + ztx*DzVz(n) )/b3 &
               )
   hTyz(i,j,k)= miu *(                           &
                ( xiz*DxVy(n) + xiy*DxVz(n) )/b1 &
               +( etz*DyVy(n) + ety*DyVz(n) )/b2 &
               +( ztz*DzVy(n) + zty*DzVz(n) )/b3 &
               )
if (W(nblk)%isx) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz1(i,j,k)/b1
   W(nblk)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txx1(i,j,k)  
   W(nblk)%hTyy1(i,j,k)= d1*( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyy1(i,j,k)  
   W(nblk)%hTzz1(i,j,k)= d1*( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tzz1(i,j,k) 
   W(nblk)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx(n) + xix*DxVy(n) ) &
               - (d1+a1)*W(nblk)%Txy1(i,j,k)  
   W(nblk)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx(n) + xix*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txz1(i,j,k)  
   W(nblk)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy(n) + xiy*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyz1(i,j,k)  
end if
if (W(nblk)%isy) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz2(i,j,k)/b2
   W(nblk)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txx2(i,j,k)  
   W(nblk)%hTyy2(i,j,k)= d2*( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyy2(i,j,k)  
   W(nblk)%hTzz2(i,j,k)= d2*( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tzz2(i,j,k) 
   W(nblk)%hTxy2(i,j,k)= d2*miu*( ety*DyVx(n) + etx*DyVy(n) ) &
               - (d2+a2)*W(nblk)%Txy2(i,j,k)  
   W(nblk)%hTxz2(i,j,k)= d2*miu*( etz*DyVx(n) + etx*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txz2(i,j,k)  
   W(nblk)%hTyz2(i,j,k)= d2*miu*( etz*DyVy(n) + ety*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyz2(i,j,k)  
end if
if (W(nblk)%isz) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz3(i,j,k)/b3
   W(nblk)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txx3(i,j,k)  
   W(nblk)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyy3(i,j,k)  
   W(nblk)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tzz3(i,j,k) 
   W(nblk)%hTxy3(i,j,k)= d3*miu*( zty*DzVx(n) + ztx*DzVy(n) ) &
               - (d3+a3)*W(nblk)%Txy3(i,j,k)  
   W(nblk)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx(n) + ztx*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txz3(i,j,k) 
   W(nblk)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy(n) + zty*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyz3(i,j,k) 
end if
end do
end do loop_xi
end do loop_eta
end do
end subroutine LxF_LyF_LzB_VHOC
subroutine LxB_LyB_LzF_VHOC
integer :: i,j,k,n,nblk
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3
#ifdef AbsVzero
   return
#endif

if (.not. freenode) return

do nblk=1,num_blk
   if (W(nblk)%nk2/=nk2) cycle
loop_eta: do j=W(nblk)%nj1,W(nblk)%nj2
loop_xi:  do i=W(nblk)%ni1,W(nblk)%ni2
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
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)
   hTxx(i,j,k)= ( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) )/b3
   hTxy(i,j,k)= miu *(                           &
                ( xiy*DxVx(n) + xix*DxVy(n) )/b1 &
               +( ety*DyVx(n) + etx*DyVy(n) )/b2 &
               +( zty*DzVx(n) + ztx*DzVy(n) )/b3 &
               )
   hTxz(i,j,k)= miu *(                           &
                ( xiz*DxVx(n) + xix*DxVz(n) )/b1 &
               +( etz*DyVx(n) + etx*DyVz(n) )/b2 &
               +( ztz*DzVx(n) + ztx*DzVz(n) )/b3 &
               )
   hTyz(i,j,k)= miu *(                           &
                ( xiz*DxVy(n) + xiy*DxVz(n) )/b1 &
               +( etz*DyVy(n) + ety*DyVz(n) )/b2 &
               +( ztz*DzVy(n) + zty*DzVz(n) )/b3 &
               )
if (W(nblk)%isx) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz1(i,j,k)/b1
   W(nblk)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txx1(i,j,k)  
   W(nblk)%hTyy1(i,j,k)= d1*( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyy1(i,j,k)  
   W(nblk)%hTzz1(i,j,k)= d1*( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tzz1(i,j,k) 
   W(nblk)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx(n) + xix*DxVy(n) ) &
               - (d1+a1)*W(nblk)%Txy1(i,j,k)  
   W(nblk)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx(n) + xix*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txz1(i,j,k)  
   W(nblk)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy(n) + xiy*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyz1(i,j,k)  
end if
if (W(nblk)%isy) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz2(i,j,k)/b2
   W(nblk)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txx2(i,j,k)  
   W(nblk)%hTyy2(i,j,k)= d2*( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyy2(i,j,k)  
   W(nblk)%hTzz2(i,j,k)= d2*( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tzz2(i,j,k) 
   W(nblk)%hTxy2(i,j,k)= d2*miu*( ety*DyVx(n) + etx*DyVy(n) ) &
               - (d2+a2)*W(nblk)%Txy2(i,j,k)  
   W(nblk)%hTxz2(i,j,k)= d2*miu*( etz*DyVx(n) + etx*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txz2(i,j,k)  
   W(nblk)%hTyz2(i,j,k)= d2*miu*( etz*DyVy(n) + ety*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyz2(i,j,k)  
end if
if (W(nblk)%isz) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz3(i,j,k)/b3
   W(nblk)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txx3(i,j,k)  
   W(nblk)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyy3(i,j,k)  
   W(nblk)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tzz3(i,j,k) 
   W(nblk)%hTxy3(i,j,k)= d3*miu*( zty*DzVx(n) + ztx*DzVy(n) ) &
               - (d3+a3)*W(nblk)%Txy3(i,j,k)  
   W(nblk)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx(n) + ztx*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txz3(i,j,k) 
   W(nblk)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy(n) + zty*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyz3(i,j,k) 
end if
end do
end do loop_xi
end do loop_eta
end do
end subroutine LxB_LyB_LzF_VHOC

subroutine LxB_LyF_LzF_VHOC
integer :: i,j,k,n,nblk
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef AbsVzero
   return
#endif

if (.not. freenode) return

do nblk=1,num_blk
   if (W(nblk)%nk2/=nk2) cycle
loop_eta: do j=W(nblk)%nj1,W(nblk)%nj2
loop_xi:  do i=W(nblk)%ni1,W(nblk)%ni2
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
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)
   hTxx(i,j,k)= ( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) )/b3
   hTxy(i,j,k)= miu *(                           &
                ( xiy*DxVx(n) + xix*DxVy(n) )/b1 &
               +( ety*DyVx(n) + etx*DyVy(n) )/b2 &
               +( zty*DzVx(n) + ztx*DzVy(n) )/b3 &
               )
   hTxz(i,j,k)= miu *(                           &
                ( xiz*DxVx(n) + xix*DxVz(n) )/b1 &
               +( etz*DyVx(n) + etx*DyVz(n) )/b2 &
               +( ztz*DzVx(n) + ztx*DzVz(n) )/b3 &
               )
   hTyz(i,j,k)= miu *(                           &
                ( xiz*DxVy(n) + xiy*DxVz(n) )/b1 &
               +( etz*DyVy(n) + ety*DyVz(n) )/b2 &
               +( ztz*DzVy(n) + zty*DzVz(n) )/b3 &
               )
if (W(nblk)%isx) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz1(i,j,k)/b1
   W(nblk)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txx1(i,j,k)  
   W(nblk)%hTyy1(i,j,k)= d1*( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyy1(i,j,k)  
   W(nblk)%hTzz1(i,j,k)= d1*( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tzz1(i,j,k) 
   W(nblk)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx(n) + xix*DxVy(n) ) &
               - (d1+a1)*W(nblk)%Txy1(i,j,k)  
   W(nblk)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx(n) + xix*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txz1(i,j,k)  
   W(nblk)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy(n) + xiy*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyz1(i,j,k)  
end if
if (W(nblk)%isy) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz2(i,j,k)/b2
   W(nblk)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txx2(i,j,k)  
   W(nblk)%hTyy2(i,j,k)= d2*( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyy2(i,j,k)  
   W(nblk)%hTzz2(i,j,k)= d2*( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tzz2(i,j,k) 
   W(nblk)%hTxy2(i,j,k)= d2*miu*( ety*DyVx(n) + etx*DyVy(n) ) &
               - (d2+a2)*W(nblk)%Txy2(i,j,k)  
   W(nblk)%hTxz2(i,j,k)= d2*miu*( etz*DyVx(n) + etx*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txz2(i,j,k)  
   W(nblk)%hTyz2(i,j,k)= d2*miu*( etz*DyVy(n) + ety*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyz2(i,j,k)  
end if
if (W(nblk)%isz) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz3(i,j,k)/b3
   W(nblk)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txx3(i,j,k)  
   W(nblk)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyy3(i,j,k)  
   W(nblk)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tzz3(i,j,k) 
   W(nblk)%hTxy3(i,j,k)= d3*miu*( zty*DzVx(n) + ztx*DzVy(n) ) &
               - (d3+a3)*W(nblk)%Txy3(i,j,k)  
   W(nblk)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx(n) + ztx*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txz3(i,j,k) 
   W(nblk)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy(n) + zty*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyz3(i,j,k) 
end if
end do
end do loop_xi
end do loop_eta
end do
end subroutine LxB_LyF_LzF_VHOC
subroutine LxF_LyB_LzB_VHOC
integer :: i,j,k,n,nblk
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef AbsVzero
   return
#endif

if (.not. freenode) return

do nblk=1,num_blk
   if (W(nblk)%nk2/=nk2) cycle
loop_eta: do j=W(nblk)%nj1,W(nblk)%nj2
loop_xi:  do i=W(nblk)%ni1,W(nblk)%ni2
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
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)
   hTxx(i,j,k)= ( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) )/b3
   hTxy(i,j,k)= miu *(                           &
                ( xiy*DxVx(n) + xix*DxVy(n) )/b1 &
               +( ety*DyVx(n) + etx*DyVy(n) )/b2 &
               +( zty*DzVx(n) + ztx*DzVy(n) )/b3 &
               )
   hTxz(i,j,k)= miu *(                           &
                ( xiz*DxVx(n) + xix*DxVz(n) )/b1 &
               +( etz*DyVx(n) + etx*DyVz(n) )/b2 &
               +( ztz*DzVx(n) + ztx*DzVz(n) )/b3 &
               )
   hTyz(i,j,k)= miu *(                           &
                ( xiz*DxVy(n) + xiy*DxVz(n) )/b1 &
               +( etz*DyVy(n) + ety*DyVz(n) )/b2 &
               +( ztz*DzVy(n) + zty*DzVz(n) )/b3 &
               )
if (W(nblk)%isx) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz1(i,j,k)/b1
   W(nblk)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txx1(i,j,k)  
   W(nblk)%hTyy1(i,j,k)= d1*( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyy1(i,j,k)  
   W(nblk)%hTzz1(i,j,k)= d1*( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tzz1(i,j,k) 
   W(nblk)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx(n) + xix*DxVy(n) ) &
               - (d1+a1)*W(nblk)%Txy1(i,j,k)  
   W(nblk)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx(n) + xix*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txz1(i,j,k)  
   W(nblk)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy(n) + xiy*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyz1(i,j,k)  
end if
if (W(nblk)%isy) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz2(i,j,k)/b2
   W(nblk)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txx2(i,j,k)  
   W(nblk)%hTyy2(i,j,k)= d2*( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyy2(i,j,k)  
   W(nblk)%hTzz2(i,j,k)= d2*( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tzz2(i,j,k) 
   W(nblk)%hTxy2(i,j,k)= d2*miu*( ety*DyVx(n) + etx*DyVy(n) ) &
               - (d2+a2)*W(nblk)%Txy2(i,j,k)  
   W(nblk)%hTxz2(i,j,k)= d2*miu*( etz*DyVx(n) + etx*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txz2(i,j,k)  
   W(nblk)%hTyz2(i,j,k)= d2*miu*( etz*DyVy(n) + ety*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyz2(i,j,k)  
end if
if (W(nblk)%isz) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz3(i,j,k)/b3
   W(nblk)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txx3(i,j,k)  
   W(nblk)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyy3(i,j,k)  
   W(nblk)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tzz3(i,j,k) 
   W(nblk)%hTxy3(i,j,k)= d3*miu*( zty*DzVx(n) + ztx*DzVy(n) ) &
               - (d3+a3)*W(nblk)%Txy3(i,j,k)  
   W(nblk)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx(n) + ztx*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txz3(i,j,k) 
   W(nblk)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy(n) + zty*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyz3(i,j,k) 
end if
end do
end do loop_xi
end do loop_eta
end do
end subroutine LxF_LyB_LzB_VHOC

subroutine LxF_LyB_LzF_VHOC
integer :: i,j,k,n,nblk
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef AbsVzero
   return
#endif

if (.not. freenode) return

do nblk=1,num_blk
   if (W(nblk)%nk2/=nk2) cycle
loop_eta: do j=W(nblk)%nj1,W(nblk)%nj2
loop_xi:  do i=W(nblk)%ni1,W(nblk)%ni2
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
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)
   hTxx(i,j,k)= ( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) )/b3
   hTxy(i,j,k)= miu *(                           &
                ( xiy*DxVx(n) + xix*DxVy(n) )/b1 &
               +( ety*DyVx(n) + etx*DyVy(n) )/b2 &
               +( zty*DzVx(n) + ztx*DzVy(n) )/b3 &
               )
   hTxz(i,j,k)= miu *(                           &
                ( xiz*DxVx(n) + xix*DxVz(n) )/b1 &
               +( etz*DyVx(n) + etx*DyVz(n) )/b2 &
               +( ztz*DzVx(n) + ztx*DzVz(n) )/b3 &
               )
   hTyz(i,j,k)= miu *(                           &
                ( xiz*DxVy(n) + xiy*DxVz(n) )/b1 &
               +( etz*DyVy(n) + ety*DyVz(n) )/b2 &
               +( ztz*DzVy(n) + zty*DzVz(n) )/b3 &
               )
if (W(nblk)%isx) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz1(i,j,k)/b1
   W(nblk)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txx1(i,j,k)  
   W(nblk)%hTyy1(i,j,k)= d1*( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyy1(i,j,k)  
   W(nblk)%hTzz1(i,j,k)= d1*( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tzz1(i,j,k) 
   W(nblk)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx(n) + xix*DxVy(n) ) &
               - (d1+a1)*W(nblk)%Txy1(i,j,k)  
   W(nblk)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx(n) + xix*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txz1(i,j,k)  
   W(nblk)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy(n) + xiy*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyz1(i,j,k)  
end if
if (W(nblk)%isy) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz2(i,j,k)/b2
   W(nblk)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txx2(i,j,k)  
   W(nblk)%hTyy2(i,j,k)= d2*( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyy2(i,j,k)  
   W(nblk)%hTzz2(i,j,k)= d2*( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tzz2(i,j,k) 
   W(nblk)%hTxy2(i,j,k)= d2*miu*( ety*DyVx(n) + etx*DyVy(n) ) &
               - (d2+a2)*W(nblk)%Txy2(i,j,k)  
   W(nblk)%hTxz2(i,j,k)= d2*miu*( etz*DyVx(n) + etx*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txz2(i,j,k)  
   W(nblk)%hTyz2(i,j,k)= d2*miu*( etz*DyVy(n) + ety*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyz2(i,j,k)  
end if
if (W(nblk)%isz) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz3(i,j,k)/b3
   W(nblk)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txx3(i,j,k)  
   W(nblk)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyy3(i,j,k)  
   W(nblk)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tzz3(i,j,k) 
   W(nblk)%hTxy3(i,j,k)= d3*miu*( zty*DzVx(n) + ztx*DzVy(n) ) &
               - (d3+a3)*W(nblk)%Txy3(i,j,k)  
   W(nblk)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx(n) + ztx*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txz3(i,j,k) 
   W(nblk)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy(n) + zty*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyz3(i,j,k) 
end if
end do
end do loop_xi
end do loop_eta
end do
end subroutine LxF_LyB_LzF_VHOC
subroutine LxB_LyF_LzB_VHOC
integer :: i,j,k,n,nblk
real(SP) :: lam,miu,lam2mu,rrho
real(SP) :: xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz
real(SP) :: rhs_Dz,lhs_Dz
real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                             DyVx,DyVy,DyVz, &
                             DzVx,DzVy,DzVz
real(SP) :: d1,d2,d3
real(SP) :: b1,b2,b3
real(SP) :: a1,a2,a3

#ifdef AbsVzero
   return
#endif

if (.not. freenode) return

do nblk=1,num_blk
   if (W(nblk)%nk2/=nk2) cycle
loop_eta: do j=W(nblk)%nj1,W(nblk)%nj2
loop_xi:  do i=W(nblk)%ni1,W(nblk)%ni2
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
   d1=Dx(i); d2=Dy(j); d3=Dz(k)
   b1=Bx(i); b2=By(j); b3=Bz(k)
   a1=Ax(i); a2=Ay(j); a3=Az(k)
   hTxx(i,j,k)= ( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTyy(i,j,k)= ( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) )/b3
   hTzz(i,j,k)= ( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) )/b1 &
               +( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) )/b2 &
               +( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) )/b3
   hTxy(i,j,k)= miu *(                           &
                ( xiy*DxVx(n) + xix*DxVy(n) )/b1 &
               +( ety*DyVx(n) + etx*DyVy(n) )/b2 &
               +( zty*DzVx(n) + ztx*DzVy(n) )/b3 &
               )
   hTxz(i,j,k)= miu *(                           &
                ( xiz*DxVx(n) + xix*DxVz(n) )/b1 &
               +( etz*DyVx(n) + etx*DyVz(n) )/b2 &
               +( ztz*DzVx(n) + ztx*DzVz(n) )/b3 &
               )
   hTyz(i,j,k)= miu *(                           &
                ( xiz*DxVy(n) + xiy*DxVz(n) )/b1 &
               +( etz*DyVy(n) + ety*DyVz(n) )/b2 &
               +( ztz*DzVy(n) + zty*DzVz(n) )/b3 &
               )
if (W(nblk)%isx) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx1(i,j,k)/b1
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy1(i,j,k)/b1
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz1(i,j,k)/b1
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy1(i,j,k)/b1
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz1(i,j,k)/b1
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz1(i,j,k)/b1
   W(nblk)%hTxx1(i,j,k)= d1*( lam2mu*xix*DxVx(n) + lam*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txx1(i,j,k)  
   W(nblk)%hTyy1(i,j,k)= d1*( lam*xix*DxVx(n) + lam2mu*xiy*DxVy(n) + lam*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyy1(i,j,k)  
   W(nblk)%hTzz1(i,j,k)= d1*( lam*xix*DxVx(n) + lam*xiy*DxVy(n) + lam2mu*xiz*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tzz1(i,j,k) 
   W(nblk)%hTxy1(i,j,k)= d1*miu*( xiy*DxVx(n) + xix*DxVy(n) ) &
               - (d1+a1)*W(nblk)%Txy1(i,j,k)  
   W(nblk)%hTxz1(i,j,k)= d1*miu*( xiz*DxVx(n) + xix*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Txz1(i,j,k)  
   W(nblk)%hTyz1(i,j,k)= d1*miu*( xiz*DxVy(n) + xiy*DxVz(n) ) &
               - (d1+a1)*W(nblk)%Tyz1(i,j,k)  
end if
if (W(nblk)%isy) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx2(i,j,k)/b2
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy2(i,j,k)/b2
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz2(i,j,k)/b2
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy2(i,j,k)/b2
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz2(i,j,k)/b2
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz2(i,j,k)/b2
   W(nblk)%hTxx2(i,j,k)= d2*( lam2mu*etx*DyVx(n) + lam*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txx2(i,j,k)  
   W(nblk)%hTyy2(i,j,k)= d2*( lam*etx*DyVx(n) + lam2mu*ety*DyVy(n) + lam*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyy2(i,j,k)  
   W(nblk)%hTzz2(i,j,k)= d2*( lam*etx*DyVx(n) + lam*ety*DyVy(n) + lam2mu*etz*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tzz2(i,j,k) 
   W(nblk)%hTxy2(i,j,k)= d2*miu*( ety*DyVx(n) + etx*DyVy(n) ) &
               - (d2+a2)*W(nblk)%Txy2(i,j,k)  
   W(nblk)%hTxz2(i,j,k)= d2*miu*( etz*DyVx(n) + etx*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Txz2(i,j,k)  
   W(nblk)%hTyz2(i,j,k)= d2*miu*( etz*DyVy(n) + ety*DyVz(n) ) &
               - (d2+a2)*W(nblk)%Tyz2(i,j,k)  
end if
if (W(nblk)%isz) then
   hTxx(i,j,k)=hTxx(i,j,k)-W(nblk)%Txx3(i,j,k)/b3
   hTyy(i,j,k)=hTyy(i,j,k)-W(nblk)%Tyy3(i,j,k)/b3
   hTzz(i,j,k)=hTzz(i,j,k)-W(nblk)%Tzz3(i,j,k)/b3
   hTxy(i,j,k)=hTxy(i,j,k)-W(nblk)%Txy3(i,j,k)/b3
   hTxz(i,j,k)=hTxz(i,j,k)-W(nblk)%Txz3(i,j,k)/b3
   hTyz(i,j,k)=hTyz(i,j,k)-W(nblk)%Tyz3(i,j,k)/b3
   W(nblk)%hTxx3(i,j,k)= d3*( lam2mu*ztx*DzVx(n) + lam*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txx3(i,j,k)  
   W(nblk)%hTyy3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam2mu*zty*DzVy(n) + lam*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyy3(i,j,k)  
   W(nblk)%hTzz3(i,j,k)= d3*( lam*ztx*DzVx(n) + lam*zty*DzVy(n) + lam2mu*ztz*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tzz3(i,j,k) 
   W(nblk)%hTxy3(i,j,k)= d3*miu*( zty*DzVx(n) + ztx*DzVy(n) ) &
               - (d3+a3)*W(nblk)%Txy3(i,j,k)  
   W(nblk)%hTxz3(i,j,k)= d3*miu*( ztz*DzVx(n) + ztx*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Txz3(i,j,k) 
   W(nblk)%hTyz3(i,j,k)= d3*miu*( ztz*DzVy(n) + zty*DzVz(n) ) &
               - (d3+a3)*W(nblk)%Tyz3(i,j,k) 
end if
end do
end do loop_xi
end do loop_eta
end do
end subroutine LxB_LyF_LzB_VHOC

subroutine abs_charac
integer i,j,k,n
real(SP) v1,v2,v3,t11,t22,t33,t12,t13,t23
real(SP) lam,miu,lam2mu,rrho,f1,f2,fct

return

if (absnode(1,1)) then
do n=1,num_blk
if (W(n)%ni1==ni1 .and. W(n)%isabs) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
   i=ni1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=rho(i,j,k)
   f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu/2.0
   v1 = Vx (i,j,k)
   v2 = Vy (i,j,k)
   v3 = Vz (i,j,k)
   t11= Txx(i,j,k)
   t22= Tyy(i,j,k)
   t33= Tzz(i,j,k)
   t12= Txy(i,j,k)
   t13= Txz(i,j,k)
   t23= Tyz(i,j,k)
   Vx (i,j,k)=0.5*(v1+t11/f1)
   Vy (i,j,k)=0.5*(v2+t12/f2)
   Vz (i,j,k)=0.5*(v3+t13/f2)
   Txx(i,j,k)=0.5*(t11+f1*v1)
   Tyy(i,j,k)=t22-fct*(t11-f1*v1)
   Tzz(i,j,k)=t33-fct*(t11-f1*v1)
   Txy(i,j,k)=0.5*(t12+f2*v2)
   Txz(i,j,k)=0.5*(t13+f2*v3)
end do
end do
end if
end do
end if

if (absnode(1,2)) then
do n=1,num_blk
if (W(n)%ni2==ni2 .and. W(n)%isabs) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
   i=ni2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=rho(i,j,k)
   f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu/2.0
   v1 = Vx (i,j,k)
   v2 = Vy (i,j,k)
   v3 = Vz (i,j,k)
   t11= Txx(i,j,k)
   t22= Tyy(i,j,k)
   t33= Tzz(i,j,k)
   t12= Txy(i,j,k)
   t13= Txz(i,j,k)
   t23= Tyz(i,j,k)
   Vx (i,j,k)=0.5*(v1-t11/f1)
   Vy (i,j,k)=0.5*(v2-t12/f2)
   Vz (i,j,k)=0.5*(v3-t13/f2)
   Txx(i,j,k)=0.5*(t11-f1*v1)
   Tyy(i,j,k)=t22-fct*(t11+f1*v1)
   Tzz(i,j,k)=t33-fct*(t11+f1*v1)
   Txy(i,j,k)=0.5*(t12-f2*v2)
   Txz(i,j,k)=0.5*(t13-f2*v3)
end do
end do
end if
end do
end if

if (absnode(2,1)) then
do n=1,num_blk
if (W(n)%nj1==nj1 .and. W(n)%isabs) then
do k=W(n)%nk1,W(n)%nk2
do i=W(n)%ni1,W(n)%ni2
   j=nj1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=rho(i,j,k)
   f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu/2.0
   v1 = Vx (i,j,k)
   v2 = Vy (i,j,k)
   v3 = Vz (i,j,k)
   t11= Txx(i,j,k)
   t22= Tyy(i,j,k)
   t33= Tzz(i,j,k)
   t12= Txy(i,j,k)
   t13= Txz(i,j,k)
   t23= Tyz(i,j,k)
   Vx (i,j,k)=0.5*(v1+t11/f2)
   Vy (i,j,k)=0.5*(v2+t12/f1)
   Vz (i,j,k)=0.5*(v3+t13/f2)
   Txx(i,j,k)=t11-fct*(t22-f1*v2)
   Tyy(i,j,k)=0.5*(t22+f1*v2)
   Tzz(i,j,k)=t33-fct*(t22-f1*v2)
   Txy(i,j,k)=0.5*(t12+f2*v1)
   Tyz(i,j,k)=0.5*(t23+f2*v3)
end do
end do
end if
end do
end if
if (absnode(2,2)) then
do n=1,num_blk
if (W(n)%nj2==nj2 .and. W(n)%isabs) then
do k=W(n)%nk1,W(n)%nk2
do i=W(n)%ni1,W(n)%ni2
   j=nj2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=rho(i,j,k)
   f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu/2.0
   v1 = Vx (i,j,k)
   v2 = Vy (i,j,k)
   v3 = Vz (i,j,k)
   t11= Txx(i,j,k)
   t22= Tyy(i,j,k)
   t33= Tzz(i,j,k)
   t12= Txy(i,j,k)
   t13= Txz(i,j,k)
   t23= Tyz(i,j,k)
   Vx (i,j,k)=0.5*(v1-t11/f2)
   Vy (i,j,k)=0.5*(v2-t12/f1)
   Vz (i,j,k)=0.5*(v3-t13/f2)
   Txx(i,j,k)=t11-fct*(t22+f1*v2)
   Tyy(i,j,k)=0.5*(t22-f1*v2)
   Tzz(i,j,k)=t33-fct*(t22+f1*v2)
   Txy(i,j,k)=0.5*(t12-f2*v1)
   Tyz(i,j,k)=0.5*(t23-f2*v3)
end do
end do
end if
end do
end if

if (absnode(3,1)) then
do n=1,num_blk
if (W(n)%nk1==nk1 .and. W(n)%isabs) then
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   k=nk1
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=rho(i,j,k)
   f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu/2.0
   v1 = Vx (i,j,k)
   v2 = Vy (i,j,k)
   v3 = Vz (i,j,k)
   t11= Txx(i,j,k)
   t22= Tyy(i,j,k)
   t33= Tzz(i,j,k)
   t12= Txy(i,j,k)
   t13= Txz(i,j,k)
   t23= Tyz(i,j,k)
   Vx (i,j,k)=0.5*(v1+t13/f2)
   Vy (i,j,k)=0.5*(v2+t23/f2)
   Vz (i,j,k)=0.5*(v3+t33/f1)
   Txx(i,j,k)=t11-fct*(t33-f1*v3)
   Tyy(i,j,k)=t22-fct*(t33-f1*v3)
   Tzz(i,j,k)=0.5*(t33+f1*v3)
   Txz(i,j,k)=0.5*(t13+f2*v1)
   Tyz(i,j,k)=0.5*(t23+f2*v2)
end do
end do
end if
end do
end if
if (absnode(3,2)) then
do n=1,num_blk
if (W(n)%nk2==nk2 .and. W(n)%isabs) then
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   k=nk2
   lam=lambda(i,j,k);miu=mu(i,j,k);lam2mu=lam+2.0*miu;rrho=rho(i,j,k)
   f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu/2.0
   v1 = Vx (i,j,k)
   v2 = Vy (i,j,k)
   v3 = Vz (i,j,k)
   t11= Txx(i,j,k)
   t22= Tyy(i,j,k)
   t33= Tzz(i,j,k)
   t12= Txy(i,j,k)
   t13= Txz(i,j,k)
   t23= Tyz(i,j,k)
   Vx (i,j,k)=0.5*(v1-t13/f2)
   Vy (i,j,k)=0.5*(v2-t23/f2)
   Vz (i,j,k)=0.5*(v3-t33/f1)
   Txx(i,j,k)=t11-fct*(t33+f1*v3)
   Tyy(i,j,k)=t22-fct*(t33+f1*v3)
   Tzz(i,j,k)=0.5*(t33-f1*v3)
   Txz(i,j,k)=0.5*(t13-f2*v1)
   Tyz(i,j,k)=0.5*(t23-f2*v2)
end do
end do
end if
end do
end if
end subroutine abs_charac

subroutine abs_tsymm
integer i,j,n,m

if (freenode) then
do n=1,num_blk
if (W(n)%nk2==nk2 .and. W(n)%isabs) then
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
   Tzz(i,j,nk2)=0.0
   Txz(i,j,nk2)=0.0
   Tyz(i,j,nk2)=0.0
do m=1,LenFD
   Tzz(i,j,nk2+m)=-Tzz(i,j,nk2-m)
   Txz(i,j,nk2+m)=-Txz(i,j,nk2-m)
   Tyz(i,j,nk2+m)=-Tyz(i,j,nk2-m)
end do
end do
end do
end if
end do
end if
end subroutine abs_tsymm
subroutine abs_vsymm
integer i,j,n,m

if (freenode) then
do n=1,num_blk
if (W(n)%nk2==nk2) then
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
do m=1,LenFD
   Vz(i,j,nk2+m)=2.0*Vz(i,j,nk2)-Vz(i,j,nk2-m)
   Vx(i,j,nk2+m)=2.0*Vx(i,j,nk2)-Vx(i,j,nk2-m)
   Vy(i,j,nk2+m)=2.0*Vy(i,j,nk2)-Vy(i,j,nk2-m)
end do
end do
end do
end if
end do
end if
end subroutine abs_vsymm

subroutine abs_extrap
integer i,j,k,n
!return

if (absnode(1,1)) then
do n=1,num_blk
if (W(n)%ni1==ni1 .and. W(n)%isabs) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=ni1-1,nx1,-1
   Vx (i,j,k)=4.0*Vx (i+1,j,k)-6.0*Vx (i+2,j,k)+4.0*Vx (i+3,j,k)-Vx (i+4,j,k)
   Vy (i,j,k)=4.0*Vy (i+1,j,k)-6.0*Vy (i+2,j,k)+4.0*Vy (i+3,j,k)-Vy (i+4,j,k)
   Vz (i,j,k)=4.0*Vz (i+1,j,k)-6.0*Vz (i+2,j,k)+4.0*Vz (i+3,j,k)-Vz (i+4,j,k)
   Txx(i,j,k)=4.0*Txx(i+1,j,k)-6.0*Txx(i+2,j,k)+4.0*Txx(i+3,j,k)-Txx(i+4,j,k)
   Tzz(i,j,k)=4.0*Tzz(i+1,j,k)-6.0*Tzz(i+2,j,k)+4.0*Tzz(i+3,j,k)-Tzz(i+4,j,k)
   Tyy(i,j,k)=4.0*Tyy(i+1,j,k)-6.0*Tyy(i+2,j,k)+4.0*Tyy(i+3,j,k)-Tyy(i+4,j,k)
   Txy(i,j,k)=4.0*Txy(i+1,j,k)-6.0*Txy(i+2,j,k)+4.0*Txy(i+3,j,k)-Txy(i+4,j,k)
   Txz(i,j,k)=4.0*Txz(i+1,j,k)-6.0*Txz(i+2,j,k)+4.0*Txz(i+3,j,k)-Txz(i+4,j,k)
   Tyz(i,j,k)=4.0*Tyz(i+1,j,k)-6.0*Tyz(i+2,j,k)+4.0*Tyz(i+3,j,k)-Tyz(i+4,j,k)
end do
end do
end do
end if
end do
end if
if (absnode(1,2)) then
do n=1,num_blk
if (W(n)%ni2==ni2 .and. W(n)%isabs) then
do k=W(n)%nk1,W(n)%nk2
do j=W(n)%nj1,W(n)%nj2
do i=ni2+1,nx2
   Vx (i,j,k)=4.0*Vx (i-1,j,k)-6.0*Vx (i-2,j,k)+4.0*Vx (i-3,j,k)-Vx (i-4,j,k)
   Vy (i,j,k)=4.0*Vy (i-1,j,k)-6.0*Vy (i-2,j,k)+4.0*Vy (i-3,j,k)-Vy (i-4,j,k)
   Vz (i,j,k)=4.0*Vz (i-1,j,k)-6.0*Vz (i-2,j,k)+4.0*Vz (i-3,j,k)-Vz (i-4,j,k)
   Txx(i,j,k)=4.0*Txx(i-1,j,k)-6.0*Txx(i-2,j,k)+4.0*Txx(i-3,j,k)-Txx(i-4,j,k)
   Tzz(i,j,k)=4.0*Tzz(i-1,j,k)-6.0*Tzz(i-2,j,k)+4.0*Tzz(i-3,j,k)-Tzz(i-4,j,k)
   Tyy(i,j,k)=4.0*Tyy(i-1,j,k)-6.0*Tyy(i-2,j,k)+4.0*Tyy(i-3,j,k)-Tyy(i-4,j,k)
   Txy(i,j,k)=4.0*Txy(i-1,j,k)-6.0*Txy(i-2,j,k)+4.0*Txy(i-3,j,k)-Txy(i-4,j,k)
   Txz(i,j,k)=4.0*Txz(i-1,j,k)-6.0*Txz(i-2,j,k)+4.0*Txz(i-3,j,k)-Txz(i-4,j,k)
   Tyz(i,j,k)=4.0*Tyz(i-1,j,k)-6.0*Tyz(i-2,j,k)+4.0*Tyz(i-3,j,k)-Tyz(i-4,j,k)
end do
end do
end do
end if
end do
end if
if (absnode(2,1)) then
do n=1,num_blk
if (W(n)%nj1==nj1 .and. W(n)%isabs) then
do k=W(n)%nk1,W(n)%nk2
do i=W(n)%ni1,W(n)%ni2
do j=nj1-1,ny1,-1
   Vx (i,j,k)=4.0*Vx (i,j+1,k)-6.0*Vx (i,j+2,k)+4.0*Vx (i,j+3,k)-Vx (i,j+4,k)
   Vy (i,j,k)=4.0*Vy (i,j+1,k)-6.0*Vy (i,j+2,k)+4.0*Vy (i,j+3,k)-Vy (i,j+4,k)
   Vz (i,j,k)=4.0*Vz (i,j+1,k)-6.0*Vz (i,j+2,k)+4.0*Vz (i,j+3,k)-Vz (i,j+4,k)
   Txx(i,j,k)=4.0*Txx(i,j+1,k)-6.0*Txx(i,j+2,k)+4.0*Txx(i,j+3,k)-Txx(i,j+4,k)
   Tzz(i,j,k)=4.0*Tzz(i,j+1,k)-6.0*Tzz(i,j+2,k)+4.0*Tzz(i,j+3,k)-Tzz(i,j+4,k)
   Tyy(i,j,k)=4.0*Tyy(i,j+1,k)-6.0*Tyy(i,j+2,k)+4.0*Tyy(i,j+3,k)-Tyy(i,j+4,k)
   Txy(i,j,k)=4.0*Txy(i,j+1,k)-6.0*Txy(i,j+2,k)+4.0*Txy(i,j+3,k)-Txy(i,j+4,k)
   Txz(i,j,k)=4.0*Txz(i,j+1,k)-6.0*Txz(i,j+2,k)+4.0*Txz(i,j+3,k)-Txz(i,j+4,k)
   Tyz(i,j,k)=4.0*Tyz(i,j+1,k)-6.0*Tyz(i,j+2,k)+4.0*Tyz(i,j+3,k)-Tyz(i,j+4,k)
end do
end do
end do
end if
end do
end if
if (absnode(2,2)) then
do n=1,num_blk
if (W(n)%nj2==nj2 .and. W(n)%isabs) then
do k=W(n)%nk1,W(n)%nk2
do i=W(n)%ni1,W(n)%ni2
do j=nj2+1,ny2
   Vx (i,j,k)=4.0*Vx (i,j-1,k)-6.0*Vx (i,j-2,k)+4.0*Vx (i,j-3,k)-Vx (i,j-4,k)
   Vy (i,j,k)=4.0*Vy (i,j-1,k)-6.0*Vy (i,j-2,k)+4.0*Vy (i,j-3,k)-Vy (i,j-4,k)
   Vz (i,j,k)=4.0*Vz (i,j-1,k)-6.0*Vz (i,j-2,k)+4.0*Vz (i,j-3,k)-Vz (i,j-4,k)
   Txx(i,j,k)=4.0*Txx(i,j-1,k)-6.0*Txx(i,j-2,k)+4.0*Txx(i,j-3,k)-Txx(i,j-4,k)
   Tzz(i,j,k)=4.0*Tzz(i,j-1,k)-6.0*Tzz(i,j-2,k)+4.0*Tzz(i,j-3,k)-Tzz(i,j-4,k)
   Tyy(i,j,k)=4.0*Tyy(i,j-1,k)-6.0*Tyy(i,j-2,k)+4.0*Tyy(i,j-3,k)-Tyy(i,j-4,k)
   Txy(i,j,k)=4.0*Txy(i,j-1,k)-6.0*Txy(i,j-2,k)+4.0*Txy(i,j-3,k)-Txy(i,j-4,k)
   Txz(i,j,k)=4.0*Txz(i,j-1,k)-6.0*Txz(i,j-2,k)+4.0*Txz(i,j-3,k)-Txz(i,j-4,k)
   Tyz(i,j,k)=4.0*Tyz(i,j-1,k)-6.0*Tyz(i,j-2,k)+4.0*Tyz(i,j-3,k)-Tyz(i,j-4,k)
end do
end do
end do
end if
end do
end if
if (absnode(3,1)) then
do n=1,num_blk
if (W(n)%nk1==nk1 .and. W(n)%isabs) then
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
do k=nk1-1,nz1,-1
   Vx (i,j,k)=4.0*Vx (i,j,k+1)-6.0*Vx (i,j,k+2)+4.0*Vx (i,j,k+3)-Vx (i,j,k+4)
   Vy (i,j,k)=4.0*Vy (i,j,k+1)-6.0*Vy (i,j,k+2)+4.0*Vy (i,j,k+3)-Vy (i,j,k+4)
   Vz (i,j,k)=4.0*Vz (i,j,k+1)-6.0*Vz (i,j,k+2)+4.0*Vz (i,j,k+3)-Vz (i,j,k+4)
   Txx(i,j,k)=4.0*Txx(i,j,k+1)-6.0*Txx(i,j,k+2)+4.0*Txx(i,j,k+3)-Txx(i,j,k+4)
   Tzz(i,j,k)=4.0*Tzz(i,j,k+1)-6.0*Tzz(i,j,k+2)+4.0*Tzz(i,j,k+3)-Tzz(i,j,k+4)
   Tyy(i,j,k)=4.0*Tyy(i,j,k+1)-6.0*Tyy(i,j,k+2)+4.0*Tyy(i,j,k+3)-Tyy(i,j,k+4)
   Txy(i,j,k)=4.0*Txy(i,j,k+1)-6.0*Txy(i,j,k+2)+4.0*Txy(i,j,k+3)-Txy(i,j,k+4)
   Txz(i,j,k)=4.0*Txz(i,j,k+1)-6.0*Txz(i,j,k+2)+4.0*Txz(i,j,k+3)-Txz(i,j,k+4)
   Tyz(i,j,k)=4.0*Tyz(i,j,k+1)-6.0*Tyz(i,j,k+2)+4.0*Tyz(i,j,k+3)-Tyz(i,j,k+4)
end do
end do
end do
end if
end do
end if
if (absnode(3,2)) then
do n=1,num_blk
if (W(n)%nk2==nk2 .and. W(n)%isabs) then
do j=W(n)%nj1,W(n)%nj2
do i=W(n)%ni1,W(n)%ni2
do k=nk2+1,nz2
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
end if
end do
end if

end subroutine abs_extrap

!--------------------------------------------------------------------}
function cal_pml_d(i,Vs,ah,nb) result(d)
integer i,nb
real(SP) Vs,d,ah
real(SP) c1,c2,c3,t
real(SP) p
real(SP) ie
ie=i
if (nb==0) then
   d=0.0
elseif (ie<0) then
   d=0.0
else
   t=3.0
   c1=8.0/15.0
   c2=-3.0/100.0
   c3=1.0/1500.0
   p=1.0
   d=t*Vs/ah *(c1+c2*nb+c3*nb**2) *(ie/real(nb,SP))**2
end if
end function cal_pml_d
function cal_pml_a(i,fc,nb) result(d)
integer i,nb
real(SP) d,fc
real(SP) ie
real(SP) dmax
ie=i
dmax=PI*fc
if (nb==0) then
   d=0.0
elseif (ie<0) then
   d=0.0
else
   d=dmax*(nb-ie)/nb
end if
end function cal_pml_a
function cal_pml_b(i,bmax,nb) result(d)
integer i,nb
real(SP) d,bmax
real(SP) ie
ie=i
if (nb==0) then
   d=1.0
elseif (ie<0) then
   d=1.0
else
   d=1.0+(bmax-1.0)*ie/nb
end if
end function cal_pml_b
function cal_pml_e(i,Vs,ah,nb) result(d)
integer i,nb
real(SP) :: Vs,d,ah
real(SP) :: ie
integer m,n
ie=i
m=(nb*ah)/(Vs*stept)
d=0.0_SP
do n=1,m
   d=d+(n*stept*Vs)**2/(nb*ah)**2
end do
d=0.8_SP/d*1.1_SP
d=exp(-d*(ie/nb)**2)
d=1.0_SP
end function cal_pml_e

subroutine coef_fdxy2fdz_abs(i,j)
  use math_mod, only : invert
  integer,intent(in) :: i,j
  real(SP),dimension(SEIS_GEO,SEIS_GEO) :: A,B,C
  real(SP) :: lam,miu,lam2mu
  real(SP) :: b1,b2,b3
  real(SP) :: e11,e12,e13,e21,e22,e23,e31,e32,e33
  integer :: k
  k=nk2
  lam=lambda(i,j,k); miu =mu(i,j,k); lam2mu=lam+2.0*miu
  b1=Bx(i); b2=By(j); b3=Bz(k)
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
       A=A/b3
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
       B=B/b1
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
       C=C/b2
       !----
       matVx2Vz(1:SEIS_GEO,1:SEIS_GEO,i,j)=matmul(A,B)
       matVy2Vz(1:SEIS_GEO,1:SEIS_GEO,i,j)=matmul(A,C)
       !----
   !MulD=A
end subroutine coef_fdxy2fdz_abs

end module abs_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
