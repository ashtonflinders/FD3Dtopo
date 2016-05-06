program seis3d_metric

! This program calculates the metric coefficient
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2007 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-14 16:30:51 -0500 (Wed, 14 Jan 2009) $
! $Revision: 510 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#include "mod_macdrp.h"

use constants_mod
use string_mod, only : string_conf
use para_mod
use math_mod
use mpi_mod
use nfseis_mod
use grid_mod
use io_mod
#ifdef MetricMPI
use mpi
#endif

implicit none

!-----------------------------------------------------------------------------

real(SP),dimension(:,:,:),allocatable :: &
     x_xi,x_eta,x_zeta,y_xi,y_eta,y_zeta,z_xi,z_eta,z_zeta
integer :: n_i,n_j,n_k
#ifdef MetricMPI
integer ierr
#endif

!-----------------------------------------------------------------------------

#ifdef MetricMPI
call MPI_INIT(ierr)
#endif

call get_conf_name(fnm_conf)

call swmpi_init(fnm_conf)
call para_init(fnm_conf)

#ifdef MetricMPI
call swmpi_cart_creat
call swmpi_datatype
call swmpi_reinit_para
#else
call swmpi_set_gindx(0,0,0)
#endif

call grid_fnm_init(fnm_conf)
call grid_alloc(iscoord=.true.,isjac=.true.,iscontra=.true.)
call alloc_covar

!-----------------------------------------------------------------------------
#ifdef MetricMPI
   n_i=thisid(1); n_j=thisid(2); n_k=thisid(3)
#else
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
#endif

#ifdef MetricMPI
if (masternode)  &
#endif
   print *, '  read in coord ...'
   call grid_coord_import(n_i,n_j,n_k)

#ifdef MetricMPI
if (masternode)  &
#endif
    print *, '  calculate metrics ...'
    call cal_e
    call cal_jac
    call cal_contravariant_e

#ifdef MetricMPI
    if (masternode)  print *, 'exchange metrics ...'
    call metric_exchange
#endif

#ifdef MetricMPI
if (masternode)  &
#endif
    print *, '  export metrics ...'
    call metric_export(n_i,n_j,n_k)

#ifndef MetricMPI
end do
end do
end do
    print *, 'exchange metrics ...'
    call metric_exchange
#endif

!-----------------------------------------------------------------------------
call dealloc_covar
call grid_dealloc

#ifdef MetricMPI
call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine alloc_covar
integer :: ierr
  allocate(  x_xi(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);   x_xi=0.0
  allocate(  y_xi(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);   y_xi=0.0
  allocate(  z_xi(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);   z_xi=0.0
  allocate( x_eta(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  x_eta=0.0
  allocate( y_eta(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  y_eta=0.0
  allocate( z_eta(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  z_eta=0.0
  allocate(x_zeta(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); x_zeta=0.0
  allocate(y_zeta(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); y_zeta=0.0
  allocate(z_zeta(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); z_zeta=0.0
end subroutine alloc_covar
subroutine dealloc_covar
  if (allocated( x_xi)) deallocate( x_xi)
  if (allocated( y_xi)) deallocate( y_xi)
  if (allocated( z_xi)) deallocate( z_xi)
  if (allocated(x_eta)) deallocate(x_eta)
  if (allocated(y_eta)) deallocate(y_eta)
  if (allocated(z_eta)) deallocate(z_eta)
  if (allocated(x_zeta)) deallocate(x_zeta)
  if (allocated(y_zeta)) deallocate(y_zeta)
  if (allocated(z_zeta)) deallocate(z_zeta)
end subroutine dealloc_covar

#ifdef MetricMPI
subroutine metric_exchange
integer,dimension(MPI_STATUS_SIZE) :: istatus
integer ierr
call MPI_SENDRECV(jac(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),411,        &
                  jac(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),411,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(jac(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),412, &
                  jac(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),412, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(jac(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),421,        &
                  jac(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),421,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(jac(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),422, &
                  jac(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),422, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(jac(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),431,        &
                  jac(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),431,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(jac(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),432, &
                  jac(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),432, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_x(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),511,        &
                  xi_x(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),511,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_x(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),512, &
                  xi_x(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),512, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_x(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),521,        &
                  xi_x(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),521,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_x(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),522, &
                  xi_x(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),522, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_x(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),531,        &
                  xi_x(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),531,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_x(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),532, &
                  xi_x(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),532, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_y(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),611,        &
                  xi_y(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),611,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_y(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),612, &
                  xi_y(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),612, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_y(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),621,        &
                  xi_y(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),621,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_y(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),622, &
                  xi_y(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),622, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_y(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),631,        &
                  xi_y(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),631,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_y(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),632, &
                  xi_y(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),632, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_z(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),711,        &
                  xi_z(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),711,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_z(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),712, &
                  xi_z(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),712, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_z(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),721,        &
                  xi_z(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),721,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_z(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),722, &
                  xi_z(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),722, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_z(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),731,        &
                  xi_z(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),731,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(xi_z(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),732, &
                  xi_z(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),732, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_x(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),811,        &
                  eta_x(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),811,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_x(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),812, &
                  eta_x(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),812, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_x(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),821,        &
                  eta_x(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),821,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_x(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),822, &
                  eta_x(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),822, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_x(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),831,        &
                  eta_x(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),831,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_x(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),832, &
                  eta_x(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),832, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_y(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),911,        &
                  eta_y(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),911,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_y(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),912, &
                  eta_y(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),912, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_y(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),921,        &
                  eta_y(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),921,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_y(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),922, &
                  eta_y(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),922, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_y(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),931,        &
                  eta_y(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),931,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_y(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),932, &
                  eta_y(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),932, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_z(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),1011,        &
                  eta_z(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1011,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_z(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1012, &
                  eta_z(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),1012, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_z(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),1021,        &
                  eta_z(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),1021,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_z(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),1022, &
                  eta_z(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),1022, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_z(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),1031,        &
                  eta_z(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),1031,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(eta_z(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),1032, &
                  eta_z(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),1032, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_x(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),1111,        &
                  zeta_x(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1111,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_x(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1112, &
                  zeta_x(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),1112, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_x(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),1121,        &
                  zeta_x(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),1121,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_x(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),1122, &
                  zeta_x(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),1122, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_x(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),1131,        &
                  zeta_x(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),1131,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_x(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),1132, &
                  zeta_x(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),1132, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_y(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),1211,        &
                  zeta_y(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1211,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_y(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1212, &
                  zeta_y(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),1212, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_y(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),1221,        &
                  zeta_y(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),1221,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_y(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),1222, &
                  zeta_y(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),1222, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_y(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),1231,        &
                  zeta_y(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),1231,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_y(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),1232, &
                  zeta_y(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),1232, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_z(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),1311,        &
                  zeta_z(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1311,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_z(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1312, &
                  zeta_z(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),1312, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_z(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),1321,        &
                  zeta_z(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),1321,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_z(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),1322, &
                  zeta_z(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),1322, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_z(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),1331,        &
                  zeta_z(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),1331,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(zeta_z(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),1332, &
                  zeta_z(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),1332, &
                  SWMPI_COMM,istatus,ierr) 
end subroutine metric_exchange
#else
subroutine metric_exchange
character (len=SEIS_STRLEN) :: filenm
integer,dimension(SEIS_GEO) :: subt
integer,dimension(SEIS_GEO) ::         &
     subs_x1,subc_x1, subs_x2,subc_x2, &
     subs_y1,subc_y1, subs_y2,subc_y2, &
     subs_z1,subc_z1, subs_z2,subc_z2
integer,dimension(LenFD) ::            &
     indx_x1,indx_x2,indx_y1,indx_y2,indx_z1,indx_z2
integer :: n_i,n_j,n_k
integer :: i,j,k

subt=(/ 1,1,1 /)
subs_x1=(/ ni2+1,ny1  ,nz1  /); subc_x1=(/ LenFD , ny    , nz    /)
subs_x2=(/ nx1  ,ny1  ,nz1  /); subc_x2=(/ LenFD , ny    , nz    /)
subs_y1=(/ nx1  ,nj2+1,nz1  /); subc_y1=(/ nx    , LenFD , nz    /)
subs_y2=(/ nx1  ,ny1  ,nz1  /); subc_y2=(/ nx    , LenFD , nz    /)
subs_z1=(/ nx1  ,ny1  ,nk2+1/); subc_z1=(/ nx    , ny    , LenFD /)
subs_z2=(/ nx1  ,ny1  ,nz1  /); subc_z2=(/ nx    , ny    , LenFD /)
indx_x1=(/ (i,i=ni1,ni1+LenFD-1) /)
indx_x2=(/ (i,i=ni2-LenFD+1,ni2) /)
indx_y1=(/ (j,j=nj1,nj1+LenFD-1) /)
indx_y2=(/ (j,j=nj2-LenFD+1,nj2) /)
indx_z1=(/ (k,k=nk1,nk1+LenFD-1) /)
indx_z2=(/ (k,k=nk2-LenFD+1,nk2) /)

do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call grid_metric_import(n_i,n_j,n_k)
   ! to x1
if (n_i>0) then
   call swmpi_change_fnm(n_i-1,n_j,n_k)
   filenm=grid_metricfnm_get(n_i-1,n_j,n_k)
   call nfseis_varput(filenm,'jac'   ,jac(indx_x1,:,:)   ,subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'xi_x'  ,xi_x(indx_x1,:,:)  ,subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'xi_y'  ,xi_y(indx_x1,:,:)  ,subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'xi_z'  ,xi_z(indx_x1,:,:)  ,subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'eta_x' ,eta_x(indx_x1,:,:) ,subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'eta_y' ,eta_y(indx_x1,:,:) ,subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'eta_z' ,eta_z(indx_x1,:,:) ,subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'zeta_x',zeta_x(indx_x1,:,:),subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'zeta_y',zeta_y(indx_x1,:,:),subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'zeta_z',zeta_z(indx_x1,:,:),subs_x1,subc_x1,subt)
end if
! to x2
if (n_i<dims(1)-1) then
   call swmpi_change_fnm(n_i+1,n_j,n_k)
   filenm=grid_metricfnm_get(n_i+1,n_j,n_k)
   call nfseis_varput(filenm,'jac'   ,jac(indx_x2,:,:)   ,subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'xi_x'  ,xi_x(indx_x2,:,:)  ,subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'xi_y'  ,xi_y(indx_x2,:,:)  ,subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'xi_z'  ,xi_z(indx_x2,:,:)  ,subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'eta_x' ,eta_x(indx_x2,:,:) ,subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'eta_y' ,eta_y(indx_x2,:,:) ,subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'eta_z' ,eta_z(indx_x2,:,:) ,subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'zeta_x',zeta_x(indx_x2,:,:),subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'zeta_y',zeta_y(indx_x2,:,:),subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'zeta_z',zeta_z(indx_x2,:,:),subs_x2,subc_x2,subt)
end if
! to y1
if (n_j>0) then
   call swmpi_change_fnm(n_i,n_j-1,n_k)
   filenm=grid_metricfnm_get(n_i,n_j-1,n_k)
   call nfseis_varput(filenm,'jac'   ,jac(:,indx_y1,:)   ,subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'xi_x'  ,xi_x(:,indx_y1,:)  ,subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'xi_y'  ,xi_y(:,indx_y1,:)  ,subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'xi_z'  ,xi_z(:,indx_y1,:)  ,subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'eta_x' ,eta_x(:,indx_y1,:) ,subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'eta_y' ,eta_y(:,indx_y1,:) ,subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'eta_z' ,eta_z(:,indx_y1,:) ,subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'zeta_x',zeta_x(:,indx_y1,:),subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'zeta_y',zeta_y(:,indx_y1,:),subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'zeta_z',zeta_z(:,indx_y1,:),subs_y1,subc_y1,subt)
end if
! to y2
if (n_j<dims(2)-1) then
   call swmpi_change_fnm(n_i,n_j+1,n_k)
   filenm=grid_metricfnm_get(n_i,n_j+1,n_k)
   call nfseis_varput(filenm,'jac'   ,jac(:,indx_y2,:)   ,subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'xi_x'  ,xi_x(:,indx_y2,:)  ,subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'xi_y'  ,xi_y(:,indx_y2,:)  ,subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'xi_z'  ,xi_z(:,indx_y2,:)  ,subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'eta_x' ,eta_x(:,indx_y2,:) ,subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'eta_y' ,eta_y(:,indx_y2,:) ,subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'eta_z' ,eta_z(:,indx_y2,:) ,subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'zeta_x',zeta_x(:,indx_y2,:),subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'zeta_y',zeta_y(:,indx_y2,:),subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'zeta_z',zeta_z(:,indx_y2,:),subs_y2,subc_y2,subt)
end if
! to k1
if (n_k>0) then
   call swmpi_change_fnm(n_i,n_j,n_k-1)
   filenm=grid_metricfnm_get(n_i,n_j,n_k-1)
   call nfseis_varput(filenm,'jac'   ,jac(:,:,indx_z1)   ,subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'xi_x'  ,xi_x(:,:,indx_z1)  ,subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'xi_y'  ,xi_y(:,:,indx_z1)  ,subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'xi_z'  ,xi_z(:,:,indx_z1)  ,subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'eta_x' ,eta_x(:,:,indx_z1) ,subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'eta_y' ,eta_y(:,:,indx_z1) ,subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'eta_z' ,eta_z(:,:,indx_z1) ,subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'zeta_x',zeta_x(:,:,indx_z1),subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'zeta_y',zeta_y(:,:,indx_z1),subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'zeta_z',zeta_z(:,:,indx_z1),subs_z1,subc_z1,subt)
end if
! to k2
if (n_k<dims(3)-1) then
   call swmpi_change_fnm(n_i,n_j,n_k+1)
   filenm=grid_metricfnm_get(n_i,n_j,n_k+1)
   call nfseis_varput(filenm,'jac'   ,jac(:,:,indx_z2)   ,subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'xi_x'  ,xi_x(:,:,indx_z2)  ,subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'xi_y'  ,xi_y(:,:,indx_z2)  ,subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'xi_z'  ,xi_z(:,:,indx_z2)  ,subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'eta_x' ,eta_x(:,:,indx_z2) ,subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'eta_y' ,eta_y(:,:,indx_z2) ,subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'eta_z' ,eta_z(:,:,indx_z2) ,subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'zeta_x',zeta_x(:,:,indx_z2),subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'zeta_y',zeta_y(:,:,indx_z2),subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'zeta_z',zeta_z(:,:,indx_z2),subs_z2,subc_z2,subt)
end if
end do
end do
end do
end subroutine metric_exchange
#endif

subroutine metric_export(n_i,n_j,n_k)
integer,intent(in) :: n_i,n_j,n_k

character (len=SEIS_STRLEN) :: filenm
integer,dimension(SEIS_GEO) :: subs,subc,subt
subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
filenm=grid_metricfnm_get(n_i,n_j,n_k)
call nfseis_grid3d_skel(filenm,nx,ny,nz,"Grid metric generated by seis3d_metric" )
call nfseis_grid3d_attput(filenm,      &
     subs,subc,subt,                   &
     ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk, &
     nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz, &
     ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,    &
     ngx1,ngx2,ngy1,ngy2,ngz1,ngz2,    &
     (/ngi1,ngi2,ngj1,ngj2,ngk1,ngk2/) )
call nfseis_grid3d_addvar(filenm,'xi_x')
call nfseis_grid3d_addvar(filenm,'xi_y')
call nfseis_grid3d_addvar(filenm,'xi_z')
call nfseis_grid3d_addvar(filenm,'eta_x')
call nfseis_grid3d_addvar(filenm,'eta_y')
call nfseis_grid3d_addvar(filenm,'eta_z')
call nfseis_grid3d_addvar(filenm,'zeta_x')
call nfseis_grid3d_addvar(filenm,'zeta_y')
call nfseis_grid3d_addvar(filenm,'zeta_z')
call nfseis_grid3d_addvar(filenm,'jac')

subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
call nfseis_varput( filenm,'xi_x', xi_x, subs,subc,subt)
call nfseis_varput( filenm,'xi_y', xi_y, subs,subc,subt)
call nfseis_varput( filenm,'xi_z', xi_z, subs,subc,subt)
call nfseis_varput( filenm,'eta_x', eta_x, subs,subc,subt)
call nfseis_varput( filenm,'eta_y', eta_y, subs,subc,subt)
call nfseis_varput( filenm,'eta_z', eta_z, subs,subc,subt)
call nfseis_varput( filenm,'zeta_x', zeta_x, subs,subc,subt)
call nfseis_varput( filenm,'zeta_y', zeta_y, subs,subc,subt)
call nfseis_varput( filenm,'zeta_z', zeta_z, subs,subc,subt)
call nfseis_varput( filenm,'jac', jac, subs,subc,subt)
end subroutine metric_export

!*************************************************************************
!*                    PART-III  metric                                   *
!*************************************************************************
subroutine cal_e
integer i,j,k
integer :: vecF(5),vecB(5)
real(SP) :: coeF(5),coeB(5)

vecF=(/ -1,0,1,2,3 /)
coeF=(/ -0.30874,-0.6326,1.2330,-0.3334,0.04168 /)
vecB=(/ -3,-2,-1,0,1 /)
coeB=(/ -0.04168,0.3334,-1.2330,0.6326,0.30874 /)

do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni2
   x_xi(i,j,k) = (                    &
      dot_product(x(i+vecB,j,k),coeB) &
     +dot_product(x(i+vecF,j,k),coeF) &
     )/steph*0.5
   y_xi(i,j,k) = (                    &
      dot_product(y(i+vecB,j,k),coeB) &
     +dot_product(y(i+vecF,j,k),coeF) &
     )/steph*0.5
   z_xi(i,j,k) = (                    &
      dot_product(z(i+vecB,j,k),coeB) &
     +dot_product(z(i+vecF,j,k),coeF) &
     )/steph*0.5
   x_eta(i,j,k) = (                   &
      dot_product(x(i,j+vecB,k),coeB) &
     +dot_product(x(i,j+vecF,k),coeF) &
     )/steph*0.5
   y_eta(i,j,k) = (                   &
      dot_product(y(i,j+vecB,k),coeB) &
     +dot_product(y(i,j+vecF,k),coeF) &
     )/steph*0.5
   z_eta(i,j,k) = (                   &
      dot_product(z(i,j+vecB,k),coeB) &
     +dot_product(z(i,j+vecF,k),coeF) &
     )/steph*0.5
   x_zeta(i,j,k) = (                  &
      dot_product(x(i,j,k+vecB),coeB) &
     +dot_product(x(i,j,k+vecF),coeF) &
     )/steph*0.5
   y_zeta(i,j,k) = (                  &
      dot_product(y(i,j,k+vecB),coeB) &
     +dot_product(y(i,j,k+vecF),coeF) &
     )/steph*0.5
   z_zeta(i,j,k) = (                  &
      dot_product(z(i,j,k+vecB),coeB) &
     +dot_product(z(i,j,k+vecF),coeF) &
     )/steph*0.5
end do
end do
end do
end subroutine cal_e

!subroutine cal_e
!integer i,j,k
!do k=nk1,nk2
!do j=nj1,nj2
!do i=ni1,ni2
!   x_xi(i,j,k) = ( &
!     m3d_FDxC1(x,i,j,k)  &
!     m3d_FDxC2(x,i,j,k)  &
!     )/steph
!   y_xi(i,j,k) = ( &
!     m3d_FDxC1(y,i,j,k)  &
!     m3d_FDxC2(y,i,j,k)  &
!     )/steph
!   z_xi(i,j,k) = ( &
!     m3d_FDxC1(z,i,j,k)  &
!     m3d_FDxC2(z,i,j,k)  &
!     )/steph
!   x_eta(i,j,k) = ( &
!     m3d_FDyC1(x,i,j,k)  &
!     m3d_FDyC2(x,i,j,k)  &
!     )/steph
!   y_eta(i,j,k) = ( &
!     m3d_FDyC1(y,i,j,k)  &
!     m3d_FDyC2(y,i,j,k)  &
!     )/steph
!   z_eta(i,j,k) = ( &
!     m3d_FDyC1(z,i,j,k)  &
!     m3d_FDyC2(z,i,j,k)  &
!     )/steph
!   x_zeta(i,j,k) = ( &
!     m3d_FDzC1(x,i,j,k)  &
!     m3d_FDzC2(x,i,j,k)  &
!     )/steph
!   y_zeta(i,j,k) = ( &
!     m3d_FDzC1(y,i,j,k)  &
!     m3d_FDzC2(y,i,j,k)  &
!     )/steph
!   z_zeta(i,j,k) = ( &
!     m3d_FDzC1(z,i,j,k)  &
!     m3d_FDzC2(z,i,j,k)  &
!     )/steph
!end do
!end do
!end do
!end subroutine cal_e
subroutine cal_jac
    real(SP),dimension(1:3) :: vec1,vec2,vec3,vecg
    integer i,j,k
    do k=nk1,nk2
    do j=nj1,nj2
    do i=ni1,ni2
       vec1=(/ x_xi(i,j,k),y_xi(i,j,k),z_xi(i,j,k) /)
       vec2=(/ x_eta(i,j,k),y_eta(i,j,k),z_eta(i,j,k) /)
       vec3=(/ x_zeta(i,j,k),y_zeta(i,j,k),z_zeta(i,j,k) /)
       call times_product(vec1,vec2,vecg)
       !jac(i,j,k)=sqrt(dot_product(vecg,vec3))
       jac(i,j,k)=dot_product(vecg,vec3)
    end do
    end do
    end do
    !call extend_metric_free(jac)
    call extend_symm(jac)
end subroutine cal_jac
subroutine cal_contravariant_e
    real(SP),dimension(1:3) :: vec1,vec2,vec3,vecg
    integer i,j,k
    do k=nk1,nk2
    do j=nj1,nj2
    do i=ni1,ni2
       vec1=(/ x_xi(i,j,k),y_xi(i,j,k),z_xi(i,j,k) /)
       vec2=(/ x_eta(i,j,k),y_eta(i,j,k),z_eta(i,j,k) /)
       vec3=(/ x_zeta(i,j,k),y_zeta(i,j,k),z_zeta(i,j,k) /)

       call times_product(vec2,vec3,vecg)
       xi_x(i,j,k)=vecg(1)/jac(i,j,k)
       xi_y(i,j,k)=vecg(2)/jac(i,j,k)
       xi_z(i,j,k)=vecg(3)/jac(i,j,k)

       call times_product(vec3,vec1,vecg)
       eta_x(i,j,k)=vecg(1)/jac(i,j,k)
       eta_y(i,j,k)=vecg(2)/jac(i,j,k)
       eta_z(i,j,k)=vecg(3)/jac(i,j,k)

       call times_product(vec1,vec2,vecg)
       zeta_x(i,j,k)=vecg(1)/jac(i,j,k)
       zeta_y(i,j,k)=vecg(2)/jac(i,j,k)
       zeta_z(i,j,k)=vecg(3)/jac(i,j,k)
    end do
    end do
    end do
   call extend_symm(  xi_x)
   call extend_symm(  xi_y)
   call extend_symm(  xi_z)
   call extend_symm( eta_x)
   call extend_symm( eta_y)
   call extend_symm( eta_z)
   call extend_symm(zeta_x)
   call extend_symm(zeta_y)
   call extend_symm(zeta_z)
end subroutine cal_contravariant_e
subroutine extend_metric_free(var)
    !replaced by extend_symm
    real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2),intent(inout) :: var
    integer k
    do k=nk2+1,nz2
       var(:,:,k)=var(:,:,2*nk2-k)
    end do
end subroutine extend_metric_free
subroutine extend_crew(w)
   real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2),intent(inout) :: w
   integer i,j,k,n
   ! x1, x2
   do k=nz1,nz2
   do j=ny1,ny2
      do n=1,LenFD
         w(ni1-n,j,k)=2.0*w(ni1,j,k)-w(ni1+n,j,k)
         w(ni2+n,j,k)=2.0*w(ni2,j,k)-w(ni2-n,j,k)
      end do
   end do
   end do
   ! y1, y2
   do k=nz1,nz2
   do i=nx1,nx2
      do n=1,LenFD
         w(i,nj1-n,k)=2.0*w(i,nj1,k)-w(i,nj1+n,k)
         w(i,nj2+n,k)=2.0*w(i,nj2,k)-w(i,nj2-n,k)
      end do
   end do
   end do
   ! z1, z2
   do j=ny1,ny2
   do i=nx1,nx2
      do n=1,LenFD
         w(i,j,nk1-n)=2.0*w(i,j,nk1)-w(i,j,nk1+n)
         w(i,j,nk2+n)=2.0*w(i,j,nk2)-w(i,j,nk2-n)
      end do
   end do
   end do
end subroutine extend_crew
subroutine extend_symm(w)
   real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2),intent(inout) :: w
   integer i,j,k,n
   ! x1, x2
   do k=nk1,nk2
   do j=nj1,nj2
      do n=1,LenFD
         w(ni1-n,j,k)=w(ni1+n,j,k)
         w(ni2+n,j,k)=w(ni2-n,j,k)
      end do
   end do
   end do
   ! y1, y2
   do k=nk1,nk2
   do i=ni1,ni2
      do n=1,LenFD
         w(i,nj1-n,k)=w(i,nj1+n,k)
         w(i,nj2+n,k)=w(i,nj2-n,k)
      end do
   end do
   end do
   ! z1, z2
   do j=nj1,nj2
   do i=ni1,ni2
      do n=1,LenFD
         w(i,j,nk1-n)=w(i,j,nk1+n)
         w(i,j,nk2+n)=w(i,j,nk2-n)
      end do
   end do
   end do
end subroutine extend_symm

!----------------- private subroutines ---------------------------------

!subroutine verify_coord_parameter(coord,hmin,hmin_idx)
!   real(SP),dimension(:,:,:),intent(in) :: coord
!   real(SP),dimension(:) :: hmin        ! minimum stephs
!   integer,dimension(:,:) :: hmin_idx    ! minimum stephs index
!   integer i,j,I1,I2,J1,J2
!   real(SP) hxmin,hymin,d
!   call scatter_IJ(I1,I2,J1,J2,comp)
!   hxmin=FD_INF; hymin=FD_INF
!   !hxmax=0.0; hymax=0.0
!   do j=J1,J2-1
!   do i=I1,I2-1
!      d=(coord(i,j,1)-coord(i+1,j,1))**2+(coord(i,j,2)-coord(i+1,j,2))**2
!      if (d<hxmin) then
!         hxmin=d; hmin_idx(1,:)=(/ i,j /)
!      end if
!      d=(coord(i,j,1)-coord(i,j+1,1))**2+(coord(i,j,2)-coord(i,j+1,2))**2
!      if (d<hymin) then
!         hymin=d; hmin_idx(2,:)=(/ i,j /)
!      end if
!   end do
!   end do
!   hmin(1)=sqrt(hxmin); hmin(2)=sqrt(hymin)
!end subroutine verify_coord_parameter

!--------------------------------------------------------

!subroutine cal_e_center(e,coord,vecNode,coeFD,ah)
!    ! user center difference to calculate coordinate partial
!    ! no use now
!    real(SP),dimension(:,:,:),intent(out) :: e
!    real(SP),dimension(:,:,:),intent(in) :: coord
!    integer i,j,I1,I2,J1,J2
!    real(SP),dimension(:) :: coeFD,ah
!    integer,dimension(:) :: vecNode
!    call scatter_IJ(I1,I2,J1,J2,comp)
!    do j=J1,J2
!    do i=I1,I2
!        e(i,j,eij2n(1,1))=dot_product(coord(i+vecNode,j,1),coeFD)/ah(1)
!        e(i,j,eij2n(2,1))=dot_product(coord(i+vecNode,j,2),coeFD)/ah(1)
!        e(i,j,eij2n(1,2))=dot_product(coord(i,j+vecNode,1),coeFD)/ah(2)
!        e(i,j,eij2n(2,2))=dot_product(coord(i,j+vecNode,2),coeFD)/ah(2)
!    end do
!    end do
!end subroutine cal_e_center

end program seis3d_metric

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
