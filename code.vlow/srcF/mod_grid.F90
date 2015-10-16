module grid_mod

! This module contains the variables for grid geometry
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

use constants_mod
use math_mod
use string_mod
use para_mod
use mpi_mod, only :    &
    set_mpi_subfix
use nfseis_mod, only : &
    nfseis_varget

implicit none

private
public ::             &
  grid_fnm_init,      &
  grid_coordfnm_get,  &
  grid_metricfnm_get, &
  grid_coord_import,  &
  grid_metric_import, &
  grid_alloc,         &
  grid_dealloc,       &
  grid_covariant,     &
  topo_alloc,         &
  topo_import

!-----------------------------------------------------------------------------
real(SP),dimension(:,:,:),allocatable,public :: &
     xi_x,xi_y,xi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac
real(SP),dimension(:,:,:),allocatable,public :: &
     x,y,z
real(SP),dimension(:,:),allocatable,public :: &
     topox,topoy,topoz
character (len=SEIS_STRLEN),public :: &
     fnm_grid_conf,                   &
     pnm_grid,fnm_topo

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

!*************************************************************************
!*                            alloc and dealloc                          *
!*************************************************************************
subroutine grid_alloc(iscoord,isjac,iscontra)
logical,optional,intent(in) :: iscoord,iscontra,isjac
integer :: ierr
if (present(iscoord)) then
  if (iscoord)  then
  allocate( x(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  x=0.0_SP
  allocate( y(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  y=0.0_SP
  allocate( z(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  z=0.0_SP
  end if
end if

if (present(iscontra)) then
  if (iscontra)  then
  allocate(  xi_x(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  xi_x=0.0_SP
  allocate(  xi_y(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  xi_y=0.0_SP
  allocate(  xi_z(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  xi_z=0.0_SP
  allocate( eta_x(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); eta_x=0.0_SP
  allocate( eta_y(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); eta_y=0.0_SP
  allocate( eta_z(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); eta_z=0.0_SP
  allocate(zeta_x(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);zeta_x=0.0_SP
  allocate(zeta_y(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);zeta_y=0.0_SP
  allocate(zeta_z(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);zeta_z=0.0_SP
  end if
end if

if (present(isjac)) then
  if (isjac)  then
  allocate(  jac(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); jac=0.0_SP
  end if
end if
end subroutine grid_alloc

subroutine grid_dealloc(iscoord,isjac,iscontra)
logical,optional,intent(in) :: iscoord,iscontra,isjac
logical :: flagcoord,flagjac,flagcontra
if ( .not. (     present(iscontra) &
            .or. present(iscoord)  &
            .or. present(isjac) ) ) then
    flagcoord=.true.; flagjac=.true.; flagcontra=.true.
else
    flagcoord=.false.; flagjac=.false.; flagcontra=.false.
end if
if (present(iscoord)) flagcoord=iscoord
if (present(iscontra)) flagcontra=iscontra
if (present(isjac)) flagjac=isjac
   
if (flagcoord)  then
  if (allocated(x)) deallocate(x)
  if (allocated(y)) deallocate(y)
  if (allocated(z)) deallocate(z)
end if

if (flagcontra) then
  if (allocated(  xi_x)) deallocate(  xi_x)
  if (allocated(  xi_y)) deallocate(  xi_y)
  if (allocated(  xi_z)) deallocate(  xi_z)
  if (allocated( eta_x)) deallocate( eta_x)
  if (allocated( eta_y)) deallocate( eta_y)
  if (allocated( eta_z)) deallocate( eta_z)
  if (allocated(zeta_x)) deallocate(zeta_x)
  if (allocated(zeta_y)) deallocate(zeta_y)
  if (allocated(zeta_z)) deallocate(zeta_z)
end if

if (flagjac) then
  if (allocated( jac)) deallocate( jac)
end if
end subroutine grid_dealloc

subroutine topo_alloc
integer :: ierr
  allocate(topox(NTPI,NTPJ),stat=ierr);  topox=0.0_SP
  allocate(topoy(NTPI,NTPJ),stat=ierr);  topoy=0.0_SP
  allocate(topoz(NTPI,NTPJ),stat=ierr);  topoz=0.0_SP
end subroutine topo_alloc

!*************************************************************************
!*                                  grid io                              *
!*************************************************************************
subroutine grid_fnm_init(fnm_conf)
character (len=*),intent(in) :: fnm_conf
integer :: fid
fid=1001
open(fid,file=trim(fnm_conf),status="old")
  call string_conf(fid,1,'GRID_CONF',2,fnm_grid_conf)
  call string_conf(fid,1,'GRID_ROOT',2,pnm_grid)
  fnm_topo='topo.nc'
close(fid)
end subroutine grid_fnm_init

function grid_coordfnm_get(n_i,n_j,n_k) result(filenm)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_grid)//'/'//'coord'//'_'//set_mpi_subfix(n_i,n_j,n_k)//'.nc'
end function grid_coordfnm_get
function grid_metricfnm_get(n_i,n_j,n_k) result(filenm)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_grid)//'/'//'metric'//'_'//set_mpi_subfix(n_i,n_j,n_k)//'.nc'
end function grid_metricfnm_get

subroutine grid_coord_import(n_i,n_j,n_k)
  integer,intent(in) :: n_i,n_j,n_k
  integer,dimension(SEIS_GEO) ::  subs,subc,subt
  character (len=SEIS_STRLEN) :: filenm
  subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  filenm=grid_coordfnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'x', x, subs,subc,subt)
  call nfseis_varget( filenm, 'y', y, subs,subc,subt)
  call nfseis_varget( filenm, 'z', z, subs,subc,subt)
end subroutine grid_coord_import

subroutine grid_metric_import(n_i,n_j,n_k)
  integer,intent(in) :: n_i,n_j,n_k
  integer,dimension(SEIS_GEO) ::  subs,subc,subt
  character (len=SEIS_STRLEN) :: filenm
  subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  filenm=grid_metricfnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'xi_x', xi_x, subs,subc,subt)
  call nfseis_varget( filenm, 'xi_y', xi_y, subs,subc,subt)
  call nfseis_varget( filenm, 'xi_z', xi_z, subs,subc,subt)
  call nfseis_varget( filenm, 'eta_x', eta_x, subs,subc,subt)
  call nfseis_varget( filenm, 'eta_y', eta_y, subs,subc,subt)
  call nfseis_varget( filenm, 'eta_z', eta_z, subs,subc,subt)
  call nfseis_varget( filenm, 'zeta_x', zeta_x, subs,subc,subt)
  call nfseis_varget( filenm, 'zeta_y', zeta_y, subs,subc,subt)
  call nfseis_varget( filenm, 'zeta_z', zeta_z, subs,subc,subt)
  call nfseis_varget( filenm, 'jac', jac, subs,subc,subt)
end subroutine grid_metric_import

subroutine topo_import
character (len=SEIS_STRLEN) :: filenm

filenm=trim(pnm_grid)//'/'//trim(fnm_topo)
call nfseis_varget(filenm,'x',topox,(/1,1/),(/NTPI,NTPJ/),(/1,1/))
call nfseis_varget(filenm,'y',topoy,(/1,1/),(/NTPI,NTPJ/),(/1,1/))
call nfseis_varget(filenm,'z',topoz,(/1,1/),(/NTPI,NTPJ/),(/1,1/))
end subroutine topo_import

!*************************************************************************
!*                                 other                                 *
!*************************************************************************

subroutine grid_covariant(i,j,k,n,ai)
integer,intent(in) :: i,j,k,n
real(SP),dimension(SEIS_GEO),intent(out) :: ai
real(SP),dimension(SEIS_GEO) :: vec1,vec2,vec3

vec1=(/   xi_x(i,j,k),  xi_y(i,j,k),  xi_z(i,j,k) /)
vec2=(/  eta_x(i,j,k), eta_y(i,j,k), eta_z(i,j,k) /)
vec3=(/ zeta_x(i,j,k),zeta_y(i,j,k),zeta_z(i,j,k) /)
if (n==1) then
   call times_product(vec2,vec3,ai)
elseif (n==2) then
   call times_product(vec3,vec1,ai)
else
   call times_product(vec1,vec2,ai)
end if
ai=ai*jac(i,j,k)
end subroutine grid_covariant

end module grid_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
