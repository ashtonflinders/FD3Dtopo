program SI_ker_pair

! This program calculates finite-frequency sensitivity kernels for pair of
! event and station.
! The core of kernel calculation coming from the code of Z Li, C Po, Z Zhang and Y Shen.
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2008 Wei ZHANG

!*****************************************************************************
!
! $Date: 2009-01-12 19:13:22 -0500 (Mon, 12 Jan 2009) $
! $Revision: 508 $
! $LastChangedBy: zhangw $
!
!*****************************************************************************

#define TomoObsConvSrc
#define VERBOSE

use constants_mod
use string_mod
use para_mod
use io_mod
use media_mod
use nfseis_mod
use src_mod
use mpi_mod
#ifdef KernelMPI
use mpi
#endif

implicit none

character (len=SEIS_STRLEN) :: pnm_obsinfo,pnm_obs
integer id_of_obs, indx_in_obs

character (len=SEIS_STRLEN) :: &
    fnm_main_conf,             &
    fnm_tomo_conf,             &
    fnm_obs,                   &
    pnm_wave,                  &
    pnm_sgt,                   &
    pnm_ker,fnm_ker
character (len=SEIS_STRLEN) :: filenm

integer,dimension(SEIS_GEO) :: ker_blk,ker_blk_n
integer ker_nt,ker_id,knt,nt1,nt2
real(SP) :: ker_dt,ker_win(2)

real(SP),dimension(:),allocatable :: Vz,Uz,T,S
character(len=SEIS_STRLEN) :: varnm_obs
real(SP) :: Vz0,Uz0,src_m0
real(SP) :: t0_sgt
integer :: n0_sgt

real(SP),dimension(:,:,:,:),allocatable :: &
    ExxS,EyyS,EzzS,ExyS,ExzS,EyzS ,        &
    ExxR,EyyR,EzzR,ExyR,ExzR,EyzR

real(kind=SP),dimension(:),allocatable :: Ka,Kb
real(SP),dimension(:),allocatable :: E11R,E22R,E33R,E12R,E13R,E23R
real(SP),dimension(:),allocatable :: E11S,E22S,E33S,E12S,E13S,E23S
real(SP),dimension(:,:,:),allocatable :: Kap,Kaq,Kbp,Kbq

integer,allocatable :: bak_subs(:,:),bak_subt(:,:),bak_sube(:,:),bak_subc(:,:)
logical,allocatable :: bak_ishere(:)

integer n_i,n_j,n_k
integer i,j,k,m,n,mt,ierr
integer,dimension(SEIS_GEO) :: &
   bsubs,bsubc,bsubt,          &
    subs, subc, subt
integer p(1)
#ifdef KernelMPI
integer,dimension(MPI_STATUS_SIZE) :: istatus
#endif

integer ncid,ncidS,ncidF
integer kapzid,kaqzid,kbpzid,kbqzid
integer TxxSid,TyySid,TzzSid,TxySid,TxzSid,TyzSid
integer TxxFzid,TyyFzid,TzzFzid,TxyFzid,TxzFzid,TyzFzid

!----------------------------------------------------------------------

#ifdef KernelMPI
call MPI_INIT(ierr)
#else
call error_except('currently, only executable in parrallel mode')
#endif

call get_conf_name(fnm_conf)

! read kernel conf
fnm_tomo_conf='TomoKernel.conf'
call read_tomo_conf(fnm_tomo_conf)
fnm_conf=fnm_main_conf

call swmpi_init(fnm_conf)
call para_init(fnm_conf)

#ifdef KernelMPI
call swmpi_cart_creat
call swmpi_reinit_para
call swmpi_datatype
#else
call swmpi_set_gindx(0,0,0)
#endif

call media_fnm_init(fnm_conf)
!call media_import

call io_init(fnm_conf)
call io_snap_read(fnm_conf)
call io_pt_read(fnm_conf)

!call alloc_seismo(ker_nt)
!!----------------------------------------------------------------------
!! seismo on receiver
!#ifdef KernelMPI
!if (seismo_on_this(pnm_obsinfo,id_of_obs,indx_in_obs, &
!                   thisid(1),thisid(2),thisid(3),n)) then
!   filenm=get_fnm_seismo(pnm_obs,thisid(1),thisid(2),thisid(3))
!   i=snap_tinv(ker_id)/pt_tinv; ! recv type ouput may be denser than snap
!   call retrieve_recvline(filenm,n,varnm_obs,Vz,i,ker_nt,i)
!   do n=0,dims(1)*dims(2)*dims(3)-1
!   if (n/=myid) then
!      call MPI_SEND(Vz,ker_nt,MPI_REAL,n,3,MPI_COMM_WORLD,ierr)
!   end if
!   end do
!else
!   call MPI_RECV(Vz,ker_nt,MPI_REAL,MPI_ANY_SOURCE,3,MPI_COMM_WORLD,istatus,ierr)
!end if
!#else
!do n_i=0,dims(1)-1
!do n_j=0,dims(2)-1
!do n_k=0,dims(3)-1
!   call swmpi_change_fnm(n_i,n_j,n_k)
!   call swmpi_set_gindx(n_i,n_j,n_k)
!if (seismo_on_this(pnm_obsinfo,id_of_obs,indx_in_obs, &
!                   n_i,n_j,n_k,n)) then
!   filenm=get_fnm_seismo(pnm_obs,n_i,n_j,n_k)
!   i=snap_tinv(ker_id)/pt_tinv; ! recv type ouput may be denser than snap
!   call retrieve_recvline(filenm,n,varnm_obs,Vz,i,ker_nt,i)
!   exit
!end if
!end do
!end do
!end do
!#endif

ker_dt=snap_tinv(ker_id)*stept
T=(/1:ker_nt/)*ker_dt

p=minloc(T,T>=ker_win(1)); nt1=p(1)
p=maxloc(T,T<=ker_win(2)); nt2=p(1)

do n=1,ker_nt
   S(n)=src_stf(n*ker_dt,frcstf_time(1),frcstf_freq(1),frcstf_id,n)
end do
print *, "integration of stf is about ", sum(S)*ker_dt
call green_stf_halfwin(ker_dt,n0_sgt,t0_sgt)
knt=nt2+n0_sgt ! for shift green function back

#ifdef TomoObsConvSrc
call cal_convolv(S,Vz,Uz,ker_dt,1,ker_nt); Vz=Uz
nt1=nt1+n0_sgt
nt2=nt2+n0_sgt
if (nt2>nt*stept/ker_dt) then
   print *, 'nt2 too large after src conv'
   print *, 'ndim, nt2=',nt*stept/ker_dt,nt2
   print *, 'src shift=',n0_sgt
   call error_except('convolv stf failed')
end if
#endif

call vel2disp(Vz,Uz,ker_nt,ker_dt)

! normalization
Vz0=(0.5*Vz(nt1)**2+0.5*Vz(nt2)**2)*ker_dt
Uz0=(0.5*Uz(nt1)**2+0.5*Uz(nt2)**2)*ker_dt
do m=nt1+1,nt2-1
   Vz0=Vz0+Vz(m)**2.0*ker_dt
   Uz0=Uz0+Uz(m)**2.0*ker_dt
end do
Vz=Vz/Vz0/src_m0
Uz=Uz/Uz0/src_m0

call save_snap_info
!-----------------------------------------------------------------------------
#ifdef KernelMPI
   n_i=thisid(1); n_j=thisid(2); n_k=thisid(3)
#else
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
#endif

call load_snap_info
call io_snap_locate(n_i,n_j,n_k)

kernode_if : if (snap_ishere(ker_id) .and. sgt_out(ker_id)) then

call ker_para_reinit
call alloc_local(ker_blk(1),ker_blk(2),ker_blk(3),knt,nt2)
call alloc_media_local(ker_blk(1),ker_blk(2),ker_blk(3))

! create kernel nc file
fnm_ker=get_fnm_snapnode_n(pnm_ker,'kernel_'//trim(varnm_obs)//'_', &
           ker_id,0,thisid(1),thisid(2),thisid(3))
call ker_skel(fnm_ker)
call nfseis_open(fnm_ker,ncid)
call nfseis_inq_varid(ncid,'Kap',kapzid)
call nfseis_inq_varid(ncid,'Kaq',kaqzid)
call nfseis_inq_varid(ncid,'Kbp',kbpzid)
call nfseis_inq_varid(ncid,'Kbq',kbqzid)

! open modeling ouput
if (snap_tcnt(ker_id)>=knt) then
   filenm=trim(pnm_wave)                        &
       //'/'//trim(io_enum('sgt_',ker_id))      &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'_n'//trim(io_out_pattern(1,5))        &
       //'.nc'
   call nfseis_open(filenm,ncidS)
   call nfseis_inq_varid(ncidS, 'Txx', TxxSid)
   call nfseis_inq_varid(ncidS, 'Tyy', TyySid)
   call nfseis_inq_varid(ncidS, 'Tzz', TzzSid)
   call nfseis_inq_varid(ncidS, 'Txy', TxySid)
   call nfseis_inq_varid(ncidS, 'Txz', TxzSid)
   call nfseis_inq_varid(ncidS, 'Tyz', TyzSid)

   filenm=trim(pnm_sgt)                         &
       //'/'//trim(io_enum('sgt_',ker_id))      &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'_n'//trim(io_out_pattern(1,5))        &
       //'.nc'
   call nfseis_open(filenm,ncidF)
   call nfseis_inq_varid(ncidF, 'Txx', TxxFzid)
   call nfseis_inq_varid(ncidF, 'Tyy', TyyFzid)
   call nfseis_inq_varid(ncidF, 'Tzz', TzzFzid)
   call nfseis_inq_varid(ncidF, 'Txy', TxyFzid)
   call nfseis_inq_varid(ncidF, 'Txz', TxzFzid)
   call nfseis_inq_varid(ncidF, 'Tyz', TyzFzid)
end if

! loop each block
do k=1,ker_blk_n(3)
do j=1,ker_blk_n(2)
do i=1,ker_blk_n(1)

#ifdef VERBOSE
   write(*,"(a7,3(i4,a1,i4.4),a9,3(i2.2))")  &
       ' block:',              &
       i,'/',ker_blk_n(1),     &
       j,'/',ker_blk_n(2),     &
       k,'/',ker_blk_n(3),     &
       ', thisid=',n_i,n_j,n_k
#endif

   bsubs=(/ (i-1)*ker_blk(1)+1,(j-1)*ker_blk(2)+1,(k-1)*ker_blk(3)+1 /)
   bsubc=(/ min(ker_blk(1),snap_subc(1,ker_id)-bsubs(1)+1),  &
            min(ker_blk(2),snap_subc(2,ker_id)-bsubs(2)+1),  &
            min(ker_blk(3),snap_subc(3,ker_id)-bsubs(3)+1) /)
   bsubt=(/ 1,1,1 /)
   subc=bsubc; subt=snap_subt(:,ker_id)*bsubt
   subs=snap_subs(:,ker_id)+(bsubs-1)*subt

if (product(bsubc)*knt/=size(ExxS)) then
   call realloc_local(bsubc(1),bsubc(2),bsubc(3),knt)
end if

! load modeling field
if (snap_tcnt(ker_id)<knt) then
   m=0; n=0; mt=knt
   do 
      if (m>=knt) exit
      n=n+1
      filenm=trim(pnm_wave)                        &
          //'/'//trim(io_enum('sgt_',ker_id))      &
          //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
          //'_n'//trim(io_out_pattern(n,5))        &
          //'.nc'
      call nfseis_open(filenm,ncidS)
      call nfseis_inq_varid(ncidS, 'Txx', TxxSid)
      call nfseis_inq_varid(ncidS, 'Tyy', TyySid)
      call nfseis_inq_varid(ncidS, 'Tzz', TzzSid)
      call nfseis_inq_varid(ncidS, 'Txy', TxySid)
      call nfseis_inq_varid(ncidS, 'Txz', TxzSid)
      call nfseis_inq_varid(ncidS, 'Tyz', TyzSid)

      filenm=trim(pnm_sgt)                         &
          //'/'//trim(io_enum('sgt_',ker_id))      &
          //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
          //'_n'//trim(io_out_pattern(n,5))        &
          //'.nc'
      call nfseis_open(filenm,ncidF)
      call nfseis_inq_varid(ncidF, 'Txx', TxxFzid)
      call nfseis_inq_varid(ncidF, 'Tyy', TyyFzid)
      call nfseis_inq_varid(ncidF, 'Tzz', TzzFzid)
      call nfseis_inq_varid(ncidF, 'Txy', TxyFzid)
      call nfseis_inq_varid(ncidF, 'Txz', TxzFzid)
      call nfseis_inq_varid(ncidF, 'Tyz', TyzFzid)
      if (m+mt>knt) mt=knt-m

      ierr=nf90_get_var(ncidS,TxxSid,ExxS(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TyySid,EyyS(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TzzSid,EzzS(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TxySid,ExyS(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TxzSid,ExzS(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TyzSid,EyzS(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))

      ierr=nf90_get_var(ncidF,TxxFzid,ExxR(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyyFzid,EyyR(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TzzFzid,EzzR(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxyFzid,ExyR(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxzFzid,ExzR(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyzFzid,EyzR(:,:,:,m+1:m+mt), &
           (/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      if (ierr /= nf90_noerr) &
      call nfseis_except(ierr,'tomo_kernel_mpi: get equavilent stress fail')
      m=m+mt
   end do
   call nfseis_close(ncidS)
   call nfseis_close(ncidF)
else
   ierr=nf90_get_var(ncidS,TxxSid,ExxS,        &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TyySid,EyyS,        &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TzzSid,EzzS,        &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TxySid,ExyS,        &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TxzSid,ExzS,        &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TyzSid,EyzS,        &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxxFzid,ExxR,       &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyyFzid,EyyR,       &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TzzFzid,EzzR,       &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxyFzid,ExyR,       &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxzFzid,ExzR,       &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyzFzid,EyzR,       &
        (/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
end if

! load media
  filenm=media_fnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'lambda', lambda, subs,subc,subt)
  call nfseis_varget( filenm, 'mu', mu, subs,subc,subt)

! convert to stain
  call stress2strain(ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,lambda,mu,bsubc)
  call stress2strain(ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,lambda,mu,bsubc)

#ifndef TomoObsConvSrc
! correct green stress caused by Fz
   call corr_green(ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,knt,n0_sgt)
#endif

! cal kernel
   call cal_kernel(Vz,Uz,              &
        ExxR,EyyR,EzzR,ExyR,ExzR,EyzR, &
        ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
        lambda,mu,                     &
        nt1,nt2,ker_dt,bsubc,          &
        Kap,Kaq,Kbp,Kbq)
! put
   call nfseis_put(ncid,kapzid,Kap,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,kaqzid,Kaq,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,kbpzid,Kbp,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,kbqzid,Kbq,bsubs,bsubc,bsubt)

end do
end do
end do

if (snap_tcnt(ker_id)>=knt) then
   call nfseis_close(ncidS)
   call nfseis_close(ncidF)
end if

   call nfseis_close(ncid)

end if kernode_if

#ifndef KernelMPI
end do
end do
end do
#endif

!-----------------------------------------------------------------------------
call media_destroy
call dealloc_local

#ifdef KernelMPI
call MPI_BARRIER(SWMPI_COMM,ierr)
call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------------------!
contains
!-----------------------------------------------------------------------!

subroutine read_tomo_conf(fnm_conf)
character (len=*),intent(in) :: fnm_conf
character (len=SEIS_STRLEN) :: fnm_src
character (len=SEIS_STRLEN) :: stf_type
integer fid,n
fid=5002
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'MAIN_CONF',2,fnm_main_conf)

call string_conf(fid,1,'OBS_INPUT_FILE',2,fnm_obs)
call string_conf(fid,1,'component_name',2,varnm_obs)
call string_conf(fid,1,'time_window',2,ker_win(1))
call string_conf(fid,1,'time_window',3,ker_win(2))

call string_conf(fid,1,'WAVE_ROOT',2,pnm_wave)
call string_conf(fid,1,'SGT_ROOT',2,pnm_sgt)
call string_conf(fid,1,'SGT_SOURCE_CONF',2,fnm_src)
call string_conf(fid,1,'sgt_source_M0',2,src_m0)

call string_conf(fid,1,'KERNEL_ROOT',2,pnm_ker)
call string_conf(fid,1,'kernel_snap_id',2,ker_id)
!call string_conf(fid,1,'kernel_snap_nt',2,ker_nt)
do n=1,SEIS_GEO
 call string_conf(fid,1,'kernel_block',n+1,ker_blk(n))
end do
close(fid)

open(fid,file=trim(fnm_src),status='old')
   call string_conf(fid,1,"number_of_force_source",2,num_force)
   call string_conf(fid,1,"force_stf_window",2,ntwin_force)
   if (num_force/=1 .or. ntwin_force/=1) then
      call error_except("SGT should only use 1 force source and 1 time window")
   end if
   call string_conf(fid,1,"force_stf_type",2,stf_type)
   frcstf_id=stf_name2id(trim(stf_type))
   call src_alloc_force(num_force,ntwin_force)
   do n=1,ntwin_force
      call string_conf(fid,1,"force_stf_timefactor",n+1,frcstf_time(n))
      call string_conf(fid,1,"force_stf_freqfactor",n+1,frcstf_freq(n))
   end do
close(fid)

open(fid,file=trim(fnm_obs),status='old')
  read(fid,*) ker_nt
  call alloc_seismo(ker_nt)
  do n=1,ker_nt
     read(fid,*) T(n),Vz(n)
  end do
close(fid)

end subroutine read_tomo_conf

subroutine ker_para_reinit
where (ker_blk==-1)
   ker_blk=snap_subc(:,ker_id)  
end where

ker_blk_n=(snap_subc(:,ker_id)+ker_blk-1)/ker_blk

end subroutine ker_para_reinit

subroutine alloc_seismo(nt)
  integer,intent(in) :: nt
  allocate(Vz(nt)); Vz=0.0
  allocate(Uz(nt)); Uz=0.0
  allocate(T (nt)); T =0.0
  allocate(S (nt)); S =0.0
end subroutine alloc_seismo

subroutine alloc_local(ki,kj,kk,kt,nt2)
  integer,intent(in) :: ki,kj,kk,kt,nt2
  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0
  allocate(ExxR(ki,kj,kk,kt)); ExxR=0.0
  allocate(EyyR(ki,kj,kk,kt)); EyyR=0.0
  allocate(EzzR(ki,kj,kk,kt)); EzzR=0.0
  allocate(ExyR(ki,kj,kk,kt)); ExyR=0.0
  allocate(ExzR(ki,kj,kk,kt)); ExzR=0.0
  allocate(EyzR(ki,kj,kk,kt)); EyzR=0.0

  allocate(Kap(ki,kj,kk   )); Kap=0.0
  allocate(Kaq(ki,kj,kk   )); Kaq=0.0
  allocate(Kbp(ki,kj,kk   )); Kbp=0.0
  allocate(Kbq(ki,kj,kk   )); Kbq=0.0

  allocate(Ka(knt),Kb(knt));Ka=0.0;Kb=0.0
  allocate(E11S(nt2)); E11S=0.0
  allocate(E22S(nt2)); E22S=0.0
  allocate(E33S(nt2)); E33S=0.0
  allocate(E12S(nt2)); E12S=0.0
  allocate(E13S(nt2)); E13S=0.0
  allocate(E23S(nt2)); E23S=0.0
  allocate(E11R(nt2)); E11R=0.0
  allocate(E22R(nt2)); E22R=0.0
  allocate(E33R(nt2)); E33R=0.0
  allocate(E12R(nt2)); E12R=0.0
  allocate(E13R(nt2)); E13R=0.0
  allocate(E23R(nt2)); E23R=0.0
end subroutine alloc_local

subroutine alloc_media_local(ki,kj,kk)
  integer,intent(in) :: ki,kj,kk
  allocate(lambda(ki,kj,kk)); lambda=0.0
  allocate(mu(ki,kj,kk)); mu=0.0
end subroutine alloc_media_local

subroutine dealloc_local
  if (allocated(Vz)) deallocate(Vz)
  if (allocated(Uz)) deallocate(Uz)
  if (allocated(T )) deallocate(T )
  if (allocated(ExxS)) deallocate(ExxS)
  if (allocated(EyyS)) deallocate(EyyS)
  if (allocated(EzzS)) deallocate(EzzS)
  if (allocated(ExyS)) deallocate(ExyS)
  if (allocated(ExzS)) deallocate(ExzS)
  if (allocated(EyzS)) deallocate(EyzS)
  if (allocated(ExxR)) deallocate(ExxR)
  if (allocated(EyyR)) deallocate(EyyR)
  if (allocated(EzzR)) deallocate(EzzR)
  if (allocated(ExyR)) deallocate(ExyR)
  if (allocated(ExzR)) deallocate(ExzR)
  if (allocated(EyzR)) deallocate(EyzR)
  if (allocated(Ka)) deallocate(Ka)
  if (allocated(Kb)) deallocate(Kb)
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)
end subroutine dealloc_local

subroutine realloc_local(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
  if (allocated(ExxS)) deallocate(ExxS)
  if (allocated(EyyS)) deallocate(EyyS)
  if (allocated(EzzS)) deallocate(EzzS)
  if (allocated(ExyS)) deallocate(ExyS)
  if (allocated(ExzS)) deallocate(ExzS)
  if (allocated(EyzS)) deallocate(EyzS)
  if (allocated(ExxR)) deallocate(ExxR)
  if (allocated(EyyR)) deallocate(EyyR)
  if (allocated(EzzR)) deallocate(EzzR)
  if (allocated(ExyR)) deallocate(ExyR)
  if (allocated(ExzR)) deallocate(ExzR)
  if (allocated(EyzR)) deallocate(EyzR)
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)

  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0
  allocate(ExxR(ki,kj,kk,kt)); ExxR=0.0
  allocate(EyyR(ki,kj,kk,kt)); EyyR=0.0
  allocate(EzzR(ki,kj,kk,kt)); EzzR=0.0
  allocate(ExyR(ki,kj,kk,kt)); ExyR=0.0
  allocate(ExzR(ki,kj,kk,kt)); ExzR=0.0
  allocate(EyzR(ki,kj,kk,kt)); EyzR=0.0

  allocate(Kap(ki,kj,kk   )); Kap=0.0 
  allocate(Kaq(ki,kj,kk   )); Kaq=0.0 
  allocate(Kbp(ki,kj,kk   )); Kbp=0.0 
  allocate(Kbq(ki,kj,kk   )); Kbq=0.0 
end subroutine realloc_local

subroutine ker_skel(filenm)
character (len=*),intent(in) :: filenm
integer,dimension(SEIS_GEO) :: subs,subc,sube
subs=snap_subs(:,ker_id); subc=snap_subc(:,ker_id);sube=snap_sube(:,ker_id)
subs(1) =out_i(subs(1));subs(2)=out_j(subs(2));subs(3)=out_k(subs(3))
sube(1) =out_i(sube(1));sube(2)=out_j(sube(2));sube(3)=out_k(sube(3))

call nfseis_grid3d_skel(filenm,subc(1),subc(2),subc(3), &
        "Finite frequency kernel of component "//trim(varnm_obs))
call nfseis_attput(filenm,'subs',subs)
call nfseis_attput(filenm,'subc',snap_subc(:,ker_id))
call nfseis_attput(filenm,'subt',snap_subt(:,ker_id))
call nfseis_attput(filenm,'sube',sube)
call nfseis_attput(filenm,'gsubs',snap_gsubs(:,ker_id))
call nfseis_attput(filenm,'gsubc',snap_gsubc(:,ker_id))
call nfseis_attput(filenm,'gsubt',snap_gsubt(:,ker_id))
call nfseis_attput(filenm,'gsube',snap_gsube(:,ker_id))
call nfseis_grid3d_addvar(filenm,'Kap')
call nfseis_grid3d_addvar(filenm,'Kaq')
call nfseis_grid3d_addvar(filenm,'Kbp')
call nfseis_grid3d_addvar(filenm,'Kbq')
end subroutine ker_skel

subroutine cal_kernel(V,U,                &
           ExxR,EyyR,EzzR,ExyR,ExzR,EyzR, &
           ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
           lam,miu,                       &
           nt1,nt2,dt,subc,Kap,Kaq,Kbp,Kbq)

integer,intent(in) :: nt1,nt2
real(SP),dimension(:),intent(in) :: V,U
real(SP),dimension(:,:,:,:),intent(in) ::         &
     ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
     ExxR,EyyR,EzzR,ExyR,ExzR,EyzR
real(SP),dimension(:,:,:),intent(in) :: miu,lam
real(SP),dimension(:,:,:),intent(out) :: Kap,Kaq,Kbp,Kbq
integer,dimension(SEIS_GEO),intent(in) :: subc
real(SP) :: dt

integer m,n,i,j,k
real(SP) Vprho,Vsrho
real(kind=SP) :: x1,x2

do k=1,subc(3)
do j=1,subc(2)
do i=1,subc(1)
   E11R=ExxR(i,j,k,1:nt2);E22R=EyyR(i,j,k,1:nt2);E33R=EzzR(i,j,k,1:nt2)
   E12R=ExyR(i,j,k,1:nt2);E13R=ExzR(i,j,k,1:nt2);E23R=EyzR(i,j,k,1:nt2)
   E11S=ExxS(i,j,k,1:nt2);E22S=EyyS(i,j,k,1:nt2);E33S=EzzS(i,j,k,1:nt2)
   E12S=ExyS(i,j,k,1:nt2);E13S=ExzS(i,j,k,1:nt2);E23S=EyzS(i,j,k,1:nt2)

do n=nt1,nt2
   x1=0.0_SP; x2=0.0_SP
do m=1,n-1
   x1=x1+1.0_SP*dt*(E11R(m)+E22R(m)+E33R(m))*(E11S(n-m)+E22S(n-m)+E33S(n-m))
   x2=x2+1.0_SP*dt*( E11R(m)*E11S(n-m)+E22R(m)*E22S(n-m)+E33R(m)*E33S(n-m)     &
               +2.0_SP*E12R(m)*E12S(n-m)+2.0_SP*E13R(m)*E13S(n-m)+2.0_SP*E23R(m)*E23S(n-m) )
end do
Ka(n)=x1; Kb(n)=2.0_SP*(x2-x1)
end do

!   Forming GSDF kernels for tau_p and tau_q.
Kap(i,j,k)=0.5*dt*(V(nt1)*Ka(nt1)+V(nt2)*Ka(nt2))
Kaq(i,j,k)=0.5*dt*(U(nt1)*Ka(nt1)+U(nt2)*Ka(nt2))
Kbp(i,j,k)=0.5*dt*(V(nt1)*Kb(nt1)+V(nt2)*Kb(nt2))
Kbq(i,j,k)=0.5*dt*(U(nt1)*Kb(nt1)+U(nt2)*Kb(nt2))
do n=nt1+1,nt2-1
   Kap(i,j,k)=Kap(i,j,k)+V(n)*Ka(n)*dt
   Kaq(i,j,k)=Kaq(i,j,k)+U(n)*Ka(n)*dt
   Kbp(i,j,k)=Kbp(i,j,k)+V(n)*Kb(n)*dt
   Kbq(i,j,k)=Kbq(i,j,k)+U(n)*Kb(n)*dt
end do

Vprho= lam(i,j,k)+2.0*miu(i,j,k); Vsrho= miu(i,j,k)
Kap(i,j,k)= Kap(i,j,k)*2.0*Vprho
Kaq(i,j,k)=-Kaq(i,j,k)*2.0*Vprho
Kbp(i,j,k)= Kbp(i,j,k)*2.0*Vsrho
Kbq(i,j,k)=-Kbq(i,j,k)*2.0*Vsrho

end do
end do
end do
end subroutine cal_kernel

subroutine cal_convolv(U,V,W,dt,nt1,nt2)
integer,intent(in) :: nt1,nt2
real(SP),dimension(:),intent(in) :: U,V
real(SP),dimension(:),intent(out) :: W
real(SP),intent(in) :: dt
integer i,j

!W=0.0
do i=nt1,nt2
   !W(i)=0.5*dt*U(0)*V(i)+0.5*dt*U(i)*V(0)
   W(i)=0.0_SP
   do j=1,i-1
      W(i)=W(i)+dt*U(j)*V(i-j)
   end do
end do
end subroutine cal_convolv

subroutine vel2disp(V,U,ntpt,dt)
real(SP),dimension(:),intent(in) :: V
real(SP),dimension(:),intent(out) :: U
integer,intent(in) :: ntpt
real(SP),intent(in) :: dt
integer n

!U=(sum(V(1:nt-1))+V(nt)*0.5)*dt
!U(0)=0.0
U(1)=0.5*V(1)*dt
do n=2,ntpt
   U(n)=U(n-1)+0.5*(V(n-1)+V(n))*dt
end do
end subroutine vel2disp

subroutine corr_green(Txx,Tyy,Tzz,Txy,Txz,Tyz,nlen,nshift)
real(SP),dimension(:,:,:,:),intent(inout) :: Txx,Tyy,Tzz,Txy,Txz,Tyz
integer,intent(in) :: nlen,nshift
Txx(:,:,:,1:nlen-nshift)=Txx(:,:,:,nshift+1:nlen)
Tyy(:,:,:,1:nlen-nshift)=Tyy(:,:,:,nshift+1:nlen)
Tzz(:,:,:,1:nlen-nshift)=Tzz(:,:,:,nshift+1:nlen)
Txy(:,:,:,1:nlen-nshift)=Txy(:,:,:,nshift+1:nlen)
Txz(:,:,:,1:nlen-nshift)=Txz(:,:,:,nshift+1:nlen)
Tyz(:,:,:,1:nlen-nshift)=Tyz(:,:,:,nshift+1:nlen)
end subroutine corr_green

subroutine stress2strain(Exx,Eyy,Ezz,Exy,Exz,Eyz,lam,miu,scl)
real(SP),dimension(:,:,:,:),intent(inout) :: Exx,Eyy,Ezz,Exy,Exz,Eyz
real(SP),dimension(:,:,:),intent(in) :: lam,miu
integer,dimension(SEIS_GEO),intent(in) :: scl
real(SP) E1,E2,E3,E0
integer i,j,k,n,n0

n0=size(Exx,4)
do n=1,n0
do k=1,scl(3)
do j=1,scl(2)
do i=1,scl(1)
   E1=(lam(i,j,k)+miu(i,j,k))/(miu(i,j,k)*(3.0*lam(i,j,k)+2.0*miu(i,j,k)))
   E2=-lam(i,j,k)/(2.0*miu(i,j,k)*(3.0*lam(i,j,k)+2.0*miu(i,j,k)))
   E3=1.0/miu(i,j,k)
   E0=E2*(Exx(i,j,k,n)+Eyy(i,j,k,n)+Ezz(i,j,k,n))
   Exx(i,j,k,n)=E0-(E2-E1)*Exx(i,j,k,n)
   Eyy(i,j,k,n)=E0-(E2-E1)*Eyy(i,j,k,n)
   Ezz(i,j,k,n)=E0-(E2-E1)*Ezz(i,j,k,n)
   Exy(i,j,k,n)=0.5*E3*Exy(i,j,k,n)
   Exz(i,j,k,n)=0.5*E3*Exz(i,j,k,n)
   Eyz(i,j,k,n)=0.5*E3*Eyz(i,j,k,n)
end do
end do
end do
end do
end subroutine stress2strain

subroutine green_stf_halfwin(dt,n0,t0)
real(SP),intent(in) :: dt
integer,intent(out) :: n0
real(SP),intent(out) :: t0

select case(frcstf_id)
case (SIG_STF_BELL)
   t0=frcstf_freq(1)/2.0_SP
case (SIG_STF_RICKER)
   t0=frcstf_time(1)
case (SIG_STF_GAUSS)
   t0=frcstf_time(1)/2.0_SP
case default
   print *, "Need to set halfwin for stf "//trim(stf_id2name(frcstf_id))
   stop 1
end select
n0=int(n0/dt+0.5)
end subroutine green_stf_halfwin

subroutine save_snap_info
allocate(bak_subs(SEIS_GEO,num_snap))
allocate(bak_sube(SEIS_GEO,num_snap))
allocate(bak_subt(SEIS_GEO,num_snap))
allocate(bak_subc(SEIS_GEO,num_snap))
allocate(bak_ishere(num_snap))

bak_subs=snap_subs
bak_subc=snap_subc
bak_subt=snap_subt
bak_sube=snap_sube
bak_ishere=snap_ishere
end subroutine save_snap_info
subroutine load_snap_info
snap_subs  =bak_subs  
snap_subc  =bak_subc  
snap_subt  =bak_subt  
snap_sube  =bak_sube  
snap_ishere=bak_ishere
end subroutine load_snap_info

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
#ifdef KernelMPI
  integer :: ierr
#endif
  print *, trim(msg)
#ifdef KernelMPI
  call MPI_ABORT(SWMPI_COMM,1,ierr)
#else
  stop 1
#endif
end subroutine error_except

end program SI_ker_pair

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
