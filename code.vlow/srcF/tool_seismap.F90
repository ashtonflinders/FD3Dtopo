program util_seismap

! This program calculates acceleration and displacement from velocity and
! their maximum distribution
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
use string_mod, only : string_conf
use math_mod
use para_mod
use mpi_mod
use nfseis_mod
use grid_mod
use media_mod
use io_mod

implicit none

integer id
logical iflag,ikeep

integer i,j,k,n,m,mt,n_i,n_j,n_k,n1,n2,n3
real t1,t2,tdelta
character (len=SEIS_STRLEN) :: fnm_nc,filenm

real(DP) :: alpha,alphasin,alphacos
real,dimension(:,:,:),allocatable :: &
     Vx,Vy,Vz,Ve,Vn,Vh,              &
     Dx,Dy,Dz,De,Dn,Dh,              &
     Ax,Ay,Az,Ae,An,Ah,              &
     VMx,VMy,VMz,VMe,VMn,VMh,VMa,    &
     DMx,DMy,DMz,DMe,DMn,DMh,DMa,    &
     AMx,AMy,AMz,AMe,AMn,AMh,AMa
integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
integer,dimension(SEIS_GEO) :: gsubs,gsubc,gsubt,gsube
integer,dimension(SEIS_GEO+1) :: subs4,subc4,subt4
integer dncid,dxid,dyid,dzid,dtid,dnid
integer ancid,axid,ayid,azid,atid,anid

!---------------------------------------------

call get_conf_name(fnm_conf)
call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)

call io_init(fnm_conf)
call io_snap_read(fnm_conf)

write(*,*) "Which snap? if the acc and disp be kept?"
read(*,*) id,ikeep

write(*,*) "alpha angle?"
read(*,*) alpha
alphasin=sin(alpha/180.0_DP*pi)
alphacos=cos(alpha/180.0_DP*pi)

n1=min(nx,snap_subc(1,id))
n2=min(ny,snap_subc(2,id))
n3=min(nz,snap_subc(3,id))
tdelta=snap_tinv(id)*stept

allocate( Vx(n1,n2,n3)); Vx=0.0
allocate( Vy(n1,n2,n3)); Vy=0.0
allocate( Vz(n1,n2,n3)); Vz=0.0
allocate( Ve(n1,n2,n3)); Ve=0.0
allocate( Vn(n1,n2,n3)); Vn=0.0
allocate( Vh(n1,n2,n3)); Vh=0.0
allocate( Ax(n1,n2,n3)); Ax=0.0
allocate( Ay(n1,n2,n3)); Ay=0.0
allocate( Az(n1,n2,n3)); Az=0.0
allocate( Ae(n1,n2,n3)); Ae=0.0
allocate( An(n1,n2,n3)); An=0.0
allocate( Ah(n1,n2,n3)); Ah=0.0
allocate( Dx(n1,n2,n3)); Dx=0.0
allocate( Dy(n1,n2,n3)); Dy=0.0
allocate( Dz(n1,n2,n3)); Dz=0.0
allocate( De(n1,n2,n3)); De=0.0
allocate( Dn(n1,n2,n3)); Dn=0.0
allocate( Dh(n1,n2,n3)); Dh=0.0
allocate(VMx(n1,n2,n3)); VMx=0.0
allocate(VMy(n1,n2,n3)); VMy=0.0
allocate(VMz(n1,n2,n3)); VMz=0.0
allocate(VMe(n1,n2,n3)); VMe=0.0
allocate(VMn(n1,n2,n3)); VMn=0.0
allocate(VMh(n1,n2,n3)); VMh=0.0
allocate(VMa(n1,n2,n3)); VMa=0.0
allocate(AMx(n1,n2,n3)); AMx=0.0
allocate(AMy(n1,n2,n3)); AMy=0.0
allocate(AMz(n1,n2,n3)); AMz=0.0
allocate(AMe(n1,n2,n3)); AMe=0.0
allocate(AMn(n1,n2,n3)); AMn=0.0
allocate(AMh(n1,n2,n3)); AMh=0.0
allocate(AMa(n1,n2,n3)); AMa=0.0
allocate(DMx(n1,n2,n3)); DMx=0.0
allocate(DMy(n1,n2,n3)); DMy=0.0
allocate(DMz(n1,n2,n3)); DMz=0.0
allocate(DMe(n1,n2,n3)); DMe=0.0
allocate(DMn(n1,n2,n3)); DMn=0.0
allocate(DMh(n1,n2,n3)); DMh=0.0
allocate(DMa(n1,n2,n3)); DMa=0.0

do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)

   call zero_local

n=0
do
   n=n+1
   fnm_nc=get_fnm_snapnode_n(pnm_out,'vel_',id,n,n_i,n_j,n_k)
   inquire(file=trim(fnm_nc),exist=iflag)
   if (.not. iflag) exit
   call nfseis_diminfo(fnm_nc,'time',mt)
   call nfseis_diminfo(fnm_nc,'I',n1)
   call nfseis_diminfo(fnm_nc,'J',n2)
   call nfseis_diminfo(fnm_nc,'K',n3)
   call nfseis_attget(fnm_nc,'subs',subs)
   call nfseis_attget(fnm_nc,'subc',subc)
   call nfseis_attget(fnm_nc,'sube',sube)
   call nfseis_attget(fnm_nc,'subt',subt)
   call nfseis_attget(fnm_nc,'gsubs',gsubs)
   call nfseis_attget(fnm_nc,'gsubc',gsubc)
   call nfseis_attget(fnm_nc,'gsube',gsube)
   call nfseis_attget(fnm_nc,'gsubt',gsubt)
   if (ikeep) then
      filenm=get_fnm_snapnode_n(pnm_out,'disp_',id,n,n_i,n_j,n_k)
      call nfseis_snap_def(filenm,dncid,dtid,tdelta,  &
           subs,subc,subt,sube,gsubs,gsubc,gsubt,gsube,  &
           title="displacement" )
      call nfseis_snap_defvar(dncid,'Dx',dxid)
      call nfseis_snap_defvar(dncid,'Dy',dyid)
      call nfseis_snap_defvar(dncid,'Dz',dzid)
      call nfseis_snap_enddef(dncid)
      filenm=get_fnm_snapnode_n(pnm_out,'acce_',id,n,n_i,n_j,n_k)
      call nfseis_snap_def(filenm,ancid,atid,tdelta,  &
           subs,subc,subt,sube,gsubs,gsubc,gsubt,gsube,  &
           title="acceleration" )
      call nfseis_snap_defvar(ancid,'Ax',axid)
      call nfseis_snap_defvar(ancid,'Ay',ayid)
      call nfseis_snap_defvar(ancid,'Az',azid)
      call nfseis_snap_enddef(ancid)
   end if
   do m=1,mt
      subs4=(/1,1,1,m/);subc4=(/n1,n2,n3,1/);subt4=(/1,1,1,1/)
      call nfseis_varget(fnm_nc,'Vx',Vx(1:n1,1:n2,1:n3), &
           subs4,subc4,subt4)
      call nfseis_varget(fnm_nc,'Vy',Vy(1:n1,1:n2,1:n3), &
           subs4,subc4,subt4)
      call nfseis_varget(fnm_nc,'Vz',Vz(1:n1,1:n2,1:n3), &
           subs4,subc4,subt4)
      call nfseis_varget(fnm_nc,'time',t2,               &
           (/      m/),(/         1/),(/      1/))
      do k=1,n3
      do j=1,n2
      do i=1,n1
         !velocity
         Ve(i,j,k)=Vx(i,j,k)*alphasin-Vy(i,j,k)*alphacos
         Vn(i,j,k)=Vx(i,j,k)*alphacos+Vy(i,j,k)*alphasin
         Vh(i,j,k)=sqrt(Vx(i,j,k)**2+Vy(i,j,k)**2)
         !displacement
         Dx(i,j,k)=Dx(i,j,k)+0.5*Vx(i,j,k)*tdelta
         Dy(i,j,k)=Dy(i,j,k)+0.5*Vy(i,j,k)*tdelta
         Dz(i,j,k)=Dz(i,j,k)+0.5*Vz(i,j,k)*tdelta
         De(i,j,k)=De(i,j,k)+0.5*Ve(i,j,k)*tdelta
         Dn(i,j,k)=Dn(i,j,k)+0.5*Vn(i,j,k)*tdelta
         Dh(i,j,k)=sqrt(Dx(i,j,k)**2+Dy(i,j,k)**2)
         !acceleration
         Ax(i,j,k)=(Vx(i,j,k)-Ax(i,j,k))/tdelta
         Ay(i,j,k)=(Vy(i,j,k)-Ay(i,j,k))/tdelta
         Az(i,j,k)=(Vz(i,j,k)-Az(i,j,k))/tdelta
         Ae(i,j,k)=(Ve(i,j,k)-Ae(i,j,k))/tdelta
         An(i,j,k)=(Vn(i,j,k)-An(i,j,k))/tdelta
         Ah(i,j,k)=sqrt(Ax(i,j,k)**2+Ay(i,j,k)**2)
         !Vmax
         if (abs(Vx(i,j,k))>VMx(i,j,k)) VMx(i,j,k)=abs(Vx(i,j,k))
         if (abs(Vy(i,j,k))>VMy(i,j,k)) VMy(i,j,k)=abs(Vy(i,j,k))
         if (abs(Vz(i,j,k))>VMz(i,j,k)) VMz(i,j,k)=abs(Vz(i,j,k))
         if (abs(Ve(i,j,k))>VMe(i,j,k)) VMe(i,j,k)=abs(Ve(i,j,k))
         if (abs(Vn(i,j,k))>VMn(i,j,k)) VMn(i,j,k)=abs(Vn(i,j,k))
         if (abs(Vh(i,j,k))>VMh(i,j,k)) then
            VMh(i,j,k)=abs(Vh(i,j,k))
            VMa(i,j,k)=xy2angle(Vx(i,j,k),Vy(i,j,k))
         end if
         if (abs(Dx(i,j,k))>DMx(i,j,k)) DMx(i,j,k)=abs(Dx(i,j,k))
         if (abs(Dy(i,j,k))>DMy(i,j,k)) DMy(i,j,k)=abs(Dy(i,j,k))
         if (abs(Dz(i,j,k))>DMz(i,j,k)) DMz(i,j,k)=abs(Dz(i,j,k))
         if (abs(De(i,j,k))>DMe(i,j,k)) DMe(i,j,k)=abs(De(i,j,k))
         if (abs(Dn(i,j,k))>DMn(i,j,k)) DMn(i,j,k)=abs(Dn(i,j,k))
         if (abs(Dh(i,j,k))>DMh(i,j,k)) then
            DMh(i,j,k)=abs(Dh(i,j,k))
            DMa(i,j,k)=xy2angle(Dx(i,j,k),Dy(i,j,k))
         end if
         if (abs(Ax(i,j,k))>AMx(i,j,k)) AMx(i,j,k)=abs(Ax(i,j,k))
         if (abs(Ay(i,j,k))>AMy(i,j,k)) AMy(i,j,k)=abs(Ay(i,j,k))
         if (abs(Az(i,j,k))>AMz(i,j,k)) AMz(i,j,k)=abs(Az(i,j,k))
         if (abs(Ae(i,j,k))>AMe(i,j,k)) AMe(i,j,k)=abs(Ae(i,j,k))
         if (abs(An(i,j,k))>AMn(i,j,k)) AMn(i,j,k)=abs(An(i,j,k))
         if (abs(Ah(i,j,k))>AMh(i,j,k)) then
            AMh(i,j,k)=abs(Ah(i,j,k))
            AMa(i,j,k)=xy2angle(Ax(i,j,k),Ay(i,j,k))
         end if
      end do
      end do
      end do
      if (ikeep) then
      call nfseis_put(dncid,dxid,Dx(1:n1,1:n2,1:n3), &
           subs4,subc4,subt4)
      call nfseis_put(dncid,dyid,Dy(1:n1,1:n2,1:n3), &
           subs4,subc4,subt4)
      call nfseis_put(dncid,dzid,Dz(1:n1,1:n2,1:n3), &
           subs4,subc4,subt4)
      call nfseis_put(dncid,dtid,t2                , &
           (/      m/),(/         1/),(/      1/))

      call nfseis_put(ancid,axid,Ax(1:n1,1:n2,1:n3), &
           subs4,subc4,subt4)
      call nfseis_put(ancid,ayid,Ay(1:n1,1:n2,1:n3), &
           subs4,subc4,subt4)
      call nfseis_put(ancid,azid,Az(1:n1,1:n2,1:n3), &
           subs4,subc4,subt4)
      call nfseis_put(ancid,atid,(t2+t1)/2.0       , &
           (/      m/),(/         1/),(/      1/))
      end if

      t1=t2
      do k=1,n3
      do j=1,n2
      do i=1,n1
         Dx(i,j,k)=Dx(i,j,k)+0.5*Vx(i,j,k)*tdelta
         Dy(i,j,k)=Dy(i,j,k)+0.5*Vy(i,j,k)*tdelta
         Dz(i,j,k)=Dz(i,j,k)+0.5*Vz(i,j,k)*tdelta
         Ax(i,j,k)=Vx(i,j,k)
         Ay(i,j,k)=Vy(i,j,k)
         Az(i,j,k)=Vz(i,j,k)
         Ae(i,j,k)=Ve(i,j,k)
         An(i,j,k)=Vn(i,j,k)
      end do
      end do
      end do
   end do
   if (ikeep) then
       call nfseis_close(dncid)
       call nfseis_close(ancid)
   end if
end do

if (n>1) then
filenm=get_fnm_snapnode_n(pnm_out,'Dmax_',id,0,n_i,n_j,n_k)
call nfseis_grid3d_skel(filenm,n1,n2,n3,title="displacement max" )
call nfseis_grid3d_addvar(filenm,'Dx')
call nfseis_grid3d_addvar(filenm,'Dy')
call nfseis_grid3d_addvar(filenm,'Dz')
call nfseis_grid3d_addvar(filenm,'De')
call nfseis_grid3d_addvar(filenm,'Dn')
call nfseis_grid3d_addvar(filenm,'Dh')
call nfseis_grid3d_addvar(filenm,'Da')
call nfseis_varput(filenm,'Dx',DMx(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Dy',DMy(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Dz',DMz(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'De',DMe(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Dn',DMn(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Dh',DMh(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Da',DMa(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))

filenm=get_fnm_snapnode_n(pnm_out,'Vmax_',id,0,n_i,n_j,n_k)
call nfseis_grid3d_skel(filenm,n1,n2,n3,title="velocity max" )
call nfseis_grid3d_addvar(filenm,'Vx')
call nfseis_grid3d_addvar(filenm,'Vy')
call nfseis_grid3d_addvar(filenm,'Vz')
call nfseis_grid3d_addvar(filenm,'Ve')
call nfseis_grid3d_addvar(filenm,'Vn')
call nfseis_grid3d_addvar(filenm,'Vh')
call nfseis_grid3d_addvar(filenm,'Va')
call nfseis_varput(filenm,'Vx',VMx(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Vy',VMy(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Vz',VMz(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Ve',VMe(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Vn',VMn(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Vh',VMh(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Va',VMa(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))

filenm=get_fnm_snapnode_n(pnm_out,'Amax_',id,0,n_i,n_j,n_k)
call nfseis_grid3d_skel(filenm,n1,n2,n3,title="acceleration max" )
call nfseis_grid3d_addvar(filenm,'Ax')
call nfseis_grid3d_addvar(filenm,'Ay')
call nfseis_grid3d_addvar(filenm,'Az')
call nfseis_grid3d_addvar(filenm,'Ae')
call nfseis_grid3d_addvar(filenm,'An')
call nfseis_grid3d_addvar(filenm,'Ah')
call nfseis_grid3d_addvar(filenm,'Aa')
call nfseis_varput(filenm,'Ax',AMx(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Ay',AMy(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Az',AMz(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Ae',AMe(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'An',AMn(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Ah',AMh(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
call nfseis_varput(filenm,'Aa',AMa(1:n1,1:n2,1:n3), &
           (/1,1,1/),(/n1,n2,n3/),(/1,1,1/))
end if

end do
end do
end do

deallocate( Vx,Vy,Vz,Vn,Ve,Vh)
deallocate( Ax,Ay,Az,An,Ae,Ah)
deallocate( Dx,Dy,Dz,Dn,De,Dh)
deallocate(VMx,VMy,VMz,VMn,VMe,VMh,VMa)
deallocate(AMx,AMy,AMz,AMn,AMe,AMh,AMa)
deallocate(DMx,DMy,DMz,DMn,DMe,DMh,DMa)

print *, 'acc and disp and max finished'

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine zero_local
Vx=0.0; Vy=0.0; Vz=0.0; Ve=0.0; Vn=0.0; Vh=0.0;
Ax=0.0; Ay=0.0; Az=0.0; Ae=0.0; An=0.0; Ah=0.0;
Dx=0.0; Dy=0.0; Dz=0.0; De=0.0; Dn=0.0; Dh=0.0;
VMx=0.0; VMy=0.0; VMz=0.0; VMe=0.0; VMn=0.0; VMh=0.0; VMa=0.0
AMx=0.0; AMy=0.0; AMz=0.0; DMe=0.0; DMn=0.0; DMh=0.0; DMa=0.0
DMx=0.0; DMy=0.0; DMz=0.0; AMe=0.0; AMn=0.0; AMh=0.0; AMa=0.0
end subroutine zero_local

function xy2angle(Lx,Ly) result (a)
real,intent(in) :: Lx,Ly
real :: a

if (abs(Lx)<SEIS_ZERO) then
   if (abs(Ly)<SEIS_ZERO) then
      a=0.0
   elseif (Ly>0.0) then
      a=90.0
   else
      a=-90.0
   end if
elseif (Lx>0) then
   if (abs(Ly)<SEIS_ZERO) then
      a=0.0
   elseif (Ly>0.0) then
      a=atan(Ly/Lx)/pi*180.0
   else
      a=atan(Ly/Lx)/pi*180.0+360.0
   end if
else
   if (abs(Ly)<SEIS_ZERO) then
      a=180.0
   else
      a=atan(Ly/Lx)/pi*180.0+180.0
   end if
end if

end function xy2angle

end program util_seismap

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
