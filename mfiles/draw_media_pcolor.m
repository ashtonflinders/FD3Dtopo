% draw_media_pcolor: Draw a cross-section of medium by using pcolor.

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Initial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date: 2008-04-27 17:31:28 -0400 (Sun, 27 Apr 2008) $
% $Revision: 469 $
% $LastChangedBy: zhangw $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

set_mfiles_path
[fnm_conf,dir_coord,dir_metric,dir_media,dir_source, ...
  dir_station,dir_out]=get_simul_path

% ----------------------- parameter -----------------------
flag_km=1;
flag_print = 0;
flag_jetwr = 0;
flag_clb=1; 
flag_title=1; 
flag_emlast=1;

% --------------------- output parameter -------------------
id = 0;

% [x,y,z] order
subs=[1,101,1];subc=[-1,1,-1];subt=[1,1,1];
scl_daspect=[1 1 1];
%scl_daspect=[5 5 1];

%varnm='rho';    %scl_caxis=[0,5000];
%varnm='mu';     %scl_caxis=[0,1e10];
%varnm='lambda'; %scl_caxis=[0,1e10];
varnm='Vp';     %scl_caxis=[0,1e10];
%varnm='Vs';     %scl_caxis=[0,1e10];

% -------------------- load data --------------------------
[snapinfo]=locate_snap(fnm_conf,0,'start',subs,'count',subc,'stride',subt);

%-- get coord data
[x,y,z]=gather_coord(snapinfo,'coorddir',dir_coord);
nx=size(x,1);ny=size(x,2);nz=size(x,3);

str_unit='m';
if flag_km
   x=x/1e3;y=y/1e3;z=z/1e3; str_unit='km';
end

switch varnm
case 'Vp'
   rho=gather_media(snapinfo,'rho','mediadir',dir_media);
   mu=gather_media(snapinfo,'mu','mediadir',dir_media);
   lambda=gather_media(snapinfo,'lambda','mediadir',dir_media);
   v=( (lambda+2*mu)./rho ).^0.5;
   v=v/1e3;
case 'Vs'
   rho=gather_media(snapinfo,'rho','mediadir',dir_media);
   mu=gather_media(snapinfo,'mu','mediadir',dir_media);
   v=( mu./rho ).^0.5;
   v=v/1e3;
case 'rho'
   v=gather_media(snapinfo,varnm,'mediadir',dir_media);
   v=v/1e3;
otherwise
   v=gather_media(snapinfo,varnm,'mediadir',dir_media);
end

% ---------------------- plot model ----------------------

hid=figure;
set(hid,'BackingStore','on');
set(hid,'renderer','zbuffer');
set(hid,'menubar','none');
set(hid,'toolbar','figure');
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf,'PaperUnits','points');
%set(gcf,'PaperPosition',[0 0 1024 768]);
%set(0, 'DefaultFigurePaperType', 'A4');

if nx==1
   if flag_emlast
      sid=pcolor(flipud(permute(squeeze(y),[2 1])), ...
                 flipud(permute(squeeze(z),[2 1])), ...
                 flipud(permute(squeeze(v),[2 1])));
   else
      sid=pcolor(permute(squeeze(y),[2 1]), ...
                 permute(squeeze(z),[2 1]), ...
                 permute(squeeze(v),[2 1]));
   end
   xlabel(['Y axis (' str_unit ')']);
   ylabel(['Z axis (' str_unit ')']);
elseif ny==1
   if flag_emlast
      %v(:,:,1)=0; v(:,:,end)=0;
      sid=pcolor(flipud(permute(squeeze(x),[2 1])), ...
                 flipud(permute(squeeze(z),[2 1])), ...
                 flipud(permute(squeeze(v),[2 1])));
   else
      sid=pcolor(permute(squeeze(x),[2 1]), ...
                 permute(squeeze(z),[2 1]), ...
                 permute(squeeze(v),[2 1]));
   end
   xlabel(['X axis (' str_unit ')']);
   ylabel(['Z axis (' str_unit ')']);
   %set(gca,'ydir','reverse');
else
   if flag_emlast
      sid=pcolor(flipud(permute(squeeze(x),[2 1])), ...
                 flipud(permute(squeeze(y),[2 1])), ...
                 flipud(permute(squeeze(v),[2 1])));
   else
      sid=pcolor(permute(squeeze(x),[2 1]), ...
                 permute(squeeze(y),[2 1]), ...
                 permute(squeeze(v),[2 1]));
   end
   xlabel(['X axis (' str_unit ')']);
   ylabel(['Y axis (' str_unit ')']);
end

% -- axis daspect --
%axis image
if exist('scl_daspect'); daspect(scl_daspect); end
axis tight

% -- colormap and colorbar
%c_spec=colormap('jetwr');
%c_spec=flipud(c_spec);
%colormap(c_spec);
if flag_jetwr; colormap(jetwr); end
if exist('scl_caxis','var'); caxis(scl_caxis); end
if flag_clb, cid=colorbar; end

% -- shading --
%shading interp;
shading flat;

if flag_title, title(varnm); end

% -------------------- save figures ------------------------
if flag_print==1
   print(gcf,'-dpng',[varnm '_ak135.png']);
end

