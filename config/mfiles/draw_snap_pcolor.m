% draw_snap_pcolor: Draw wavefield snapshot by using pcolor.

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
flag_emlast=1;
flag_jetwr=1;

flag_print = 0;
flag_avi=0;

varnm='Vz'; taut=0.5;
scl_daspect=[1 1 1];
%scl_daspect=[10 10 1];
scl_caxis=[-0.5 0.5];

%id = 1; subs=[1,1,1];subc=[-1,-1,1];subt=[1,1,1];
id = 2; subs=[1,1,1];subc=[-1,1,-1];subt=[1,1,1];
%id = 3; subs=[1,1,1];subc=[1,-1,-1];subt=[1,1,1];
%id = 4; subs=[1,1,165];subc=[-1,-1,1];subt=[1,1,1];
n1=100; n2=5000; dn=100;
%n1=10; n2=n1; dn=1;

% -------------------- load data --------------------------
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
%-- get coord data
[x,y,z]=gather_coord(snapinfo,'coorddir',dir_coord);
nx=size(x,1);ny=size(x,2);nz=size(x,3);

str_unit='m';
if flag_km
   x=x/1e3;y=y/1e3;z=z/1e3; str_unit='km';
end

% ----------------------- plot figure -----------------------------------
% -- create new window --
   hid=figure;
   set(hid,'BackingStore','on');
   set(hid,'renderer','zbuffer');
   %set(hid,'menubar','none');
   set(gcf, 'PaperPositionMode', 'manual');
   set(gcf,'PaperUnits','points')
   set(gcf,'PaperPosition',[0 0 1024 768])

if flag_avi
   aviid = avifile(['snap_' num2str(id,'%3.3i') '_' varnm '.avi']);
end

% -- time loop --

for nlayer=n1:dn:n2

[v,t]=gather_snap(snapinfo,id,nlayer,varnm,'outdir',dir_out);

disp([ '  draw ' num2str(nlayer) 'th layer(t=' num2str(t) ')']);

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

%axis image
%shading interp;
shading flat;
if exist('scl_caxis'); caxis(scl_caxis); end
if exist('scl_daspect'); daspect(scl_daspect); end

if flag_jetwr; colormap(jetwr); end

colorbar('vert')

titlestr=['Snapshot of ' varnm ' at ' ...
          '{\fontsize{16}{\bf ' ...
          num2str(double(t),'%07.3f') ...
          '}}s'];
title(titlestr)

drawnow
pause(taut);

if flag_print==1
   fnm_out=[varnm '_ndim',num2str(nlayer,'%5.5i')];
   set(gca,'FontName','FixedWidth');
   print(gcf,'-dpng',[fnm_out '.png']);
end

if flag_avi==1
   F = getframe(gca);
   aviid = addframe(aviid,F);
end

end

if flag_avi==1
   aviid = close(aviid);
end

