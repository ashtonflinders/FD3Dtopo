% draw_media_print: Print figure into different elements.

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Draft

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date: 2008-04-27 17:31:28 -0400 (Sun, 27 Apr 2008) $
% $Revision: 469 $
% $LastChangedBy: zhangw $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmt='.tiff';fmtdrv='-dtiff';
%fmt='.eps';fmtdrv='-depsc2';
%fmt='.png';fmtdrv='-dpng';

%fig_dir='fig_700x300_gtopo30_sed1_crust2_ak135_USGS2x2smooth_model'
%fig_dir='fig_700x300_gtopo30_sed1_crust2_ak135_USGS2x2smooth_no3d_model'
%fig_dir='fig_700x300_gtopo30_sed1_crust2_ak135_USGS2x2smooth_flat2_model'

daspect([2,2,1]);
fig_dir='fig_700x300_gtopo30_sed1_crust2_ak135_USGS2x2smooth_model/3D221'
%daspect([5,5,1]);
%fig_dir='fig_700x300_gtopo30_sed1_crust2_ak135_USGS2x2smooth_model/3D551'
%daspect([1,1,1]);
%fig_dir='fig_700x300_gtopo30_sed1_crust2_ak135_USGS2x2smooth_model/3D111'

if ~isdir(fig_dir); mkdir(fig_dir); end

saveas(gcf,[fig_dir '/' 'Vs_3D.fig']);

shading interp;
caxis([1.2,4.6]);
xlim([-300,400]);
ylim([-200,100]);
set(gca,'xtick',[-300,-150,50,200,400]);

% ---------------------------------------------------------

% color
colormap(jet);
print(fmtdrv,[fig_dir '/' 'Vs_3D_color' fmt]);
set(gca,'visible','off');
print(fmtdrv,[fig_dir '/' 'Vs_3D_color_objonly' fmt]);

% gray
colormap(gray);
set(gca,'visible','on');
print(fmtdrv,[fig_dir '/' 'Vs_3D_gray' fmt]);
set(gca,'visible','off');
print(fmtdrv,[fig_dir '/' 'Vs_3D_gray_objonly' fmt]);

% colorbar east
figure;
colormap(jet);
caxis([1.2,4.6]);
cid=colorbar('location','east','YTick',[1.2,3,4.6]);
set(gca,'visible','off');
print(fmtdrv,[fig_dir '/' 'Vs_3D_color_clb_east' fmt]);
colorbar off

colormap(gray);
caxis([1.2,4.6]);
cid=colorbar('location','east','YTick',[1.2,3,4.6]);
print(fmtdrv,[fig_dir '/' 'Vs_3D_gray_clb_east' fmt]);

% colorbar south
figure;
colormap(jet);
caxis([1.2,4.6]);
cid=colorbar('location','south','XTick',[1.2,3,4.6]);
set(gca,'visible','off');
print(fmtdrv,[fig_dir '/' 'Vs_3D_color_clb_south' fmt]);
colorbar off

colormap(gray);
caxis([1.2,4.6]);
cid=colorbar('location','south','XTick',[1.2,3,4.6]);
print(fmtdrv,[fig_dir '/' 'Vs_3D_gray_clb_south' fmt]);
