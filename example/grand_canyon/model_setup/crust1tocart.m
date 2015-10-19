% A. Flinders
%
% 12/02/2013
%
% This is a modified version of Y.Shen's crust2tocart.m code. It is updated
% to use crust1.0 over crust2.0 (e.g. 1x1 deg opposed to 2x2 deg). There
% has been a renaming convention when transitioning to crust1, and it
% should be noted that the crustal maps now start at index 1, and not index
% 0 (as they did in crust2.0).
%
% Crust1.0 does not extract a topographic file as crust2.0 did, but both
% however are based on etopo1. The topography/bathymetry layer in crust1.0
% is bd-2; where bd-1 (the first layer) is the topogrpahy plus seasurface
%
% This script makes the previously used sedmap2cart.m redundant, as
% crust1.0 incorporates the 1x1 degree sediment thickness maps previously
% distributed seperatly from crust2.0
%
% Similarly, this script makes redunant the previously used etopo2cart.m 
% for the same reasons.
%
% Lastly, this script has the option to take in a high res DEM to replace
% the surface topography of etopo1
%
% 04/28/2015
%
% Discovered awhile several months ago, we can not use crust1.0 as it is
% for small local scale studies with high frequency sources. The upper low
% velocity zones of crust1.0 end up trapping energy in the upper layers
% causing significant ringing.
%
% Several attempts have been made to get aorund this, including an
% isovelocity model, and setting the upper 3 layers to the same velocity.
% Unfortunatly due to the small scale, it looks like this causes the Sv and
% Rayleigh wave to coincide, possibly causing issues with the kernel
% calculations (near source high sensitiviy, dipole pattern).


clear all; close all; clc

%% higher res. topo DEM ----------------------------
if 1
    DEM_file='./grand_canyon_topo.grd'
    dem_out = 'model.hiresDEM.nc'
end
%

remove_icewater = 1;
remove_upper = 0 ;
flat_isovel = 0 ;

plot_global = 0;
mirror = 0;

complete_isovel = 0;
    S_iso=3200;
    P_iso=6000;
    D_iso=2500;
    
% Directories
%-------------------------- output model name -----------------------------
fnm_out='model.crust1.nc';

%-------------------------- mfile directories -----------------------------
addpath(genpath('../../../config/mfiles'))

%------------------------ Crust1.0 directories ----------------------------
pnm_map='../../../config/crust1/crust1-maps/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get boundary coordinates
define_xy_coords
vx=x; vy=y; clear x y
[mx,my]=meshgrid(vx,vy);

lat1=floor(minlat);
lat2=ceil(maxlat);
lon1=floor(minlon);
lon2=ceil(maxlon);

%--------------------------- read in crust1.0 -----------------------------
% Crust1 coordinates 
vlat=-89.5:89.5; nlat=length(vlat);
vlon=-179.5:179.5; 
nlon=length(vlon);
nlayer=9;

% topography and bathymetry (e.g. bottom of water but top of ice)!
fid=fopen([pnm_map '/' 'map-bd2'],'r');
topo=fscanf(fid,'%f',[nlon,nlat]); topo=topo';
fclose(fid);

% topography and seasurface
fid=fopen([pnm_map '/' 'map-bd1'],'r');
bd1=fscanf(fid,'%f',[nlon,nlat]); bd1=bd1';
fclose(fid);

% read interface/layer boundaries, we go from bd2-bd9 (bd#=n+1)
for n=1:nlayer-1
    idt=fopen([pnm_map '/' 'map-bd' num2str(n+1)],'r');
    Lb{n}=fscanf(idt,'%f',[nlon,nlat]); Lb{n}=Lb{n}'; %#ok<*SAGROW>
    if n==1, Lh{n}=bd1-Lb{n}; else Lh{n}=Lb{n-1}-Lb{n}; end
    fclose(idt);
end

for n=1:nlayer
    idvp=fopen([pnm_map '/' 'map-vp' num2str(n)],'r');
    idvs=fopen([pnm_map '/' 'map-vs' num2str(n)],'r');
    iddp=fopen([pnm_map '/' 'map-ro' num2str(n)],'r');
    Vp{n}=fscanf(idvp,'%f',[nlon,nlat]); Vp{n}=Vp{n}';
    Vs{n}=fscanf(idvs,'%f',[nlon,nlat]); Vs{n}=Vs{n}';
    Dp{n}=fscanf(iddp,'%f',[nlon,nlat]); Dp{n}=Dp{n}';
    fclose(idvp);fclose(idvs);fclose(iddp);
end

topo=flipdim(topo,1);
bd1=flipdim(bd1,1);
for n=1:nlayer-1
    Lb{n}=flipdim(Lb{n},1);
    Lh{n}=flipdim(Lh{n},1);
end

for n=1:nlayer
    Vp{n}=flipdim(Vp{n},1);
    Vs{n}=flipdim(Vs{n},1);
    Dp{n}=flipdim(Dp{n},1);
end

% convert to m
topo=topo*1e3;
bd1=bd1*1e3;
for n=1:nlayer-1
    Lb{n}=Lb{n}*1e3;Lh{n}=Lh{n}*1e3;
end
for n=1:nlayer
    Vp{n}=Vp{n}*1e3; Vs{n}=Vs{n}*1e3; Dp{n}=Dp{n}*1e3;
end

%% cut crust1.0 to boundaries 
j1=find(vlat>=lat1); j1=j1(1); j2=find(vlat<=lat2); j2=j2(end);
i1=find(vlon>=lon1); i1=i1(1); i2=find(vlon<=lon2); i2=i2(end);

vlat([1:j1-2,j2+2:end])=[];
vlon([1:i1-2,i2+2:end])=[];
topo([1:j1-2,j2+2:end],:)=[]; topo(:,[1:i1-2,i2+2:end])=[];
bd1([1:j1-2,j2+2:end],:)=[]; bd1(:,[1:i1-2,i2+2:end])=[];
for n=1:nlayer-1
    Lb{n}([1:j1-2,j2+2:end],:)=[]; Lb{n}(:,[1:i1-2,i2+2:end])=[];
    Lh{n}([1:j1-2,j2+2:end],:)=[]; Lh{n}(:,[1:i1-2,i2+2:end])=[];
end
for n=1:nlayer
    Vp{n}([1:j1-2,j2+2:end],:)=[]; Vp{n}(:,[1:i1-2,i2+2:end])=[];
    Vs{n}([1:j1-2,j2+2:end],:)=[]; Vs{n}(:,[1:i1-2,i2+2:end])=[];
    Dp{n}([1:j1-2,j2+2:end],:)=[]; Dp{n}(:,[1:i1-2,i2+2:end])=[];
end
%
%% convert to Cartesian
[lon,lat]=meshgrid(vlon,vlat);
[xlon,ylat]=geo2cart(lat,lon,lat0,lon0,alpha,x0,y0);

% interpolate to equal spacing
Z0=griddata(xlon,ylat,topo,mx,my);
BD1=griddata(xlon,ylat,bd1,mx,my);

for n=1:nlayer-1
    Z{n}=griddata(xlon,ylat,Lb{n},mx,my,'cubic');
    H{n}=griddata(xlon,ylat,Lh{n},mx,my,'cubic');
    H{n}(find(H{n}<0))=0;
end

for n=1:nlayer
    P{n}=griddata(xlon,ylat,Vp{n},mx,my,'cubic');
    S{n}=griddata(xlon,ylat,Vs{n},mx,my,'cubic');
    D{n}=griddata(xlon,ylat,Dp{n},mx,my,'cubic');
end

for n=1:nlayer
    P_near{n}=griddata(xlon,ylat,Vp{n},mx,my,'nearest');
    S_near{n}=griddata(xlon,ylat,Vs{n},mx,my,'nearest');
    D_near{n}=griddata(xlon,ylat,Dp{n},mx,my,'nearest');
end
%
%% replace less than nearest with nearest values, and zeros with means
for n=1:nlayer
    jP = find(P{n} < P_near{n});
    P{n}(jP)=P_near{n}(jP);
    
    unique_val=unique(P_near{n});
    unique_val(find(unique_val == 0))=[];
    if ~isempty(unique_val)
        jP = find(P{n} == 0 | P{n} < min(unique_val));
        P{n}(jP)=mean(unique_val);
    end
    
    jS = find(S{n} < S_near{n});
    S{n}(jS)=S_near{n}(jS);
    
    unique_val=unique(S_near{n});
    unique_val(find(unique_val == 0))=[];
    if ~isempty(unique_val)
        jS = find(S{n} == 0 | S{n} < min(unique_val));
        S{n}(jS)=mean(unique_val);
    end
    
    jD = find(D{n} < D_near{n});
    D{n}(jD)=D_near{n}(jD);
    
    unique_val=unique(D_near{n});
    unique_val(find(unique_val == 0))=[];
    if ~isempty(unique_val)
        jD = find(D{n} == 0 | D{n} < min(unique_val));
        D{n}(jD)=mean(unique_val);
    end
end
% check for nans
for n=1:nlayer-1
    if ~isempty(find(isnan(Z{n})) | find(isnan(H{n})))
       error('nan in Z or H');
    end
end
for n=1:nlayer
    if ~isempty(find(isnan(P{n})) | find(isnan(S{n})) |  find(isnan(D{n})) )
       error('nan in Vp Vs or Dp');
    end
end
% remove water and ice layer
if remove_icewater
    for n=1:2
        P{n}=[];
        S{n}=[];
        D{n}=[];
        H{n}=[];
        Z{n}=[];
    end
    Hnew=H(~cellfun(@isempty,H)); H=Hnew;
    Znew=Z(~cellfun(@isempty,Z)); Z=Znew;
    Pnew=P(~cellfun(@isempty,P)); P=Pnew;
    Snew=S(~cellfun(@isempty,S)); S=Snew;
    Dnew=D(~cellfun(@isempty,D)); D=Dnew;
end
nlayer=length(P);



%--------------------- bring in high res DEM and cut ----------------------
if exist('DEM_file')
    [lon_DEM,lat_DEM,Z_DEM]=grdread2(DEM_file);

    jd1=find(lat_DEM>=lat1); jd1=jd1(1); jd2=find(lat_DEM<=lat2); jd2=jd2(end);
    id1=find(lon_DEM>=lon1); id1=id1(1); id2=find(lon_DEM<=lon2); id2=id2(end);

    lat_DEM([1:jd1-10,jd2+10:end])=[];
    lon_DEM([1:id1-10,id2+10:end])=[];
    Z_DEM([1:jd1-10,jd2+10:end],:)=[]; Z_DEM(:,[1:id1-10,id2+10:end])=[];
    Z_DEM=double(Z_DEM);

    % convert to Cartesian
    [LON_DEM,LAT_DEM]=meshgrid(lon_DEM,lat_DEM);
    [X_DEM,Y_DEM]=geo2cart(LAT_DEM,LON_DEM,lat0,lon0,alpha,x0,y0);

    % interpolate to equal spacing
    DEM_SAM=griddata(X_DEM,Y_DEM,Z_DEM,mx,my);
    
    % make sure Z{1} and Z{2} layer depths are beneath the DEM topography; 
    % issues will have likely arrisen given the local resolution of crust1.0 not
    % acurrately representing the depths of lowlands etc...
    for n=1:nlayer-1
        check_inv=find(Z{n}>DEM_SAM);
        Znew{n}=Z{n}; Znew{n}(check_inv)=DEM_SAM(check_inv);
        Hnew{n}=H{n}; Hnew{n}(check_inv)=0;
    end
    
    % do the same check but for the other layers (e.g. make sure layers are
    % under each other
    for n=2:nlayer-1
        check_inv=find(Znew{n}>Znew{n-1});
        Znew{n}(check_inv)=Znew{n-1}(check_inv);
        Hnew{n}(check_inv)=0;
    end
    
    % Adjust the first layer thickness to account for the new topography
    % H{1} is thickness of layer between surface (topo) and first layer
    Hnew{1}=DEM_SAM-Znew{1}+Hnew{1};
    
    if mirror == 1
        % this will make the layers mirror the topograpghy
        Znew{1}=Znew{1}+(DEM_SAM-Z0);
        % Adjust the layer depths to reflect the new topography (propagate through)
        for n=2:nlayer-1
            Znew{n}=Znew{n}+(DEM_SAM-Z0);
            Hnew{n}=Znew{n-1}-Znew{n}
        end
    end
    
    % output surface DEM in projected coordinates
    if 1
        nc_create_empty(dem_out);
        nc_add_dimension(dem_out,'x',nx);
        nc_add_dimension(dem_out,'y',ny);

        var.Nctype='float';var.Attribute=[];
        var.Name='x';var.Dimension={'x'};nc_addvar(dem_out,var);
        var.Name='y';var.Dimension={'y'};nc_addvar(dem_out,var);
        var.Name='topo';var.Dimension={'y','x'};nc_addvar(dem_out,var);

        nc_varput(dem_out,'x',vx');
        nc_varput(dem_out,'y',vy');
        nc_varput(dem_out,'topo',DEM_SAM);
    end
    
    Z0=DEM_SAM;
    % have to change this if we have areas with lots of water....
    BD1=DEM_SAM;
    Z=Znew;
    H=Hnew;
    
end

%% additional controls
if flat_isovel == 1
    for n=1:nlayer
            P{n}=ones(size(P{n}))*mean(mean(P{n}));
            S{n}=ones(size(S{n}))*mean(mean(S{n}));
            D{n}=ones(size(D{n}))*mean(mean(D{n}));
    end
    
    for n=1:nlayer-1
        if n == 1
            Z{n}=ones(size(Z{n}))*min(min(Z{n}));
            H{n}=Z0-Z{n};
        else
            Z{n}=ones(size(Z{n}))*mean(mean(Z{n}));
            H{n}=Z{n-1}-Z{n};
        end
    end
end

if complete_isovel == 1
    for n=1:nlayer
        P{n}=ones(size(P{n}))*P_iso;
	    S{n}=ones(size(S{n}))*S_iso;
        D{n}=ones(size(D{n}))*D_iso;
    end         
end

%% remove layers of minimum thickness
for n=1:nlayer-1;
    if max(max(H{n})) < 100
        H{n}=[];
        Z{n}=[];
        P{n}=[];
        S{n}=[];
        D{n}=[];
    end
end
Hnew=H(~cellfun(@isempty,H)); H=Hnew;
Znew=Z(~cellfun(@isempty,Z)); Z=Znew;
Pnew=P(~cellfun(@isempty,P)); P=Pnew;
Snew=S(~cellfun(@isempty,S)); S=Snew;
Dnew=D(~cellfun(@isempty,D)); D=Dnew;

nlayer=length(P);


%% new plot
if 0
    [MX,MY,MZind]=meshgrid(vx,vy,1:size(S,2));

    half_x=size(MX,1)/2;
    half_y=size(MX,2)/2;

    S_3d=cat(3,S{:});
    P_3d=cat(3,P{:});
    Z_3d=cat(3,BD1,Z{:});
    H_3d=cat(3,H{:});
    D_3d=cat(3,D{:});

    S_clim=[min(min(min(S_3d))) max(max(max(S_3d)))];
    P_clim=[min(min(min(P_3d))) max(max(max(P_3d)))];

    figure()
    clf
    subplot(231)
    set(gcf,'renderer','zbuffer');
    surface(mx,my,DEM_SAM,'EdgeColor','none');
    title('Hi-Res Topography/Bathymetry');
    caxis([min(min(DEM_SAM)) max(max(DEM_SAM))])
    axis tight
    colorbar
    subplot(232)
    h=pcolor(permute(MY(:,half_y,:),[1 3 2]),permute(Z_3d(:,half_y,:),[1 3 2]),permute(S_3d(:,half_y,:),[1 3 2]));
    set(h,'EdgeColor','none')
    %ylim([-15 4]*1e3 )
    colorbar
    caxis(S_clim)
    subplot(233)
    h=pcolor(permute(MX(half_x,:,:),[2 3 1]),permute(Z_3d(half_x,:,:),[2 3 1]),permute(S_3d(half_x,:,:),[2 3 1]));
    set(h,'EdgeColor','none')
    %ylim([-15 4]*1e3 )
    colorbar
    caxis(S_clim)
    subplot(234)
    surf(mx,my,DEM_SAM,'EdgeColor','none');
    hold on
    for n=1:size(Z,2)
        surf(MX(:,:,n),MY(:,:,n),Z{n},'EdgeColor','none')
    end
    subplot(235)
    h=pcolor(permute(MY(:,half_y,:),[1 3 2]),permute(Z_3d(:,half_y,:),[1 3 2]),permute(P_3d(:,half_y,:),[1 3 2]));
    set(h,'EdgeColor','none')
    %ylim([-15 4]*1e3 )
    colorbar
    caxis(P_clim)
    subplot(236)
    h=pcolor(permute(MX(half_x,:,:),[2 3 1]),permute(Z_3d(half_x,:,:),[2 3 1]),permute(P_3d(half_x,:,:),[2 3 1]));
    set(h,'EdgeColor','none')
    %ylim([-15 4]*1e3 )
    colorbar
    caxis(P_clim)
    colormap(flipud(jet))
end
%% 
% output surface DEM in projected coordinates, just for later visulization
if exist('dem_out')
    nc_create_empty(dem_out);
    nc_add_dimension(dem_out,'x',nx);
    nc_add_dimension(dem_out,'y',ny);

    var.Nctype='float';var.Attribute=[];
    var.Name='x';var.Dimension={'x'};nc_addvar(dem_out,var);
    var.Name='y';var.Dimension={'y'};nc_addvar(dem_out,var);
    var.Name='topo';var.Dimension={'y','x'};nc_addvar(dem_out,var);

    nc_varput(dem_out,'x',vx');
    nc_varput(dem_out,'y',vy');
    nc_varput(dem_out,'topo',DEM_SAM);
end

%%

C1Vp=permute(cat(3,P{:}),[3 1 2]);
C1Vs=permute(cat(3,S{:}),[3 1 2]);
C1Dp=permute(cat(3,D{:}),[3 1 2]);
C1H=permute(cat(3,H{:}),[3 1 2]);
C1Z=permute(cat(3,Z{:}),[3 1 2]);

nlayer=size(C1Vp,1)-1;
Vp=zeros([nlayer,size(C1Vp,2),size(C1Vp,3),2]);
Vs=zeros([nlayer,size(C1Vp,2),size(C1Vp,3),2]);
Dp=zeros([nlayer,size(C1Vp,2),size(C1Vp,3),2]);

for n=1:nlayer
    Vp(n,:,:,1)=C1Vp(n,:,:); %velocity at the top and (next) bottom of the layer
    Vp(n,:,:,2)=C1Vp(n,:,:);
    Vs(n,:,:,1)=C1Vs(n,:,:);
    Vs(n,:,:,2)=C1Vs(n,:,:);
    Dp(n,:,:,1)=C1Dp(n,:,:);
    Dp(n,:,:,2)=C1Dp(n,:,:);
end

% extend down to 120 km
C1Z(nlayer-1,:,:)=-120e3*ones(size(Vp,2),size(Vp,3));
C1H(nlayer-1,:,:)=C1Z(nlayer-2,:,:)-C1Z(nlayer-1,:,:);
%-------------------------  conf --------------------------------------------

QsF0=1.0; QsINF=1000.0;
Qs=ones(size(Vp))*1e3;

Vp_poly_d=ones(1,nlayer);
Vs_poly_d=ones(1,nlayer);
rho_poly_d=ones(1,nlayer);
Qs_poly_d=ones(1,nlayer);

%---------------------------- create nc file -------------------------------------
if 1
    nc_create_empty(fnm_out);
    nc_add_dimension(fnm_out,'x',nx);
    nc_add_dimension(fnm_out,'y',ny);
    nc_add_dimension(fnm_out,'layer',nlayer);
    nc_add_dimension(fnm_out,'side',2);
    
    nc_attput(fnm_out,nc_global,'QsF0',QsF0);
    nc_attput(fnm_out,nc_global,'QsINF',QsINF);
    
    var.Nctype='float';var.Attribute=[];
    var.Name='x';var.Dimension={'x'};nc_addvar(fnm_out,var);
    var.Name='y';var.Dimension={'y'};nc_addvar(fnm_out,var);
    var.Name='Vp_poly_d';var.Dimension={'layer'};nc_addvar(fnm_out,var);
    var.Name='Vs_poly_d';var.Dimension={'layer'};nc_addvar(fnm_out,var);
    var.Name='rho_poly_d';var.Dimension={'layer'};nc_addvar(fnm_out,var);
    var.Name='Qs_poly_d';var.Dimension={'layer'};nc_addvar(fnm_out,var);
    var.Name='thickness';var.Dimension={'layer','y','x'};nc_addvar(fnm_out,var);
    var.Name='Vp'; var.Dimension={'layer','y','x','side'};nc_addvar(fnm_out,var);
    var.Name='Vs'; var.Dimension={'layer','y','x','side'};nc_addvar(fnm_out,var);
    var.Name='rho'; var.Dimension={'layer','y','x','side'};nc_addvar(fnm_out,var);
    var.Name='Qs'; var.Dimension={'layer','y','x','side'};nc_addvar(fnm_out,var);
    
    nc_varput(fnm_out,'x',vx);
    nc_varput(fnm_out,'y',vy);
    nc_varput(fnm_out,'Vp_poly_d',Vp_poly_d');
    nc_varput(fnm_out,'Vs_poly_d',Vs_poly_d');
    nc_varput(fnm_out,'rho_poly_d',rho_poly_d');
    nc_varput(fnm_out,'Qs_poly_d',Qs_poly_d');
    
    nc_varput(fnm_out,'thickness',C1H);
    nc_varput(fnm_out,'Vp',Vp);
    nc_varput(fnm_out,'Vs',Vs);
    nc_varput(fnm_out,'rho',Dp);
    nc_varput(fnm_out,'Qs',Qs);
    
    disp('finished creating')
end


