clear all

MFILE_ROOT='./';
path([MFILE_ROOT '/fun-spool'],path);
path([MFILE_ROOT '/saclab'],path);
path([MFILE_ROOT '/fileexchange'],path);

%pnm=pwd; indx=strfind(pnm,'/'); event_name=pnm(indx(end)+1:end);
event_name='seismo.3000x2400_etopo2_sed1_crust2_ak135_USGS5x5smooth';

%fnm_nc='2005.01.10.18.47.synthetic.interp2.cubic.nc';
%fnm_nc='2005.01.10.18.47.synthetic.interp2.linear.nc';
%fnm_nc='2005.01.10.18.47.synthetic.recv.nearest.nc';

%fnm_nc=[event_name '.synthetic.interp2.linear.nc'];
fnm_nc=[event_name '.synthetic.nearest.nc'];

pnm_sac = [event_name '_sac']; mkdir(pnm_sac);

dinfo = nc_getdiminfo(fnm_nc,'number_of_station'); NSTAT=dinfo.Length;

for n=1:NSTAT
    t=nc_varget(fnm_nc,'time');
    Vx=nc_varget(fnm_nc,'Vx',[n-1,0],[1,-1]);
    Vy=nc_varget(fnm_nc,'Vy',[n-1,0],[1,-1]);
    Vz=nc_varget(fnm_nc,'Vz',[n-1,0],[1,-1]);
    lat=nc_varget(fnm_nc,'latitude',[n-1],[1]);
    lon=nc_varget(fnm_nc,'longitude',[n-1],[1]);
    snm=nc_varget(fnm_nc,'station_name',[n-1,0],[1,-1]);

    t=[0 t'];Vx=[0 Vx];Vy=[0 Vy];Vz=[0 Vz];

    Sx=bsac(t,Vx); Sy=bsac(t,Vy); Sz=bsac(t,Vz);
    ch(Sx,'STLA',lat,'STLO',lon);
    ch(Sy,'STLA',lat,'STLO',lon);
    ch(Sz,'STLA',lat,'STLO',lon);

    wsac([pnm_sac '/' event_name '.' strtrim(snm) '.SHN.SAC'],Sx);
    wsac([pnm_sac '/' event_name '.' strtrim(snm) '.SHE.SAC'],Sy);
    wsac([pnm_sac '/' event_name '.' strtrim(snm) '.SHZ.SAC'],Sz);
end
