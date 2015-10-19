clear all
addpath(genpath('../../../config/mfiles/'))
% load x,y,z
fnm_topo='model.hiresDEM.nc'
x=nc_varget(fnm_topo,'x'); nx=length(x);
y=nc_varget(fnm_topo,'y'); ny=length(y);
z0=nc_varget(fnm_topo,'topo');


z=[-10, -30, -60]*1e3;
g=[20 20 15];

ng=sum(g)
disp(['number of cells = ' num2str(ng)]);
disp(['number of vertical grid points = ' num2str(ng+1)]); 

nlayer=length(z);

fmt_vmap='%15.7f %15.7f';
for n=1:nlayer+1; fmt_vmap=[fmt_vmap ' %16.7f']; end
fmt_vmap=[fmt_vmap '\n']; %output format: x,y,z(1:nlayer)

fid=fopen('SeisGrid.vmap.dat','w');

fprintf(fid,'%i %i %i\n',nx,ny,nlayer+1);
for j=1:ny
    for i=1:nx
       fprintf(fid, fmt_vmap, x(i), y(j), ...
        [z0(j,i) z0(j,i)+z(1) z(2:nlayer)]);
    end
end

