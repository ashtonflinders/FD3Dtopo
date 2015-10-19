% define the x and y coordinates to be used by other matlab codes

minlon=-112.63; maxlon=-111.76; minlat=35.84;  maxlat=36.52;
% center of the xy coordinates in latitude and longitude
lat0=(minlat+maxlat)/2; lon0=(minlon+maxlon)/2; 
x0=0; y0=0; alpha=90;
% 

delta_x=0.5; delta_y=0.5; % 500 m horizontal grid spacing
% find nx and ny that can be divided by an integer (# of blks in x and y)
ny=round((maxlat-minlat)*111/delta_x/4)*4
nx=round((maxlon-minlon)*111*cos(lat0*pi/180)/delta_y/4)*4

fid = fopen('model_par.txt','wt');
fprintf(fid,'%s %7.2f\n', 'minlon = ', minlon);
fprintf(fid,'%s %7.2f\n', 'maxlon = ', maxlon);
fprintf(fid,'%s %7.2f\n', 'minlat = ', minlat);
fprintf(fid,'%s %7.2f\n', 'maxlat = ', maxlat);
fprintf(fid,'%s %8.3f\n','lon0 = ', lon0);
fprintf(fid,'%s %8.3f\n','lat0 = ', lat0);
fprintf(fid,'%s %4.2f\n','dx = ',delta_x);
fprintf(fid,'%s %4.2f\n','dy = ',delta_y);
fprintf(fid,'%s %5.0f\n','nx = ',nx);
fprintf(fid,'%s %5.0f\n','ny = ',ny);

x=[-nx/2:1:nx/2-1]*delta_x*1e3;
y=[-ny/2:1:ny/2-1]*delta_y*1e3;

if nx ~= length(x) | ny ~= length(y)
display('warning: length of coordinate vector is not equal to the desired length')
end

