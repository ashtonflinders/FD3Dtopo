%close all
clear all

% calls a mini script shared by other codes to define the x and y coordinates
define_xy_coords

fmt_out='%f\n';
% x
fid=fopen('SeisGrid.coordx.dat','w');
fprintf(fid,'%i\n',nx);
fprintf(fid,fmt_out,x);
fclose(fid);
% y
fid=fopen('SeisGrid.coordy.dat','w');
fprintf(fid,'%i\n',ny);
fprintf(fid,fmt_out,y);
fclose(fid);

