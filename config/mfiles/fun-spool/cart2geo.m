function [lat,lon]=cart2geo(x,y,x0,y0,lat0,lon0,alpha)
%
% cart2geo: Convert Cartesian coordinate to geographyic coordinate
%
% Usage: [lat,lon]=cart2geo(x,y,x0,y0,lat0,lon0,alpha)
%

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Added help information, but uncomplete.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date: 2008-04-27 17:31:28 -0400 (Sun, 27 Apr 2008) $
% $Revision: 469 $
% $LastChangedBy: zhangw $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=alpha/180*pi;

D=111319.5;

lat=lat0+( (x-x0)*cos(a)+(y-y0)*sin(a) )/D;
lon=lon0+( (x-x0)*sin(a)-(y-y0)*cos(a) )/D./cos(lat/180*pi);
