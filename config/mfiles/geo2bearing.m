function [b] = geo2bearing(lat1,lon1,lat2,lon2)
% given the lat and long of two points on earth in degrees, output lat and lon of gcp path
% http://www.movable-type.co.uk/scripts/latlong.html
R=6371; %earth mean radius
lat1=lat1*pi/180; lon1=lon1*pi/180;
lat2=lat2*pi/180; lon2=lon2*pi/180;
%d=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))*R;
dlon=lon2-lon1;
b=atan2(sin(dlon)*cos(lat2),cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon));
b=b*180/pi;


