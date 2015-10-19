function [latm, lonm] = geo2gcp(lat1,lon1,lat2,lon2)
% given the lat and long of two points on earth in degrees, output lat and lon of gcp path
R=6371; %earth mean radius
lat1=lat1*pi/180; lon1=lon1*pi/180;
lat2=lat2*pi/180; lon2=lon2*pi/180;
%d=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))*R;
dlon=lon2-lon1;
Bx=cos(lat2)*cos(dlon);
By=cos(lat2)*sin(dlon);
latm=atan2(sin(lat1)+sin(lat2),sqrt((cos(lat1)+Bx)^2 + By*By));
lonm=lon1+atan2(By,cos(lat1)+Bx);
latm=latm*180/pi;
lonm=lonm*180/pi;

