function [lat2, lon2] = geo2gcp(lat1,lon1,d,theta)
% given the lat and lon of a point on earth in degrees, distance and azimuth 
% of the path, calculate lat and lon of the end point
R=6371; %earth mean radius
lat1=lat1*pi/180; lon1=lon1*pi/180;
theta=theta*pi/180;
lat2=asin(sin(lat1)*cos(d/R)+cos(lat1)*sin(d/R)*cos(theta));
lon2=lon1+atan2(sin(theta)*sin(d/R)*cos(lat1),cos(d/R)-sin(lat1)*sin(lat2));
lat2=lat2*180/pi;
lon2=lon2*180/pi;

