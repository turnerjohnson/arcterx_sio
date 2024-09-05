function [lon,lat] = do2ll(dist,offset,lon0,lat0,lon1,lat1)

% Calculates lat,lon given dist, offset start and end of line lon0, lat0,
% lon1, lat1
%
% D. Rudnick, 28 May 2018

% Set parameters
nm2km = 1.852;          %nautical miles to kilometers
deg2min = 60;           %degrees to minutes
deg2rad = pi/180;       %degrees to radians
deg2km=deg2min*nm2km;   %degrees to km

% Calculate angle of new coordinate system relative to east
dyy = (lat1-lat0);
dxx = cos(1/2*(lat1+lat0)*deg2rad)*(lon1-lon0);
theta = atan2(dyy,dxx);

% Calculate dx, dy of lon, lat in new coordinate system by rotating
zhat=dist+1i*offset;
z=zhat*exp(1i*theta);

dx=real(z);
dy=imag(z);

% Calculate lon, lat 
lat=dy/deg2km+lat0;
lon=dx/deg2km/cos(1/2*(lat1+lat0)*deg2rad)+lon0;
