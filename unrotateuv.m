function bindata = unrotateuv(bindata,lon0,lat0,lon1,lat1)

% bindata = calcualongacross(bindata,line)
% Calculates along and across shore velocity given the line.

% adapted from calcdistfromshore
% D. Rudnick, 19 Oct 2012
% D. Rudnick, 16 May 2017 - added support for udop

% Set parameters
deg2min = 60;
deg2rad = pi/180;

% Calculate angle of new coordinate system relative to east
dyy = (lat1-lat0);
dxx = cos(1/2*(lat1+lat0)*deg2rad)*(lon1-lon0);
theta = -atan2(dyy,dxx);

% Calculate along and across shore velocity
if isfield(bindata,'uacross')
w=(bindata.uacross+1i*bindata.ualong)*exp(-1i*theta);
bindata.u=real(w);
bindata.v=imag(w);
end

%udop, if exists
if isfield(bindata,'udopacross')
   w=(bindata.udopacross+1i*bindata.udopalong)*exp(-1i*theta);
   bindata.udop=real(w);
   bindata.vdop=imag(w);
end

%geov if it exists
if isfield(bindata,'geou')
   w=(bindata.geou+1i*bindata.geov)*exp(-1i*theta);
   bindata.geoull=real(w);
   bindata.geovll=imag(w);
end
