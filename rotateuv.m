function bindata = rotateuv(bindata,lon0,lat0,lon1,lat1)

% Calculates along and across shore velocity given the line.
% bindata = calcualongacross(bindata,line)

% adapted from calcdistfromshore
% D. Rudnick, 19 Oct 2012
% D. Rudnick, 16 May 2017 - added support for udop

% Set parameters
deg2min = 60;
deg2rad = pi/180;

% Calculate angle of new coordinate system relative to east
dyy = (lat1-lat0);
dxx = cos(1/2*(lat1+lat0)*deg2rad)*(lon1-lon0);
theta = atan2(dyy,dxx);

% Calculate along and across shore velocity
w=(bindata.u+1i*bindata.v)*exp(-1i*theta);
bindata.ualong=imag(w);
bindata.uacross=real(w);

%udop, if exists
if isfield(bindata,'udop')
   w=(bindata.udop+1i*bindata.vdop)*exp(-1i*theta);
   bindata.udopalong=imag(w);
   bindata.udopacross=real(w);
end

