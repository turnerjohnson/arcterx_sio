function map = computeGeoVxy_reflevel(map,nreflev)

% Compute the geostrophic velocity field.
% Input parameters:
    % 'map':  'map' or 'avmap'
%
%
% D. Rudnick, 6 Apr 2020, xy version - variational method
% D. Rudnick, 13 Jun 2022, xy version - reference level


% Constants
g = 9.81;       % gravity, m/s^2
rho0 = 1027;    % reference density, kg/m^3

f=sw_f(mean(map.lat(:)));

dx=(map.x(2)-map.x(1));  %in km  
dy=(map.y(2)-map.y(1));  %in km
dz=-10;   %negative because depth is positive downward

R=(g./rho0./f).*cumtrapz(map.sigmastable,1)*dz/1000; %another fix for km
Rref=R-R(nreflev,:,:,:);

[Rx,~,Ry]=gradient(Rref,dx,dz,dy,1);

map.geou=Ry;
map.geov=-Rx;
