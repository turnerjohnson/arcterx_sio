function mapstruct = computePVxy(mapstruct)

% Compute the geostrophic velocity field and add it to the map structure.
% %% potential vorticity??????????????????????????
% Input parameters:
    % 'mapstruct':  'map' or 'avmap'
%
%
% D. Rudnick, June 3, 2019

% Constants
% rho0 = 1027;    % reference density, kg/m^3

[nz,~,~,nsecs] = size(mapstruct.rho);

dx=1000*(mapstruct.x(2)-mapstruct.x(1));    
dy=1000*(mapstruct.y(2)-mapstruct.y(1));
dz = mapstruct.depth(1)-mapstruct.depth(2);         % vertical spacing, m. depth is positive down

f=repmat(sw_f(mapstruct.lat),1,1,nsecs,nz);
f=permute(f,[4 1 2 3]);
mapstruct.f=f;

[dsigmadx,dsigmadz,dsigmady] = gradient(mapstruct.sigma,dx,dz,dy,1);
[dvdx,dvdz,dvdy] = gradient(mapstruct.udopalong,dx,dz,dy,1);
[dudx,dudz,dudy] = gradient(mapstruct.udopacross,dx,dz,dy,1);

mapstruct.pv = -((f+dvdx-dudy).*dsigmadz-dvdz.*dsigmadx+dudz.*dsigmady)./(mapstruct.rho+1000);
mapstruct.pvvert = -((f+dvdx-dudy).*dsigmadz)./(mapstruct.rho+1000);
mapstruct.pvhoriz = (dvdz.*dsigmadx-dudz.*dsigmady)./(mapstruct.rho+1000);
mapstruct.ro = (dvdx-dudy)./f;
mapstruct.div = (dudx+dvdy)./f;
