function mapstruct = computePVxy_geo(mapstruct)

% Compute the geostrophic velocity field and add it to the map structure.
% Input parameters:
    % 'mapstruct':  'map' or 'avmap'
%
%
% D. Rudnick, June 3, 2019

% Constants
rho0 = 1027;    % reference density, kg/m^3

[nz,~,~,nsecs] = size(mapstruct.sigmastable);

dx=1000*(mapstruct.x(2)-mapstruct.x(1));    
dy=1000*(mapstruct.y(2)-mapstruct.y(1));
dz = mapstruct.depth(1)-mapstruct.depth(2);         % vertical spacing, m. depth is positive down

f=repmat(sw_f(mapstruct.lat),1,1,nsecs,nz);
f=permute(f,[4 1 2 3]);
mapstruct.f=f;

[dsigmadx,dsigmadz,dsigmady] = gradient(mapstruct.sigmastable,dx,dz,dy,1);
[dvdx,dvdz,dvdy] = gradient(mapstruct.geov,dx,dz,dy,1);
[dudx,dudz,dudy] = gradient(mapstruct.geou,dx,dz,dy,1);

mapstruct.pvgeo=-((f+dvdx-dudy).*dsigmadz-dvdz.*dsigmadx+dudz.*dsigmady)./rho0;
mapstruct.pvvertgeo=-((f+dvdx-dudy).*dsigmadz)./rho0;
mapstruct.pvhorizgeo=(dvdz.*dsigmadx-dudz.*dsigmady)./rho0;
mapstruct.rogeo=(dvdx-dudy)./f;
mapstruct.divgeo=(dudx+dvdy)./f;
