function mapout = computeageouv(map)

% Compute ageostrophic horiztonal velocity given mapw including w.
%
% D. Rudnick, 22 Jun 2020 

% Constants
g = 9.81;       % gravity, m/s^2
rho0 = 1027;    % reference density, kg/m^3

f=sw_f(mean(map.lat(:))); % Coriolis

mapout = map; %initialize return array

% grid spacings (meters)
dx=(map.x(2)-map.x(1))*1000;
dy=(map.y(2)-map.y(1))*1000;
dz=-10;   %negative because depth is positive downward

% gradients of sigma, time?
[dsigmadx,dsigmadz,dsigmady,~] = gradient(map.sigmastable,dx,dz,dy,86400);
N2 = -g/rho0/(f*f)*dsigmadz; % Brunt Vaisala
[dN2dx,~,dN2dy] = gradient(N2,dx,dz,dy,1); % gradients of N2

% gradients of geostorphic velocities 
[dugdx,~,dugdy,~]=gradient(map.geou,dx,dz,dy,86400);
[dvgdx,~,dvgdy,~]=gradient(map.geov,dx,dz,dy,86400);

% gradient of calculated vertical velocity
[dwdx,~,dwdy,~]=gradient(map.w,dx,dz,dy,86400);

% PV vectors
Qx=2*g/rho0/(f*f)*(dugdx.*dsigmadx+dvgdx.*dsigmady);
Qy=2*g/rho0/(f*f)*(dugdy.*dsigmadx+dvgdy.*dsigmady);

dN2wdx = dN2dx.*map.w + N2.*dwdx;
dN2wdy = dN2dy.*map.w + N2.*dwdy;

% integrate to get ageostrophic velocities
mapout.ua = -flipud(cumtrapz(flipud(dN2wdx-Qx),1)*dz);  %negative sign and flipud to integrate from bottom up
mapout.va = -flipud(cumtrapz(flipud(dN2wdy-Qy),1)*dz);
