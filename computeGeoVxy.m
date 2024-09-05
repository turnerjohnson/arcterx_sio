function map = computeGeoVxy(map)

% Compute the geostrophic velocity field.
% Input parameters:
    % 'map':  'map' or 'avmap'
%
%
% D. Rudnick, 6 Apr 2020, xy version - variational method








% Constants
g = 9.81;       % gravity, m/s^2
rho0 = 1027;    % reference density, kg/m^3
nzmax=40;       % This is the standard value: max depth levels

f=sw_f(mean(map.lat(:))); %Coriolis

[nz,nx,ny,ntime] = size(map.rho); %density field dims

dx=(map.x(2)-map.x(1));  %x grid spacing in km  
dy=(map.y(2)-map.y(1));  %y grid spacing in km
dz=-10;   %z grid spacing, negative because depth is positive downward

R=(g./rho0./f).*cumtrapz(map.sigmastable,1)*dz/1000; % Rossby radius calculation
Rwig=squeeze(mean(R(1:nzmax,:,:,:),1)); %depth averages Rossby Radius
zeta=1000*map.ro.*map.f; %relative vorticity in units of km
zetawig=squeeze(mean(zeta(1:nzmax,:,:,:),1)); % depth avg. relative vort.
uwig=squeeze(mean(map.udopacross(1:nzmax,:,:,:),1)); % depth-avg across-track velocity
vwig=squeeze(mean(map.udopalong(1:nzmax,:,:,:),1)); % depth avg along-track velocity
[~,~,Dsq]=laplacian([nx ny]); % Laplacian operator
Dsq=-Dsq/dx/dy;             % scaled Laplacian operator

% boundary conditions: y-dir, N/S
indsouth=sub2ind([nx ny],(2:nx-1)',ones(nx-2,1)); % south boundary indices
vsouth=ones(nx-2,1)/dy;% south boundary values
indsouthin=sub2ind([nx ny],(2:nx-1)',2*ones(nx-2,1)); % inner south boundary indices
vsouthin=-vsouth; % inner south boundary values
indnorth=sub2ind([nx ny],(2:nx-1)',ny*ones(nx-2,1)); % north boundary indices
vnorth=vsouth; % north boundary values
indnorthin=sub2ind([nx ny],(2:nx-1)',(ny-1)*ones(nx-2,1)); % inner north boundary indices
vnorthin=vsouthin; % inner north boundary values.

% boundary conditions: x-direction E/W
indwest = sub2ind([nx ny], ones(ny,1), (1:ny)'); % west boundary indices
vwest = ones(ny,1) / dx; % west bound values
indwestin = sub2ind([nx ny], 2*ones(ny,1), (1:ny)'); %inner west boundary indices
vwestin = -vwest; % inner west boundary values
indeast = sub2ind([nx ny], nx*ones(ny,1), (1:ny)'); % east boundary indices
veast = vwest; % east bound values
indeastin = sub2ind([nx ny], (nx-1)*ones(ny,1), (1:ny)'); % inner east indices
veastin = vwestin; % inner east values

%apply BCs to to Dsq (Laplacian Operator)
indcol=[indsouth; indsouthin; indnorth; indnorthin; indwest; indwestin; indeast; indeastin];
indrow=[indsouth; indsouth; indnorth; indnorth; indwest; indwest; indeast; indeast];
indDsq=sub2ind([nx*ny nx*ny],indrow,indcol);
vs=[vsouth; vsouthin; vnorth; vnorthin; vwest; vwestin; veast; veastin];
Dsq(indrow,:)=0;
Dsq(indDsq)=vs;
drow=[indsouth; indnorth; indwest; indeast];

% modify Laplacian: a weird little try to change the magnitude
Dsq(2,:)=0;
Dsq(2,2)=1;

map.psi=NaN(nx,ny,ntime); %initialize streamfunction

% solve for streamfunction
for n=1:ntime
   Rwign=squeeze(Rwig(:,:,n)); % extract Rossby radius for time n
   DsqRwign=Dsq*Rwign(:); % apply Laplacian to Rossby radius
   zetawign=squeeze(zetawig(:,:,n)); % extract vorticity to time n
   zetawign=zetawign(:); 
   d=DsqRwign+zetawign; % RHS
   % apply BCs
   d(drow)=DsqRwign(drow)+[uwig(2:nx-1,1,n); -uwig(2:nx-1,ny,n); -vwig(1,:,n)'; vwig(nx,:,n)']; 
   % the weird try again
   d(2)=0; % weird try at specific condition?
   % weird try
   psin=Dsq\d; % solving for streamfunction!
   map.psi(:,:,n)=reshape(psin,nx,ny); %store
end

% Calculate geostrophic velocities
[psiy,psix]=gradient(map.psi,dy,dx,1); % gradients of streamfunction
psix=repmat(psix,1,1,1,nz); % replicate matrix
psix=permute(psix,[4 1 2 3]);
psiy=repmat(psiy,1,1,1,nz);
psiy=permute(psiy,[4 1 2 3]);

[Rx,~,Ry]=gradient(R,dx,dz,dy,1); % gradients of Rossby radius

map.geou=-psiy+Ry; % geostrophic u/across-track velocity
map.geov=psix-Rx; % geostrophic v/along-track velocity

% [rhox,~,rhoy] = gradient(map.sigmastable,1000*dx,1000*dz,dy,1);      % density gradients to km also
% 
% dvdz = -(g/rho0).*rhox./f;                     % geostrophic shear
% vshear = cumtrapz(dvdz,1)*dz;                     % integral of geostrophic shear
% map.geovf = vshear + repmat(map.ualong-nanmean(vshear),nz,1,1,1);
% 
% dudz = (g/rho0).*rhoy./f;                     % geostrophic shear
% ushear = cumtrapz(dudz,1)*dz;                     % integral of geostrophic shear
% map.geouf = ushear + repmat(map.uacross-nanmean(ushear),nz,1,1,1);
