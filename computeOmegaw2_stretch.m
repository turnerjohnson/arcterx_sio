function mapout = computeOmegaw2_stretch(map)

% Compute omega vertical velocity.
%
% D. Rudnick, 6 Apr 2020, xy version - variational method
% D. Rudnick, 4 Sep 2020, make the sparse matrix in one fell swoop, way
% faster
% D. Rudnick 8 Apr 2021, stretch vertically

% Constants
g = 9.81;       % gravity, m/s^2
rho0 = 1027;    % reference density, kg/m^3
nstretch=100;    % vertical stretch factor

mapout.x=map.x;
mapout.y=map.y;
mapout.depth=[0; map.depth]; % add surface to depth array
mapout.time=map.time;
mapout.lat=map.lat;
mapout.lon=map.lon;

f=sw_f(mean(map.lat(:))); % Coriolis

[~,nx,ny,ntime] = size(map.sigmastable);  % dimensions from sigmastable, nz defined after stretch
nzout=length(mapout.depth); % number of depth levels to return

dx=(map.x(2)-map.x(1))*1000; % calculate grid spacing
dy=(map.y(2)-map.y(1))*1000; 
dz=-10;   %negative because depth is positive downward

sigmastable=cat(1,map.sigmastable(1,:,:,:),map.sigmastable); % add surface layer to sigmastable
geou=cat(1,map.geou(1,:,:,:),map.geou); % add surface layer to geo u
geov=cat(1,map.geov(1,:,:,:),map.geov); % add surf layer to geo v

mapout.sigmastable=sigmastable;
mapout.geou=geou;
mapout.geov=geov;

[dsigmadx,dsigmadz,dsigmady]=gradient(sigmastable,dx,dz,dy,1); % calculate sigma gradients
d2sigmadx2=secondderivative(sigmastable,2,dx); % 2nd deriv in x, Dan function
d2sigmady2=secondderivative(sigmastable,3,dy); % 2nd deriv in y
[~,~,d2sigmadxdy]=gradient(dsigmadx,dx,dz,dy,1); % d squared /dxdy

[dugdx,~,dugdy]=gradient(geou,dx,dz,dy,1); % calculate gradients of geo u
d2ugdx2=secondderivative(geou,2,dx);       % 2nd x deriv
d2ugdy2=secondderivative(geou,3,dy);       % 2nd y deriv

[dvgdx,~,dvgdy]=gradient(geov,dx,dz,dy,1); % gradients in geo v
d2vgdx2=secondderivative(geov,2,dx);       % 2nd x deriv
d2vgdy2=secondderivative(geov,3,dy);       % 2nd y deriv

N2=-g/rho0/(f*f)*dsigmadz;  % calculate buoyancy frequency
N2(end+1:end+nstretch,:,:,:)=0;  % add stretch layer with N2=0
% N2(end+1:end+nstretch,:,:,:)=repmat(N2(end-5,:,:,:),[nstretch 1 1 1]);  % here is the stretch
[dN2dx,~,dN2dy] = gradient(N2,dx,dz,dy,1); % calculate gradients of N2
d2N2dx2 = secondderivative(N2,2,dx); % 2nd x deriv N2
d2N2dy2 = secondderivative(N2,3,dy); % 2nd y deriv N2

% forcing term
force = 2*g/rho0/(f*f)*((d2ugdx2+d2ugdy2).*dsigmadx + (d2vgdx2+d2vgdy2).*dsigmady + ...
   (dugdy+dvgdx).*d2sigmadxdy + dugdx.*d2sigmadx2 +dvgdy.*d2sigmady2);
force(end+1:end+nstretch,:,:,:) = 0;  % stretched layers have 0 forcing
% set forcing to 0 at boundaries
force(1,:,:,:)=0;
force(end,:,:,:)=0;
force(:,1,:,:)=0;
force(:,end,:,:)=0;
force(:,:,1,:)=0;
force(:,:,end,:)=0;
force = reshape(force,[],ntime);

nzp1 = size(N2,1); % vertical dimension size, to account for the stretch
mapout.w = NaN(nzp1,nx,ny,ntime); %initialize array to return
ntot = nzp1*nx*ny; % total number grid points

% set up grid indices for matrix creation
sz=[nzp1 nx ny];
nz=nzp1-1;
nzm1=nz-1;
onesnzm1=ones(nzm1,1);
iz0=(2:nz)';
izp1=iz0+1;
izm1=iz0-1;
dz2=dz*dz;
dx2=dx*dx;
dy2=dy*dy;

vzp1=1/dz2*onesnzm1;
vzm1=vzp1;

% initialize arrays for matrix creation
for itime=1:ntime 
% for itime=29:29
   kk=zeros(7*nzm1*(nx-2)*(ny-2)+2*nx*ny+2*nx*nz+2*ny*nz,1);
   jj=kk;
   vv=kk;
   ilast=0;

   %Dirichlet for w all around
   for iy=1:ny
      for ix=1:nx
         kk1=sub2ind(sz,1,ix,iy);
         kk2=sub2ind(sz,nzp1,ix,iy);
         nval=2;
         kk(ilast+1:ilast+nval)=[kk1; kk2];
         jj(ilast+1:ilast+nval)=[kk1; kk2];
         vv(ilast+1:ilast+nval)=[1; 1];
         ilast=ilast+nval;
      end
   end
   for iy=1:ny
      for iz=2:nzp1
         kk1=sub2ind(sz,iz,1,iy);
         kk2=sub2ind(sz,iz,nx,iy);
         nval=2;
         kk(ilast+1:ilast+nval)=[kk1; kk2];
         jj(ilast+1:ilast+nval)=[kk1; kk2];
         vv(ilast+1:ilast+nval)=[1; 1];
         ilast=ilast+nval;
      end
   end
   for iz=2:nzp1
      for ix=2:nx-1
         kk1=sub2ind(sz,iz,ix,1);
         kk2=sub2ind(sz,iz,ix,ny);
         nval=2;
         kk(ilast+1:ilast+nval)=[kk1; kk2];
         jj(ilast+1:ilast+nval)=[kk1; kk2];
         vv(ilast+1:ilast+nval)=[1; 1];
         ilast=ilast+nval;
      end
   end

   % interior
   for iy=2:ny-1
      iy0=iy*onesnzm1;
      iyp1=(iy+1)*onesnzm1;
      iym1=(iy-1)*onesnzm1;
      for ix=2:nx-1
         ix0=ix*onesnzm1;
         ixp1=(ix+1)*onesnzm1;
         ixm1=(ix-1)*onesnzm1;
         
         j0=sub2ind(sz,iz0,ix0,iy0);
         jzp1=sub2ind(sz,izp1,ix0,iy0);
         jzm1=sub2ind(sz,izm1,ix0,iy0);
         
         jxp1=sub2ind(sz,iz0,ixp1,iy0);
         jxm1=sub2ind(sz,iz0,ixm1,iy0);
         
         jyp1=sub2ind(sz,iz0,ix0,iyp1);
         jym1=sub2ind(sz,iz0,ix0,iym1);
         
         k=repmat(j0,7,1);
         
         v0=-2/dz2+d2N2dx2(iz0,ix,iy,itime)+d2N2dy2(iz0,ix,iy,itime)-2*(1/dx2+1/dy2)*N2(iz0,ix,iy,itime);
         vxp1=N2(iz0,ix,iy,itime)/dx2+dN2dx(iz0,ix,iy,itime)/dx;
         vxm1=N2(iz0,ix,iy,itime)/dx2-dN2dx(iz0,ix,iy,itime)/dx;
         vyp1=N2(iz0,ix,iy,itime)/dy2+dN2dy(iz0,ix,iy,itime)/dy;
         vym1=N2(iz0,ix,iy,itime)/dy2-dN2dy(iz0,ix,iy,itime)/dy;
         
         nval=7*nzm1;
         kk(ilast+1:ilast+nval)=k;
         jj(ilast+1:ilast+nval)=[j0; jzp1; jzm1; jxp1; jxm1; jyp1; jym1];
         vv(ilast+1:ilast+nval)=[v0; vzp1; vzm1; vxp1; vxm1; vyp1; vym1];
         ilast=ilast+nval;
      end
   end
   % create sparse matrix
   S=sparse(kk(1:ilast),jj(1:ilast),vv(1:ilast),ntot,ntot);
   wn=S\force(:,itime); % Solve for w
   mapout.w(:,:,:,itime)=reshape(wn,nzp1,nx,ny);
end

mapout.w(nzout+1:end,:,:,:)=[]; % remove extra layers?
