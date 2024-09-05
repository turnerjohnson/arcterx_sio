function mapstruct = addderivedvarsxy(mapstruct,ctd)

% Computes derived fields and adds them to the map
% structure (mapstruct).
%
%
% K.Zaba: Jun14,2016
% D. Rudnick, July 28, 2016
% D. Rudnick, 10 May 2019 - Calypso version
% D. Rudnick, 27 Feb 2019 - Calypso xy

t = mapstruct.t;
s = mapstruct.s;
[nz,nx,ny,nt] = size(t);

% Initialize derived fields
[xx,yy]=meshgrid(mapstruct.x,mapstruct.y);
[mapstruct.lon,mapstruct.lat]=do2ll(xx',yy',ctd.x0,ctd.y0,ctd.x1,ctd.y1);
p=sw_pres(mapstruct.depth*ones(1,nx*ny),ones(nz,1)*mapstruct.lat(:)');
mapstruct.theta = nan(nz,nx,ny,nt);
mapstruct.rho   = nan(nz,nx,ny,nt);
mapstruct.sigma = nan(nz,nx,ny,nt);

% Loop through sections, compute derived fields: theta, rho, sigma
for iit=1:nt
   ss=s(:,:,:,iit);
   ss=ss(:);
   tt=t(:,:,:,iit);
   tt=tt(:);
   theta=sw_ptmp(ss,tt,p(:),0);
   mapstruct.theta(:,:,:,iit)=reshape(theta,nz,nx,ny);
   rho=sw_dens(ss,tt,p(:))-1000;
   mapstruct.rho(:,:,:,iit)=reshape(rho,nz,nx,ny);
   sigma=sw_pden(ss,tt,p(:),0)-1000;
   mapstruct.sigma(:,:,:,iit)=reshape(sigma,nz,nx,ny);
end
