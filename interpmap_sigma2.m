function sigint=interpmap_sigma2(mapstruct,sigsig)
%function to interpolate to surfaces from mapped data.
%

sigma=mapstruct.sigmastable;
p=mapstruct.depth;

nsig=length(sigsig);
[nz,nx,ny,ntime]=size(sigma);

sigint.sigma=sigsig(:);
sigint.lat=mapstruct.lat;
sigint.lon=mapstruct.lon;
sigint.x=mapstruct.x;
sigint.y=mapstruct.y;
sigint.time=mapstruct.time;

vars=fieldnames(mapstruct);
k=0;

sigma=reshape(sigma,nz,[]);
ntot=nx*ny*ntime;
sigint.depth=nan(nsig,ntot);

for kk=1:length(vars)
   if ~strncmp('sigma',vars{kk},5) && ndims(mapstruct.(vars{kk})) == 4 && size(mapstruct.(vars{kk}),1) ~= 1
      sigint.(vars{kk})=nan(nsig,ntot);
      k=k+1;
      vsig{k}=vars{kk};
   end
end



for ii=1:ntot
   for jj=1:nsig
      j1=max(find(sigma(:,ii) <= sigsig(jj)));
      if ~isempty(j1) && j1 ~= nz
         coef=(sigsig(jj)-sigma(j1,ii))/(sigma(j1+1,ii)-sigma(j1,ii));
         for k=1:length(vsig)
            vtmp=reshape(mapstruct.(vsig{k}),nz,[]);
            sigint.(vsig{k})(jj,ii)=vtmp(j1,ii)+coef*(vtmp(j1+1,ii)-vtmp(j1,ii));
         end
         sigint.depth(jj,ii)=p(j1)+coef*(p(j1+1)-p(j1));
      end
   end
end

sigint.depth=reshape(sigint.depth,nsig,nx,ny,ntime);
for k=1:length(vsig)
   sigint.(vsig{k})=reshape(sigint.(vsig{k}),nsig,nx,ny,ntime);
end


