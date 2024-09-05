function map=computeStableDensity(map)
% function computeStableDensity(map) finds the profiles closest to the
% mapped potential density profiles and adds the result to map as
% sigmastable
%
% D. Rudnick, 18 March 2020

[nz,nx,ny,ntime]=size(map.sigma);

sigma=reshape(map.sigma,nz,[]);

A=diag(ones(nz,1),0)-diag(ones(nz-1,1),1);
A(nz,:)=[];
b=zeros(nz-1,1);
C=eye(nz);
options=optimset('display','off');
numfix=0;
for n=1:size(sigma,2)
   if any(A*sigma(:,n) > 0) % if unstable (density not monotonic?)
      numfix=numfix+1;      % use least squares to find closest stable profile.
      sigma(:,n)=lsqlin(C,sigma(:,n),A,b,[],[],[],[],[],options);
   end
end
disp(numfix)
map.sigmastable=reshape(sigma,nz,nx,ny,ntime);
