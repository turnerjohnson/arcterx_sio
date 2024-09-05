function sf=structurefunction(ctd,iminemax,iz,dtimemax,numbins)
% sf=structurefunction(ctd,iz) calculates structure function given ctd
% structure and the depth index iz and maxdtime in hours
%

% get only my data
ii=ctd.missid <= iminemax;
x=ctd.dist(ii);
y=ctd.offset(ii);
time=ctd.time(ii);
u=ctd.udopacross(:,ii);
v=ctd.udopalong(:,ii);

% calculate deltas
dx=x-x';
dy=y-y';
dtime=time-time';
du=u(iz,:)'-u(iz,:);
dv=v(iz,:)'-v(iz,:);
dr=sqrt(dx.^2+dy.^2);

% rotate the velocities
theta = atan2(dy,dx);
w=(du+1i*dv).*exp(-1i*theta);
dul=real(w);
dut=imag(w);

%plot locaations of profiles
subplot(1,3,1);
cla;
hold on;
for nn=1:iminemax
   plot(ctd.lon(ctd.missid == nn),ctd.lat(ctd.missid == nn),'.');
end
set(gca,'dataaspectratio',[1 cosd(mean(ctd.lat(ii))) 1]);
hold off;
xlabel('Longitude');
ylabel('Latitude');

% plot histogram of separations
subplot(1,3,2);
jj=dtime > 0 & dtime < dtimemax*3600;
bin=linspace(log10(min(dr(jj))),log10(max(dr(jj))),numbins);
binr=10.^bin;
histogram(dr(jj),binr)
set(gca,'xscale','log')
xlabel('\Deltar (km)');
ylabel('Count');
title([num2str(ctd.depth(iz)) ' m'])

% calculate and plot structure function
subplot(1,3,3);
kk=discretize(dr(jj),binr);
binu3=zeros(1,length(binr)-1);
binu2=binu3;
binu=binu3;
duljj=dul(jj);
for n=1:length(binu)
   binu3(n)=nanmean(duljj(kk == n).^3);
   binu2(n)=nanmean(duljj(kk == n).^2);
   binu(n)=nanmean(duljj(kk == n));
end
plot(runmean(binr',2),binu3,'-x');
set(gca,'xscale','log');
grid on;
xlabel('\Deltar (km)');
ylabel('<\Deltau_l^3> (m^3/s^3)')
title(['\Deltat < ' num2str(dtimemax) ' h'])

sf.binu3=binu3;
sf.binr=binr;
sf.binu2=binu2;
sf.binu=binu;
