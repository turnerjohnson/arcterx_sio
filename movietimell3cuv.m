function movietimell3cuv(map,var1str,whichuv,iz,cont1,nskip,ctd,dtime)

ntend=length(map.time);

map=unrotateuv(map,ctd.x0,ctd.y0,ctd.x1,ctd.y1);

if isfield(map,'sigmastable')
   cont2=20:0.1:30;
   cont2h=20:0.2:30;
   var2='sigmastable';
   z = map.depth;
   units = ' m';
else
   cont2=0:10:1000;
   cont2h=0:20:1000;
   var2='depth';
   z = map.sigma;
   units = ' kg/m^3';
end

dtimeut=dtime*86400;

v = VideoWriter('ARCTERX','MPEG-4');
v.FrameRate=5;
open(v);

for n=1:ntend
   for nplot=1:3
      subplot(1,3,nplot);
      if any(any(nanmean(map.(var1str{nplot})(iz(nplot),:,:,n),3)))
         pxyll_xyuv(map,var1str{nplot},var2,whichuv{nplot},iz(nplot),n,cont1{nplot},cont2,cont2h,nskip);
         title([num2str(z(iz(nplot))) units ', ' char(ut2ds(map.time(n),'mm/dd'))]);
      end
      ii=ctd.time >= map.time(n)-dtimeut & ctd.time <= map.time(n)+dtimeut;
      hold on; plot(ctd.lon(ii),ctd.lat(ii),'.','color',[0.3 0.3 0.3]); hold off
   end
   drawnow;
   frame = getframe(gcf);
   writeVideo(v,frame);
   clf;
end

close(v);
