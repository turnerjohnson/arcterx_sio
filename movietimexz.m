function movietimexz(map,var1str,iy,cont1)

ntend=length(map.time);

if isfield(map,'sigmastable')
   cont2=20:0.25:30;
   cont2h=20:1:30;
   var2='sigmastable';
else
   cont2=0:10:400;
   cont2h=0:20:400;
   var2='depth';
   z = map.sigma;
   units = ' kg/m^3';
end

v = VideoWriter('ARCTERX','MPEG-4');
v.FrameRate=5;
open(v);

for n=1:ntend
   if any(any(nanmean(map.(var1str)(:,:,iy,n),2)))
      pxz_xy(map,var1str,var2,iy,n,cont1,cont2,cont2h);
   end
   drawnow;
   frame = getframe(gcf);
   writeVideo(v,frame);
   clf;
end

close(v);
