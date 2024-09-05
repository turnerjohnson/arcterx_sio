function plotfaketracks(ctd,ndivestail,side)

dtime=ctd.time(2)-ctd.time(1);
timemax=max(ctd.time);
ndives=timemax./dtime;

v = VideoWriter('FakeTracks','MPEG-4');
v.FrameRate=5;
open(v);

for n=1:ndives
   for nn=1:max(ctd.missid)
      jj=ctd.missid == nn & ctd.time <= n*dtime & ctd.time >= (n-ndivestail)*dtime;
      hold on;
      plot(ctd.x(jj),ctd.y(jj),'-x');
   end
   set(gca,'dataaspectratio',[1 1 1],'xlim',[0 side],'ylim',[0 side]);
   hold off
   drawnow;
   frame = getframe(gcf);
   writeVideo(v,frame);
   clf;
end