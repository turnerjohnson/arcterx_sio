function makeplots(dep)

tbar=5:30;
sbar=34:0.05:35;
n2bar=0:1:12;
ubar=-1:0.1:1;
sigbar=20:0.5:30;
sigbarh=20:30;
plotspath='spray_plots';
[~,~]=mkdir(plotspath);

d=cell(length(dep),1);
for n=1:length(dep)
   d{n}=load(dep(n).name);
end

pvela2bin(d{1}.bindata,':',d{2}.bindata,':',d{3}.bindata,':',d{4}.bindata,':');
filename=fullfile(plotspath,'arcterx_vel.png');
print('-dpng',filename);

for n=1:length(dep)
   bindata=d{n}.bindata;
   sigint=interp_sigma(bindata,20:0.1:30);
   for m=1:length(dep(n).dive)
      if length(dep(n).dive{m}) == 1
         nn=dep(n).dive{m}:length(d{n}.bindata.time);
      else
         nn=dep(n).dive{m};
      end
      deppath=fullfile(plotspath,dep(n).name);
      [~,~]=mkdir(deppath);
      psection(bindata,'t','sigma',nn,tbar,sigbar,sigbarh,'dist');
      filename=fullfile(deppath,[dep(n).name '_' num2str(min(nn)) '-' num2str(max(nn)) '_t.png']);
      set(gca,'ylim',[0 400]);
      titl=get(gca,'Title');
      titl.String=['Spray ' num2str(d{n}.satdata.sn) ', ' titl.String];
      print('-dpng',filename);
      
      psection(bindata,'s','sigma',nn,sbar,sigbar,sigbarh,'dist');
      filename=fullfile(deppath,[dep(n).name '_' num2str(min(nn)) '-' num2str(max(nn)) '_s.png']);
      set(gca,'ylim',[0 400]);
      titl=get(gca,'Title');
      titl.String=['Spray ' num2str(d{n}.satdata.sn) ', ' titl.String];
      print('-dpng',filename);
      
      [~,dsigmadz]=gradient(bindata.sigma,bindata.depth(2)-bindata.depth(1));
      dsigmadz(dsigmadz < 0)=0;
      bindata.buoyfreq=sqrt(9.8/1027*dsigmadz)/(2*pi)*3600;
      psection(bindata,'buoyfreq','sigma',nn,n2bar,sigbar,sigbarh,'dist');
      filename=fullfile(deppath,[dep(n).name '_' num2str(min(nn)) '-' num2str(max(nn)) '_buoyfreq.png']);
      set(gca,'ylim',[0 400]);
      hc=get(gcf,'children');
      hc(1).YLabel.String='Buoyancy Frequency (cycles/h)';
      titl=get(gca,'Title');
      titl.String=['Spray ' num2str(d{n}.satdata.sn) ', ' titl.String];
      print('-dpng',filename);
      
      jj=find(~all(isnan(bindata.abs)));
      kk=intersect(jj,nn);
      
%       psection(bindata,'udop','sigma',kk,ubar,sigbar,sigbarh,'dist');
%       filename=fullfile(deppath,[dep(n).name '_' num2str(min(nn)) '-' num2str(max(nn)) '_udop.png']);
%       set(gca,'ylim',[0 400]);
%       print('-dpng',filename);
%       
%       psection(bindata,'vdop','sigma',kk,ubar,sigbar,sigbarh,'dist');
%       filename=fullfile(deppath,[dep(n).name '_' num2str(min(nn)) '-' num2str(max(nn)) '_vdop.png']);
%       set(gca,'ylim',[0 400]);
%       print('-dpng',filename);
      
      bindata.abs(bindata.abs > 90)=nan;
      psection(bindata,'abs','sigma',kk,[],sigbar,sigbarh,'dist');
      filename=fullfile(deppath,[dep(n).name '_' num2str(min(nn)) '-' num2str(max(nn)) '_abs.png']);
      set(gca,'ylim',[0 400]);
      titl=get(gca,'Title');
      titl.String=['Spray ' num2str(d{n}.satdata.sn) ', ' titl.String];
      print('-dpng',filename);
      
      psection(bindata,'fl','sigma',nn,[],sigbar,sigbarh,'dist');
      filename=fullfile(deppath,[dep(n).name '_' num2str(min(nn)) '-' num2str(max(nn)) '_fl.png']);
      set(gca,'ylim',[0 400]);
      titl=get(gca,'Title');
      titl.String=['Spray ' num2str(d{n}.satdata.sn) ', ' titl.String];
      print('-dpng',filename);
      
      pths(bindata,nn,'time');
      filename=fullfile(deppath,[dep(n).name '_' num2str(min(nn)) '-' num2str(max(nn)) '_ths.png']);
      titl=get(gca,'Title');
      titl.String=['Spray ' num2str(d{n}.satdata.sn) ', ' titl.String];
      print('-dpng',filename);
      
      psigint(sigint,'theta',nn,20:30,'dist');
      filename=fullfile(deppath,[dep(n).name '_' num2str(min(nn)) '-' num2str(max(nn)) '_sigint.png']);
      titl=get(gca,'Title');
      titl.String=['Spray ' num2str(d{n}.satdata.sn) ', ' titl.String];
      print('-dpng',filename);
   end
end
