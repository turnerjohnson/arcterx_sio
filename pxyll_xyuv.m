function [c1,h1,c2,h2] = pxyll_xyuv(mapstruct,var1str,var2str,whichuv,iz,itime,cont1,cont2,cont2h,nskip)

% Plot section with x=lon, y=lat, color=var1, contour=var2, quiver/arrows=velocity.
% input parameters:
    % 'mapstruct':  'map','avmap', or 'anommap'
    % 'var1str':    variable to be plotted with a color contour
    % 'var2str':    variable to be plotted with a black contour
    % 'whichuv':    velocity type: "geo", "obs", or "ageo"
    % 'iz':         z indices to average over to make plot
    % 'itime':      time indices to average over to make plot
    % 'cont1':      contour vector for var1 (if null, then matlab decides
    %               contour levels)
    % 'cont2':      contour vector for var2 (if null, then matlab decides
    %               contour levels)
    % 'cont2h':     heavy/labeled contour vector for var2 (if null, no
    %               heavy contours)
    % 'nskip':      interval to skip velocity quiver indices with (1:nskip:right bound)
% 
% K.Zaba: Jun15,2016
% D. Rudnick, July 29, 2016
% D. Rudnick, 23 April 2020

% Data
x = mapstruct.lon;
y = mapstruct.lat;

% [yy,xx]=meshgrid(y,x);
xx=x;
yy=y;

uscl=0.5;
coslat=cosd(mean(y(:)));

if isvector(mapstruct.depth)
    z = mapstruct.depth;
    units = ' m';
else
    z = mapstruct.sigma;
    units = ' kg/m^3';
end

if ndims(mapstruct.(var1str)) == 4
   var1 = squeeze(mapstruct.(var1str)(iz,:,:,itime));
   if ut2dn(mapstruct.time(end))-ut2dn(mapstruct.time(1)) < 366
      dnfmt='mm/dd';
   else
      dnfmt='yyyy/mm/dd';
   end
   if isscalar(itime)
      dnstr = char(ut2ds(mapstruct.time(itime),dnfmt));
   elseif ndims(mapstruct.(var1str)) == 4
      dnstr=[char(ut2ds(mapstruct.time(itime(1)),dnfmt)) ' - ' char(ut2ds(mapstruct.time(itime(end)),dnfmt))];
   end
end

% Plot 
[c1,h1] = contourf(x,y,var1,cont1,'linestyle','none');
if isempty(cont1)
   set(h1,'levellistmode','auto');
   cont1=get(h1,'levellist');
   dc=cont1(end)-cont1(end-1);
   if any(strcmp(var1str,{'geov','geou','u','v','udop','vdop','udopalong','udopacross','ro','w','ua','va','rogeo'})) || (isfield(mapstruct,'isanom') && mapstruct.isanom) || strncmp(var1str,'w',1)
      if abs(cont1(end)) > abs(cont1(1))
         cont1=-cont1(end):dc/2:cont1(end)+eps;
      else
         cont1=cont1(2)-dc:dc/2:-cont1(2)+dc+eps;
      end
   else
      cont1=cont1(2)-dc:dc/2:cont1(end)+dc+eps;
   end
   cla('reset');
   [c1,h1] = contourf(x,y,var1,cont1,'linestyle','none');
end

hold on;

if ~isempty(var2str)
   if ndims(mapstruct.(var2str)) == 4
      var2 = squeeze(mapstruct.(var2str)(iz,:,:,itime));
   else
      var2=mapstruct.(var2str);
   end
   cont2(ismember(cont2,cont2h)) = [];     % contour levels
   [c2,h2] = contour(x,y,var2,cont2,'k');
   [ch,hh] = contour(x,y,var2,cont2h,'k','linewidth',2);
   clabel(ch,hh);
   if isempty(cont2)
      set(h2,'levellistmode','auto');
   end
end

set(gca,'xlim',[min(x(:)) max(x(:))],'ylim',[min(y(:)) max(y(:))],'clim',[cont1(1) cont1(end)],'dataaspectratio',[1 coslat 1]);

if strcmp(whichuv,'geo')
   u=squeeze(mapstruct.geoull(iz,:,:,itime));
   v=squeeze(mapstruct.geovll(iz,:,:,itime));
   utitl='Geostrophic velocity';
   uscale=uscl;
   ulab='1 m/s';
elseif strcmp(whichuv,'obs')
   u=squeeze(mapstruct.udop(iz,:,:,itime));
   v=squeeze(mapstruct.vdop(iz,:,:,itime));
   utitl='Observed velocity';
   uscale=uscl;
   ulab='1 m/s';
elseif strcmp(whichuv,'ageo')
   u=squeeze(mapstruct.ua(iz,:,:,itime));
   v=squeeze(mapstruct.va(iz,:,:,itime));
   utitl='Ageostrophic velocity';
   uscl=150;
   uscale=uscl/10;
   ulab='0.1 m/s';
else
   error('geo or obs or ageo')
end

titl=[num2str(z(iz)) units ', ' dnstr];

quiver(xx(1:nskip:end,1:nskip:end),yy(1:nskip:end,1:nskip:end),uscl*u(1:nskip:end,1:nskip:end),uscl*coslat*v(1:nskip:end,1:nskip:end),0,'color','k');
quiver(min(x(:))+0.1,max(y(:))-0.1,uscale,0,0,'color','k');
text(min(x(:))+0.1,max(y(:))-0.1,ulab,'VerticalAlignment','Top');
text(min(x(:))+0.15+uscale,max(y(:))-0.1,utitl,'HorizontalAlignment','left','VerticalAlignment','middle');

xlabel('Longitude');
ylabel('Latitude');
title(titl);

extlab = var1str;
if strcmp(var1str,'t')
   extlab = 'Temperature (\circC)';
   cmap=jet(64);
elseif strcmp(var1str,'s')
   extlab = 'Salinity (psu)';
   cmap=jet(64);
elseif strcmp(var1str,'geou')
   extlab='Across-front Geostrophic Velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'geov')
   extlab='Along-front Geostrophic Velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'udop')
   extlab='Eastward Velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'vdop')
   extlab='Northward Velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'udopalong')
   extlab='Along-front Velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'udopacross')
   extlab='Across-front Velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strncmp(var1str,'w',1)
   extlab='Vertical Velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'theta')
   extlab='Potential Temperature (\circC)';
   cmap=jet(64);
elseif strcmp(var1str,'sigma')
   extlab='Potential Density (kg/m^3)';
   cmap=flipud(jet(64));
elseif strcmp(var1str,'rho')
   extlab='Density (kg/m^3)';
   cmap=flipud(jet(64));
elseif strcmp(var1str,'depth')
   extlab='Depth (m)';
   cmap=jet(64);
elseif strcmp(var1str,'fl')
   extlab='Chlorophyll fluorescence';
   cmap=jet(64);
elseif strcmp(var1str,'ro')
   extlab='\zeta/f';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'pv')
   extlab='Potential vorticity (m^{-1}s^{-1})';
   cmap=jet(64);
elseif strcmp(var1str,'rogeo')
   extlab='Geostrophic \zeta/f';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'pvgeo')
   extlab=' Geostrophic potential vorticity (m^{-1}s^{-1})';
   cmap=jet(64);
elseif strcmp(var1str,'abs')
   extlab='Acoustic backscatter (dB)';
   cmap=jet(64);
else
   cmap=jet(64);
end
if (isfield(mapstruct,'isanom') && mapstruct.isanom)
      cmap=redblue(length(cont1)-1);
      extlab=['Anomaly of ' extlab];
end
colormap(gca,cmap);
hb=extbar('v',cont1,extlab);
colormap(hb,cmap);

hold off;

end