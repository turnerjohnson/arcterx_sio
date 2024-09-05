function [c1,h1,c2,h2] = pxz_xy(mapstruct,var1str,var2str,iy,itime,cont1,cont2,cont2h)

% Plot section with x=dist, y=depth, color=var1, contour=var2.
% input parameters:
    % 'mapstruct':  'map','avmap', or 'anommap'
    % 'var1str':    variable to be plotted with a color contour
    % 'var2str':    variable to be plotted with a black contour
    % 'itime':      time indices to average over to make plot
    % 'cont1':      contour vector for var1 (if null, then matlab decides
    %               contour levels)
    % 'cont2':      contour vector for var2 (if null, then matlab decides
    %               contour levels)
    % 'cont2h':     heavy/labeled contour vector for var2 (if null, no
    %               heavy contours)
% 
%
% K.Zaba: Jun15,2016
% D. Rudnick, July 29, 2016
% D. Rudnick, 23 April 2020

% Data
x = mapstruct.x;
if isvector(mapstruct.depth)
    z = mapstruct.depth;
    ylim = [0 mapstruct.depth(end)];
    ylab='Depth (m)';
else
    z = mapstruct.sigma;
    ylim = [mapstruct.sigma(1) mapstruct.sigma(end)];
    ylab='Potential Density (kg/m^3)';
end

if ndims(mapstruct.(var1str)) == 4
   var1 = squeeze(mapstruct.(var1str)(:,:,iy,itime));
   if ut2dn(mapstruct.time(end))-ut2dn(mapstruct.time(1)) < 366
      dnfmt='mm/dd';
   else
      dnfmt='yyyy/mm/dd';
   end
   if isscalar(itime)
      dnstr = char(ut2ds(mapstruct.time(itime),dnfmt));
   elseif ndims(mapstruct.(var1str)) == 3
      dnstr=[char(ut2ds(mapstruct.time(itime(1)),dnfmt)) ' - ' char(ut2ds(mapstruct.time(itime(end)),dnfmt))];
   end
else
   var1=mapstruct.(var1str);
   dnstr='Mean';
end
titl=[num2str(mapstruct.y(iy)) ' km, ' dnstr];

% Plot 
[c1,h1] = contourf(x,z,var1,cont1,'linestyle','none');
if isempty(cont1)
   set(h1,'levellistmode','auto');
   cont1=get(h1,'levellist');
   dc=cont1(end)-cont1(end-1);
   if any(strcmp(var1str,{'geov','geou','u','v','udop','vdop','udopalong','udopacross','ro','w','ua','va','rogeo'})) || (isfield(mapstruct,'isanom') && mapstruct.isanom)
      if abs(cont1(end)) > abs(cont1(1))
         cont1=-cont1(end):dc/2:cont1(end)+eps;
      else
         cont1=cont1(2)-dc:dc/2:-cont1(2)+dc+eps;
      end
   else
      cont1=cont1(2)-dc:dc/2:cont1(end)+dc+eps;
   end
   cla('reset');
   [c1,h1] = contourf(x,z,var1,cont1,'linestyle','none');
end

hold on;

if ~isempty(var2str)
   if ndims(mapstruct.(var2str)) == 4
      var2 = squeeze(mapstruct.(var2str)(:,:,iy,itime));
   else
      var2=mapstruct.(var2str);
   end
   cont2(ismember(cont2,cont2h)) = [];     % contour levels
   [c2,h2] = contour(x,z,var2,cont2,'k');
   [ch,hh] = contour(x,z,var2,cont2h,'k','linewidth',2);
   clabel(ch,hh);
   if isempty(cont2)
      set(h2,'levellistmode','auto');
   end
end

set(gca,'ydir','rev','xlim',[min(x) max(x)],'ylim',ylim,'clim',[cont1(1) cont1(end)]);
xlabel('Distance (km)')
ylabel(ylab)
title(titl)

% % Plot Bathymetry
% if isvector(mapstruct.depth)
%    topo=load(['topo' mapstruct.line(1:2)]);
%    area(topo.dist,topo.topo,'basevalue',ylim(2),'facecolor','k');
% end

extlab = var1str;
if strcmp(var1str,'t')
   extlab = 'Temperature (\circC)';
   cmap=jet(64);
elseif strcmp(var1str,'s')
   extlab = 'Salinity (psu)';
   cmap=jet(64);
elseif strcmp(var1str,'geov')
   extlab='Along-front Geostrophic Velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'geou')
   extlab='Across-front Geostrophic Velocity (m/s)';
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
elseif strcmp(var1str,'ua')
   extlab='Across-front ageostrophic velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'va')
   extlab='Along-front ageostrophic velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'udopacross')
   extlab='Across-front Velocity (m/s)';
   cmap=redblue(length(cont1)-1);
elseif strcmp(var1str,'w')
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