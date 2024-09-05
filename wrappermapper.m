function [A,ctd,config,map,mapem,mapemsig,mapw,mapwem,mapwemsig]=wrappermapper(year)
% function [A,ctd,config,map,mapem,mapw,mapwem]=wrappermapper(year)

% function to do the mapping thing

% anisotropic 

% cell array of input variables
vars = {'ualong','uacross','t','s','udopalong','udopacross','fl','abs'};

vars2d = {'t','s','udop','vdop','udopalong','udopacross','fl','abs'};
    % vars={'u','v','t','s'};
    % vars2d={'t','s'};
errthresh = 0.3;
    % nreflev = 65;

% allmeanxy(): combine missions, rotate velocity coordinate system
    % in: 'year': string (dev purposes), 'vars': variable list
    % outputs:  'ctd': (structure: joined missions' processed data)
            %   'A': annual mean field of CTD
[A,ctd] = allmeanxy(year,vars2d);

% mapxyfaster(): 
% config: info about constants
% map 
[config,map] = mapxyfaster(ctd,A,vars);

save(['map' num2str(year)],'-v7.3'); %? why save map here? what does v7.3 mean?

map = addderivedvarsxy(map,ctd);
map = computePVxy(map);
map = computeStableDensity(map);
map = computeGeoVxy(map);
    % map=computeGeoVxy_reflevel(map,nreflev);

%calculate w
mapw = computeOmegaw2_stretch(map);
mapw = computeageouv(mapw);
mapw = computePVxy_geo(mapw);

% error mask: where no data, error is large
mapem = errmask(map,errthresh);
mapwem = errmaskw(mapw,map,errthresh);

% add sigma 
mapemsig = interpmap_sigma2(mapem,22:0.1:27);
mapwemsig = interpmap_sigma2(mapwem,22:0.1:27);
save(['map' num2str(year)],'-v7.3'); %
