function ctd = combineMissionsxy(year,vars,depthmin)

% Concatenate data from all missions along a particular Calypso line. 
% Use variables in cell array vars. 

% ctd=combineMissionsxy(year,vars,depthmin)

% input parameters:
    % 'year':       2018 or 2019 or 2022
    % 'vars':       cell array of variables
    % 'depthmin':   minimum depth of profile to include

% K.Zaba, May22,2014
% D. Rudnick, September 8, 2016 - update to use cugn.txt
% D. Rudnick, March 12, 2018 - update to use vars and data in binmat
% D. Rudnick, April 22, 2018 - roger version
% D. Rudnick, 10 May 2019 - Calypso version
% D. Rudnick, 31 May 2019 - Calypso 2019 version
% D. Rudnick, 10 February 2020 - xy version
% D. Rudnick, 1 March 2021 - add mission identifier
% D. Rudnick, 25 March 2021 - add minimum depth
% D. Rudnick, 1 December 2023 - ARCTERX version

% CTD Variables
ctd1D = {'time','lon','lat','dist','offset', ...
         'timeu','lonu','latu','distu','offsetu','u','v'};
ctd2D = vars;
nctd1D = length(ctd1D); 
nctd2D = length(ctd2D); 
depthLevels = (10:10:1000)'; 
nz = length(depthLevels);
nzmin = round(depthmin/10);

% line stuff 
ctd.y0=19.5;
ctd.y1=19.5;
ctd.x0=141.75;
ctd.x1=142;
ctd.time0=dn2ut(datenum(2023,5,30));
ctd.time1=dn2ut(datenum(2023,8,1));
ctd.depthmin=depthmin;
%ctd.missions={'23503301','23503701','23504501','23505601','osu684','osu685','osu686','sg141','sg526','sg528','sg687','8797','8799','8800','8801','8802'};
      % dan's spray (naming scheme: yymm serial depl of month)  %osu gliders, seagliders
ctd.missions={'23503301','23503701','23504501','23505601','osu684','osu685','osu686','sg141','sg526','sg528','sg687'};


% Initialize CTD Output 
ctd.depth = depthLevels;
ctd.missid = [];
for iictd1D = 1:nctd1D
    ctd.(ctd1D{iictd1D}) = [];
end
for iictd2D = 1:nctd2D
    ctd.(ctd2D{iictd2D}) = [];
end

for imiss=1:length(ctd.missions)
   filename=[ctd.missions{imiss} '_bin.mat'];
   if exist(filename,'file')
      load(filename,'bindata');
   else
      load(ctd.missions{imiss},'bindata');
   end
   
   [bindata.dist,bindata.offset] = ll2do(bindata.lon,bindata.lat,ctd.x0,ctd.y0,ctd.x1,ctd.y1);
   
   if isfield(bindata,'lonu') %Spray
      [bindata.distu,bindata.offsetu] = ll2do(bindata.lonu,bindata.latu,ctd.x0,ctd.y0,ctd.x1,ctd.y1);
   elseif length(bindata.depth) == 25 %kluge for SOLO-IIs
      bindata.distu=bindata.dist;
      bindata.offsetu=bindata.offset;
      bindata.lonu=bindata.lon;
      bindata.latu=bindata.lat;
      bindata.timeu=bindata.time;
      bindata.u=nan(size(bindata.time));
      bindata.v=nan(size(bindata.time));
   else %Slocum
      bindata.distu=bindata.dist;
      bindata.offsetu=bindata.offset;
      bindata.lonu=bindata.lon;
      bindata.latu=bindata.lat;
      bindata.timeu=bindata.time;
   end
      
   %find indices of complete profiles at least to the minimum depth and greater than
   %ctd.time0, also get rid of profiles that are all nan
%    [~,nn]=sort(~isnan(bindata.t));
%    kk=nn(end,:)';
%    jj=kk >= nzmin & bindata.time > ctd.time0 & ~all(isnan(bindata.t))';
   jj=all(~isnan(bindata.t(1:nzmin,:)))' & bindata.time > ctd.time0 & bindata.time < ctd.time1;
   
   ctd.missid = [ctd.missid; imiss*ones(size(bindata.time(jj)))];
   
   for iictd1D=1:nctd1D
      ctd.(ctd1D{iictd1D}) = [ctd.(ctd1D{iictd1D}); bindata.(ctd1D{iictd1D})(jj)];
   end
   
   nt=length(bindata.time(jj));
   for iictd2D=1:nctd2D
      tmp=nan(nz,nt);
      if isfield(bindata,ctd2D{iictd2D})
         nzz=min(nz,size(bindata.(ctd2D{iictd2D}),1));
         tmp(1:nzz,:)=bindata.(ctd2D{iictd2D})(1:nzz,jj);
      end
      ctd.(ctd2D{iictd2D})=[ctd.(ctd2D{iictd2D}) tmp];
   end
      
   clear('bindata');
end
