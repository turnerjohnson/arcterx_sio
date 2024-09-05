function [config,map] = mapxyfaster(ctd,A,vars)
% Wrapper function for computing a map of Calypso data.
% The map is 3d in x, y, and t at each depth
% The resulting grids are 4D: <ndepth,nx,ny,ntime>. Input parameters are:
    % ctd: ctd data in ctd structure
    % annual cycle structure (A)
    % vars:   input variables like 'u','v','t','s','udop','vdop'
    %         for multiple variables, 'vari' needs to be a cell array, 
    %         e.g. {'u','v','t','s'}
% This script consists of two loops: the outer loop iterates through all
% input variables and the inner loop iterates through all depth levels. The
% output structures are:
    % 'config':     record of parameters used in the mapping computation
    % 'anommap':    objective maps of anomalies from the annual cycle
% In addition, the 'anommap' structure includes error maps for each 
% variable and the depth/spatial/temporal grids. 
% 
%
% K.Zaba May13,2016 
% D. Rudnick 10 May 2019 - Calypso version
% D. Rudnick, 12 Feb 2020 - xy version
% D. Rudnick, 23 Jun 2022 - a try at making things faster by using
% previously calculated error maps

% Parameters
Lx = 100;                        % length scale, in km in x_across dir ('x')
Ly = 100;                        % length scale, in km in x_along dir ('y')
Lt = 20;                        % time scale, in days
gauss = [1 -1/Lx^2 -1/Ly^2 -1/Lt^2];    % [amp, -1/Lx^2, -1Ly^2 -1/Lt^2]
% grid spacing: 10km
xg=(-250:25:250)';                     % x grid of map (km)
yg=(-250:25:250)';                     % y grid of map (km)
dt = 1;                        % temporal resolution of map (days)
noise = 0.1;                    % noise to signal ratio

nx=length(xg);
ny=length(yg);
nz=length(ctd.depth);

% Variables
if ~iscell(vars)
    vars = {vars};
end
nvari = length(vars);
variStr1D = {'u','v','ualong','uacross'};          % 1D variables
variStr2D = {'t','s','udop','vdop','udopalong','udopacross','fl','abs'};  % 2D variables
ii1DVar=ismember(vars,variStr1D);  n1D=length(find(ii1DVar));   % 1D variables
ii2DVar=ismember(vars,variStr2D);  n2D=length(find(ii2DVar));   % 2D variables
vars = [vars(ii1DVar) vars(ii2DVar)];                           
varitype = [repmat({'1D'},[1,n1D]) repmat({'2D'},[1,n2D])];

% Generate temporal grid
dvstart = datevec(ut2dn(min(ctd.time)));
dvend = datevec(ut2dn(max(ctd.time)));
tstart = dn2ut(datenum(dvstart(1),dvstart(2),dvstart(3))); % first date to do
tend = dn2ut(datenum(dvend(1),dvend(2),dvend(3)));     % last date to do
timeg=(tstart:(dt*86400):tend)';
% timeg=(dn2ut(datenum(2019,4,24)):(dt*86400):dn2ut(datenum(2019,4,25)))';
ntime=length(timeg);

% Configuration Parameters
timestamp = java.lang.System.currentTimeMillis/1000;    % timestamp for map computation, in UTC
config = struct('Lx',Lx,'Ly',Ly,'Lt',Lt,'noise',noise,'timestamp',timestamp);

% Initialize output structure map
map.missions = ctd.missions;                           % mission ID
map.x = xg;                          % x grid
map.y = yg;                          % x grid
map.time = timeg;                          % temporal grid
map.depth = ctd.depth;                       % depth levels
for ivar = 1:nvari
    if ismember(vars{ivar},variStr1D)
        grdvari = NaN*zeros(1,nx,ny,ntime);              % 1D data
    elseif ismember(vars{ivar},variStr2D)
        grdvari = NaN*zeros(nz,nx,ny,ntime);       % 2D data
    end
    map.(vars{ivar})     = grdvari;          % anomaly
    map.err.(vars{ivar}) = grdvari;          % error
end
map = orderfields(map,['missions','depth','time','x','y',vars,'err']);

% 3D grid for map
[xgxg,ygyg,timegtimeg]=meshgrid(xg,yg,timeg);

% Loop through variables
for ivar=1:nvari
    
    vari0 = vars{ivar};
    disp(['Computing maps for vari: ',vari0])
    
    % Covariance matrices (A,B) only computed once for each
    % variable type (1D and 2D).  
    if ivar == 1 %first 1D variable
       x = ctd.distu;
       y = ctd.offsetu;
       time = ctd.timeu;
       ndepth = 1;
       [E,G] = computeCovMtrxxy(x,y,time,xgxg(:),ygyg(:),timegtimeg(:),gauss,noise);
    elseif ivar == n1D+1  % first 2D variable
       x = ctd.dist;
       y = ctd.offset;
       time = ctd.time;
       ndepth = length(ctd.depth);
       [E,G] = computeCovMtrxxy(x,y,time,xgxg(:),ygyg(:),timegtimeg(:),gauss,noise);
    end
    
    % Data
    if strcmp(varitype{ivar},'1D')
        data = ctd.(vari0)';
    elseif strcmp(varitype{ivar},'2D')
        data = ctd.(vari0);
    end
       
    % Loop through depth levels
    for idepth = 1:ndepth
        
        % Data
        data0 = data(idepth,:)';
        
        % Decide whether the error map needs to be computed.  Search
        % through all previous variables of the same type and depths to see
        % if any had the same data locations/times
        if ivar ==1             % first 1D variable
           compErr = true;
        elseif ivar <= n1D      % second or higher 1D variable
           compErr = false;
           varErr = vars{1};
           idepErr = 1;
        elseif ivar == n1D+1    % first 2D variable
           compErr = true;
           if idepth > 1        % loop though previous depths if not first depth
              for itrydep = 1:idepth-1  % next try previous depths of this variable
                 datap = ctd.(vari0)(itrydep,:)';
                 if all(~xor(~isnan(datap),~isnan(data0)))
                    compErr = false;
                    varErr = vari0;
                    idepErr = itrydep;
                    break
                 end
              end
           end
        elseif ivar > n1D+1     % second or higher 2D variable
           compErr = true;
           for itryvar = n1D+1:ivar-1 % first try previous 2D variables
              for itrydep = 1:ndepth
                 datap = ctd.(vars{itryvar})(itrydep,:)';
                 if all(~xor(~isnan(datap),~isnan(data0)))
                    compErr = false;
                    varErr = vars{itryvar};
                    idepErr = itrydep;
                    break
                 end
              end
              if compErr == false
                 break
              end
           end
           for itrydep = 1:idepth-1  % next try previous depths of this variable
              datap = ctd.(vari0)(itrydep,:)';
              if all(~xor(~isnan(datap),~isnan(data0)))
                 compErr = false;
                 varErr = vari0;
                 idepErr = itrydep;
                 break
              end
           end
        end
                
        % Remove nans
        igood = (~isnan(x) & ~isnan(y) & ~isnan(time) & ~isnan(data0));
        x0 = x(igood);
        y0 = y(igood);
        time0 = time(igood);
        data0 = data0(igood); 
        
        if length(data0) > 1
           % Remove mean
           meandata0 = meanmodel(A,vari0,idepth,time0,x0,y0);
           data0 = data0-meandata0;
           
           % Generate 'map'
           [mappy,mse] = generateMapxy(E(igood,igood), G(:,igood),data0,gauss,compErr);
           mappy=reshape(mappy+meanmodel(A,vari0,idepth,timegtimeg(:),xgxg(:),ygyg(:)),size(xgxg));
           map.(vari0)(idepth,:,:,:) = permute(mappy,[2 1 3]);
           if compErr
              mse=reshape(mse,size(xgxg));
              map.err.(vari0)(idepth,:,:,:) = permute(mse,[2 1 3]);
           else
              map.err.(vari0)(idepth,:,:,:) = map.err.(varErr)(idepErr,:,:,:);
           end
           
           if ndepth>1
              disp(['    ',num2str(map.depth(idepth)),' m',', compErr=',num2str(compErr),', varErr=',varErr,', idepErr=',num2str(idepErr)]);
           end
        end
    end
end
