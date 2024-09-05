function bindata=binsolo(data,pmin,pstep,pmax,pd,exclude)
% function bindata=binspray(data,pmin,pstep,pmax,pd,exclude) makes the structure
% bindata of data binnned in pressure on the grid [pmin:pstep:pmax].  The
% string pd can be either 'p' or 'd' to bin in pressure or depth. The
% optional string exclude can be 'none','bad','questionable' to tell the
% routine which points to exclude. The default is 'bad'.
% 

% Dan Rudnick, September 1, 2008
% Robert Todd, 15 September 2008, added flag support
%              27 October 2008, added time of creation field (bintime)
%              7 July 2009, added timeu,latu,lonu,u,v fields
%              8 July 2009, bindata.u,v set to NaN where quality flag not Gps_Good
% Dan Rudnick, March 2, 2011, added oxygen, and logic to tell whether fl
% and ox are present
% Dan Rudnick, November 27, 2013, adapted from binspray


% Define flags
Ctd_Sensor_Off = 9;
Ctd_Bad=7;
Ctd_Questionable=3;
Gps_Good = 0;

bindata.time=data.time(:,2);
bindata.lat=data.lat(:,2);
bindata.lon=data.lon(:,2);

if pd == 'p'
   bindata.p=(pmin:pstep:pmax)';
   pstr='p';
elseif pd =='d'
   bindata.depth=(pmin:pstep:pmax)';
   pstr='depth';
else
   error('pd must be ''p'' (pressure) or ''d'' (depth)');
end

if nargin == 6
    switch exclude(1)
        case 'n'
            maxflag = Ctd_Sensor_Off;
        case 'b'
            maxflag = Ctd_Bad;
        case 'q'
            maxflag = Ctd_Questionable;
        otherwise
            error('exclude must be ''none'', ''bad'' or ''questionable''');
    end
else
    maxflag = Ctd_Bad;
end

np=length(bindata.(pstr));
nt=length(bindata.time);
bindata.t=nan(np,nt);
bindata.s=nan(np,nt);
bindata.theta=nan(np,nt);
bindata.sigma=nan(np,nt);
bindata.rho=nan(np,nt);

for n=1:nt
   try
      if ~isempty(data.(pstr){n}) % then there are data to bin
         ibin=round((data.(pstr){n}-pmin)/pstep)+1;
         for m=1:np
            iit = (ibin == m) & data.qual.t{n} < maxflag;
            iis = (ibin == m) & data.qual.s{n} < maxflag;
            ii = iit & iis;
            bindata.t(m,n)=nanmean(data.t{n}(iit));
            bindata.s(m,n)=nanmean(data.s{n}(iis));
            bindata.theta(m,n)=nanmean(data.theta{n}(ii));
            bindata.sigma(m,n)=nanmean(data.sigma{n}(ii));
            bindata.rho(m,n)=nanmean(data.rho{n}(ii));
         end
      end
   catch err
      display(err.identifier);
   end
end

bindata.bintime = round(dn2ut(now));
