function A=meanfit(ctd,year0,vars)
% function A=meanfit_c(ctd,vars)
% Calculates mean structure A given ctd structure 
% vars: cell list of variables
%
% D. Rudnick, 10 February 2020

% number of levels and number of bins in x
nlev=size(ctd.t,1);

% initialize
A.missions=ctd.missions;
A.year0=year0;
A.depth=ctd.depth;
A.u.constant=nan;
A.u.time=nan;
A.u.x=nan;
A.u.y=nan;
A.v.constant=nan;
A.v.time=nan;
A.v.x=nan;
A.v.y=nan;

% get time, x, y
time=ut2dn(ctd.time)-datenum(year0,1,1);
x=ctd.dist;
y=ctd.offset;

%do each variable
for ii2=1:length(vars)
   A.(vars{ii2}).constant=nan(nlev,1);
   A.(vars{ii2}).time=nan(nlev,1);
   A.(vars{ii2}).x=nan(nlev,1);
   A.(vars{ii2}).y=nan(nlev,1);
   t=ctd.(vars{ii2});
   
   %do each level
   for n=1:nlev
      kk=~isnan(t(n,:));
      tlev=t(n,kk)';
      timelev=time(kk);
      xlev=x(kk);
      ylev=y(kk);
      G=[ones(size(timelev)) timelev xlev ylev];
      mm=G\tlev;
      A.(vars{ii2}).constant(n)=mm(1);
      A.(vars{ii2}).time(n)=mm(2);
      A.(vars{ii2}).x(n)=mm(3);
      A.(vars{ii2}).y(n)=mm(4);
   end
end

% now do u and v
kk=~isnan(ctd.u);
tlev=ctd.u(kk);
timelev=time(kk);
xlev=x(kk);
ylev=y(kk);
G=[ones(size(timelev)) timelev xlev ylev];
mm=G\tlev;
A.u.constant=mm(1);
A.u.time=mm(2);
A.u.x=mm(3);
A.u.y=mm(4);

tlev=ctd.v(kk);
mm=G\tlev;
A.v.constant=mm(1);
A.v.time=mm(2);
A.v.x=mm(3);
A.v.y=mm(4);

% and ualong and uacross if necessary
if isfield(ctd,'ualong')
   tlev=ctd.ualong(kk);
   mm=G\tlev;
   A.ualong.constant=mm(1);
   A.ualong.time=mm(2);
   A.ualong.x=mm(3);
   A.ualong.y=mm(4);
   
   tlev=ctd.uacross(kk);
   mm=G\tlev;
   A.uacross.constant=mm(1);
   A.uacross.time=mm(2);
   A.uacross.x=mm(3);
   A.uacross.y=mm(4);
end
