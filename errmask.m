function mapmask=errmask(map,errthresh)

% Mask map using error mask
%
% K.Zaba; Aug8,2016

mapmask=map;
mapmask.errthresh=errthresh;

vars = {'t','s','u','v','udop','vdop','fl','abs','ualong','uacross','udopalong','udopacross'};
for k=1:length(vars)
    var0 = vars{k};
    if isfield(mapmask,var0)
       mask = mapmask.err.(var0)>errthresh;
       mapmask.(var0)(mask)=nan;
%        if strcmp(var0,'u') && isfield(mapmask,'uacross')
%           mapmask.uacross(mask) = nan;
%        elseif strcmp(var0,'v') && isfield(mapmask,'ualong')
%           mapmask.ualong(mask) = nan;
%        end
%        if strcmp(var0,'udop') && isfield(mapmask,'udopacross')
%           mapmask.udopacross(mask) = nan;
%        elseif strcmp(var0,'vdop') && isfield(mapmask,'udopalong')
%           mapmask.udopalong(mask) = nan;
%        end
    end
end

errts = max(mapmask.err.t,mapmask.err.s);
varsder = {'theta','rho','sigma','sigmastable'};
for k=1:length(varsder)
    var0 = varsder{k};
    if isfield(mapmask,var0)
        mask = errts>errthresh;
        mapmask.(var0)(mask)=nan;
    end
end

if isfield(mapmask,'udop')
   errtsu=max(errts,mapmask.err.udop);
elseif isfield(mapmask,'udopalong')
   errtsu=max(errts,mapmask.err.udopalong);
else
   errtsu=errts;
end
varstsu={'geou','geov','pv','pvvert','pvhoriz','w','ro','div'};
for k=1:length(varstsu)
    var0 = varstsu{k};
    if isfield(mapmask,var0)
        mask = errtsu>errthresh;
        mapmask.(var0)(mask)=nan;
    end
end

