function mapmask=errmaskw(mapw,map,errthresh)

% Mask map using error mask special for w
%

mapmask=mapw;
mapmask.errthresh=errthresh;

errt=cat(1,map.err.t(1,:,:,:),map.err.t);
errs=cat(1,map.err.s(1,:,:,:),map.err.s);
if isfield(mapmask,'udop')
   erru=cat(1,map.err.udop(1,:,:,:),map.err.udop);
elseif isfield(mapmask,'udopalong')
   erru=cat(1,map.err.udopalong(1,:,:,:),map.err.udopalong);
else %no ADCP data, so just call the error zero so ts error wins
   erru=zeros(size(errt));
end


errts = max(errt,errs);
varsder = {'sigmastable'};
for k=1:length(varsder)
    var0 = varsder{k};
    if isfield(mapmask,var0)
        mask = errts>errthresh;
        mapmask.(var0)(mask)=nan;
    end
end

errtsu=max(errts,erru);
varstsu={'geou','geov','w','ua','va','pvgeo','pvvertgeo','pvhorizgeo','rogeo','divgeo'};
for k=1:length(varstsu)
    var0 = varstsu{k};
    if isfield(mapmask,var0)
        mask = errtsu>errthresh;
        mapmask.(var0)(mask)=nan;
    end
end

varsu = {'udopacross','udopalong'};
for k=1:length(varsu)
    var0 = varsu{k};
    if isfield(mapmask,var0)
        mask = erru>errthresh;
        mapmask.(var0)(mask)=nan;
    end
end
