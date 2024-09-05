function [E,G] = computeCovMtrxxy(x,y,time,xgxg,ygyg,timegtimeg,gauss,noise)

% Compute covariance matricies. Outputs:
    % E: data covariance matrix
    % G: model-data covariance matrix
%
% K.Zaba, May22,2014

% Parameters
nmaxfac = 0.9;  % maximum fraction of nonzero elements in covariance matrices
covmin = 1e-3;

% construct data covariance matrix
npt = length(x);
nmax = round(nmaxfac*npt^2);
rows = zeros(nmax,1);
cols = zeros(nmax,1);
vals = zeros(nmax,1);
rows(1:npt) = 1:npt;
cols(1:npt) = 1:npt;
vals(1:npt) = noise;
ilast = npt;
for ipt = 1:npt % loop through obs points to directly construct sparse matrix
    xx = x(ipt)-x(ipt:npt);
    yy = y(ipt)-y(ipt:npt);
    timetime = (time(ipt)-time(ipt:npt))/86400;
    dd = gauss(1)*exp(gauss(2)*xx.^2+gauss(3)*yy.^2+gauss(4)*timetime.^2); % covar of ipt with others
    ii = find(dd >= covmin);
    jj = ipt*ones(size(ii));
    dd = dd(ii);
    ii = ii+ipt-1;
    tmprows = [ii; jj(2:end)];
    tmpcols = [jj; ii(2:end)];
    tmpvals = [dd; dd(2:end)];
    nval = length(tmprows);
    rows(ilast+1:ilast+nval) = tmprows;
    cols(ilast+1:ilast+nval) = tmpcols;
    vals(ilast+1:ilast+nval) = tmpvals;
    ilast = ilast+nval;
    if ilast > nmax
        fprintf('Warning: preallocated matrix size exceeded at step %g\n',ipt)
    end
end
rows = rows(1:ilast);
cols = cols(1:ilast);
vals = vals(1:ilast);
E = sparse(rows,cols,vals,npt,npt);

% construct model-data covariance matrix
nmap=length(xgxg);
nmax = round(nmaxfac*npt*nmap);
rows = zeros(nmax,1);
cols = zeros(nmax,1);
vals = zeros(nmax,1);
ilast = 0;
for ipt = 1:nmap % loop through map points to directly construct sparse matrix
    xx = xgxg(ipt)-x;
    yy = ygyg(ipt)-y;
    timetime = (timegtimeg(ipt)-time)/86400;
    dd = gauss(1)*exp(gauss(2)*xx.^2+gauss(3)*yy.^2+gauss(4)*timetime.^2); % covar of ipt with others
    ii = find(dd>=covmin);
    jj = ipt*ones(size(ii));
    dd = dd(ii);
    tmprows = jj;
    tmpcols = ii;
    tmpvals = dd;
    nval = length(tmprows);
    rows(ilast+1:ilast+nval) = tmprows;
    cols(ilast+1:ilast+nval) = tmpcols;
    vals(ilast+1:ilast+nval) = tmpvals;
    ilast = ilast+nval;
    if ilast > nmax
        fprintf('Warning: preallocated matrix size exceeded at step %g\n',ipt)
    end
end
rows = rows(1:ilast);
cols = cols(1:ilast);
vals = vals(1:ilast);
G = sparse(rows,cols,vals,nmap,npt);
