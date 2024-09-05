function [map,mse] = generateMapxy(A,B,d,gauss,varargin)

% Generate objective map.
    % A: data-data covariance
    % B: model-data covariance
    % d: detrended data, nans removed
    % gauss: 4 parameters of Gaussian covariance matrix
    % varargin: domse=true calculates the mean square errorl domse=false
        % does not calculate the mean square error and is much faster and
        % less memory intensive; if domse is not provided, the mean square
        % error is calculated by default
% 
% Dan Rudnick, 14 Feb 2020

% Parse Inputs
if isempty(varargin)
    domse = true;
else
    domse = varargin{1};
end

% make map and mse, if desired
C = A\d;
map = B*C;
if domse
%    mse = diag(1-B*(A\B')/gauss(1)); % this step is very slow
   mse=1-sum(B'.*(A\B'))/gauss(1);
else
   mse=[];
end