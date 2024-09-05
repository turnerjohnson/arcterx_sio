function [A,ctd]=allmeanxy(year,vars)
% function [A,ctd]=allmeanxy(year,vars)
% Calculate mean for Calpyso given year.
% This is a wrapper function with lots of special numbers that could
% change.
%
% D. Rudnick, July 28, 2016
% D. Rudnick, April 20, 2018 - Roger version
% D. Rudnick 10 May 2019 - Calypso version
% D. Rudnick, 10 February 2020 - xy version
% D. Rudnick, 25 March 2021 - added depthmin
% D. Rudnick, 29 April 2021 - added vars as an argument

% vars={'t','s','udop','vdop','fl','abs'};
depthmin=200;

ctd=combineMissionsxy(year,vars,depthmin);
%add udopalong, udopacross
ctd=rotateuv(ctd,ctd.x0,ctd.y0,ctd.x1,ctd.y1);

A=meanfit(ctd,year,vars);






