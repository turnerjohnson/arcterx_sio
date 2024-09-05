function model=meanmodel(A,vari,idepth,time,x,y)
% function meanmod applies the mean model in A to variable vari at depth
% idepth and time, x, y
%
% Dan Rudnick, 14 Feb 2020

timer=ut2dn(time)-datenum(A.year0,1,1);

model = A.(vari).constant(idepth) + A.(vari).time(idepth)*timer + A.(vari).x(idepth)*x + A.(vari).y(idepth)*y;
