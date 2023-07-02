% plot output from example 2, to verify it lies on the surface of an
% ellipsoid
%
% created june 17, 2022 
%

clear;

datafile = 'example2.txt';

xall = load(datafile);
nx = size(xall,1);

scatter3(xall(:,1),xall(:,2),xall(:,3));