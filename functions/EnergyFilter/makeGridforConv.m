%% makeGridforConv
% Create a square grid of equally space points extending up to radius0
% Typically the function is used to create the Cartesian coordinate grid
% for a Gaussian convolution kernel
%
%% Syntax
% |[x,y] = makeGridforConv(radius0,step)|
%
%
%% Description
% |[x,y] = makeGridforConv(radius0,step)| Description
%
%
%% Input arguments
% |radius0| - _SCALAR_ - The grid will extend from -radius0 to +radius0
%
% |step| -_SCALAR_- Distance between points in the grid
%
%
%% Output arguments
%
% |x| - SCALAR MATRIX_ - The X coordinate at point X(x,y)
%
% |y| - SCALAR MATRIX_ - The Y coordinate at point Y(x,y)
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [x,y] = makeGridforConv(radius0,step)

    %grid for gaussian convolution
    x = single(-radius0+step/2 : step : radius0-step/2);
    y = single(-radius0+step/2 : step : radius0-step/2);
    [x,y] = ndgrid(x,y);
end
