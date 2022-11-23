%% makeListOfCoordinates
% Create a list with the (x,y) coordinates of all the points in a square matrix when
% the coordinate points on the X and Y axis are given
%
%% Syntax
% |res = help_header(im1,im2)|
%
%
%% Description
% |res = help_header(im1,im2)| Description
%
%
%% Input arguments
% |x_axis| - _SCALAR VECTOR_ - |x_axis(x)| X coordinate of the x-th point in the square matrix
%
% |y_axis| - _SCALAR VECTOR_ - |y_axis(y)| Y coordinate of the y-th point in the square matrix
%
%% Output arguments
%
% |pts| - _SCALAR MATRIX_ -  |pts(i,:)=[x,y]| Coordinates (as defined by |x_axis| and |y_axis|) of the i-th point in the grid
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [pts , Xc , Yc] = makeListOfCoordinates(x_axis , y_axis)

  [Yc,Xc] = meshgrid(y_axis , x_axis);
  pts = [Xc(:) , Yc(:)]; %Linearise Xc; convert it into a column vector; placed column vector Yc necxt to it.

end
