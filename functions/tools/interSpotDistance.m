%% interSpotDistance
% Compute the distance between spots
%
%% Syntax
% |DistMat = interSpotDistance(spot) |
%
%
%% Description
% |DistMat = interSpotDistance(spot) | Description
%
%
%% Input arguments
% |spot| - _SCALAR MATRIX_ - |spot(i) = [x,y]| x,y position of the i-th spot
%
%
%% Output arguments
%
% |DistMat| - _SCALAR MATRIX_ - DistMat(i,j) distance between the i-th and j-th spot
%
% |Dx| - _SCALAR MATRIX_ - Dx(i,j) distance along the X axis between the i-th and j-th spot
%
% |Dy| - _SCALAR MATRIX_ - Dy(i,j) distance along the Y axis between the i-th and j-th spot
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [DistMat , Dx , Dy]= interSpotDistance(spot)

spotX = spot(:,1);
spotY = spot(:,2);

X1 = repmat(spotX,1,length(spotX));
Y1 = repmat(spotY,1,length(spotY));

X2 = X1';
Y2 = Y1';

DistMat = sqrt((X2-X1).^2 + (Y2-Y1).^2);
Dx = sqrt((X2-X1).^2);
Dy = sqrt((Y2-Y1).^2);

end
