%% isPtinGrd
% Determine whether a point is within the limits of a rectangular grid
%
%% Syntax
% |flag = isPtinGrd(Pt,GridSize)|
%
%
%% Description
% |flag = isPtinGrd(Pt,GridSize)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y]
%
% |GridSize| -_SCALAR MATRIX_- Dimension of the grid on which the spot are placed |GridSize(:,1) = [minX , minY]| |GridSize(:,2) = [maxX , maxY]|
%
%% Output arguments
%
% |flag| - _BOOLEAN_ - True = Pt is inside grid. false = Pt is outside grid
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function flag = isPtinGrd(Pt,GridSize)
  flag = Pt(1)>= GridSize(1,1) && Pt(2)>=GridSize(1,2) && Pt(1)<=GridSize(2,1) && Pt(2)<=GridSize(2,2);

end
