%% getDistances2
% Compute the square of the distance between a point and each points in a list of points
%
%% Syntax
% |Dist2 = getDistances2(ListOfPoints , RefPoint)|
%
%
%% Description
% |Dist2 = getDistances2(ListOfPoints , RefPoint)| Description
%
%
%% Input arguments
% |ListOfPoints(i,:)=[x,y]| Position iof the i-th point in the list
%
% |RefPoint=[x,y]| Position of the point
%
%
%% Output arguments
%
% |Dist2(i)| Distance square of the i-th point in the list to the new point
%
% |idx| -_SCALAR VECTOR_- Indexes of all the points closest to RefPoint
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Dist2 , idx]= getDistances2(ListOfPoints , RefPoint)
  if (length(ListOfPoints)==0)
    Dist2 = [];
    return;
  end
  Dist2 = (ListOfPoints(:,1) - RefPoint(1)).^2 + (ListOfPoints(:,2) - RefPoint(2) ).^2;
  minD = min(Dist2); %find the distance to the spots closest to P1
  idx = find(Dist2 == minD); %Find the indices of the spots closest to P1
end
