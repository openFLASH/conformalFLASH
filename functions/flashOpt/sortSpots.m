%% sortSpots
% Sort the spots so that the neighbourgh spots are delivered at similar time
%
%% Syntax
% |sortedSpots = sortSpots(spot,Nmaps)|
%
%
%% Description
% |sortedSpots = sortSpots(spot,Nmaps)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y,z,w ,i]
%
% |Nmaps| -_STRUCTURE_-
% * |Nmaps.NeighbourghMap| -_SCALAR VECTOR_- |NeighbourghMap(d,i)| d=# of delivered spot; i=# of impacted spot; |NeighbourghMap(d,i)|  = 0 if there is an impact
% * |Nmaps.NeighbourghWeightMap| -_SCALAR VECTOR_- |NeighbourghWeightMap(d,i)| fraction of the dose of spot d (# of delivered spot) that is also delivered at spot i (i=# of impacted spot)
% * |Nmaps.NeighbourghTimeMap| -_SCALAR MATRIX_- |NeighbourghTimeMap(d,i)| Time (ms) required to go from spot d to spot i
%
%% Output arguments
% |sortedSpots| - _SCLAR MATRIX_ - Ordered spots. The i-th spot to deliver is spot(i,:) = [x,y,z,w ,i]
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function sortedSpots = sortSpots(spot,Nmaps)
  NbSpots = length(spot); %# of spots
  [~,posFirst]=max(Nmaps.NeighbourghMap); %Coordinate of the first neighbourgh
  [~,posLast]=max(flip(Nmaps.NeighbourghMap));
  posLast = NbSpots - posLast +1; %Coordinate of the last neighbourgh
  posAv = (posFirst + posLast)./2;

  %sort the sequence of the spots
  %[~,index] = sortrows([posAv ; posFirst ; posLast]'); %sort on the centre of the 1s
  [~,index] = sortrows([posFirst ; posLast]'); %sort on the position of first 1s

  sortedSpots = spot(index,:); %sort the spot in the chosen order
end
