%% scan2Neighbor
% Order the spots so as to scan from one spot to one of its neighbourgh.
% Start at the spot with the least neightbourgh and then move From
% closest neighbourgh to closest neightbourgh
%
%% Syntax
% |OrderedSpot =  scan2Neighbor(spot , sigmaAtIso , BDL)|
%
%
%% Description
% |OrderedSpot =  scan2Neighbor(spot , sigmaAtIso , BDL)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y]
%
% |sigmaAtIso| - _SCLAR_ - Sigma (m) of tha Gaussian beam profile of one PBS spot
%
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
%
%% Output arguments
%
% |OrderedSpot| - _SCLAR MATRIX_ - The |spot| in the ordered sequence. i-th spot to deliver is OrderedSpot(i,:) = [x,y]
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function OrderedSpot =  scan2Neighbor(spot , sigmaAtIso , BDL)

  spotSequence = [];
  LeftSpots = 1:size(spot,1); %List of the indices of spots not yet selected

  Nmaps = getTopologicalMaps(spot , BDL , sigmaAtIso); %Get the topological maps
  Neighbourgh = Nmaps.NeighbourghMap;
  [~ , idxSpot] = min(sum(Neighbourgh,1)); %Find the spot with least neightbough
  spotSequence = idxSpot; %The first spot in the sequence is the spot with least neightbourgh
  LeftSpots(idxSpot)=[]; %Remove the selected spot from the list of remaining spots

  while (~isempty(LeftSpots))
      %Loop untill all spots are added to the sequence
      [~, idxNext] = min(Nmaps.NeighbourghTimeMap(idxSpot,LeftSpots)); %find the closest neighbourgh
      spotSequence = [spotSequence , LeftSpots(idxNext)];
      idxSpot = LeftSpots(idxNext);
      LeftSpots(idxNext)=[]; %Remove the selected spot from the list of remaining spots
  end

OrderedSpot = spot(spotSequence,:); %The ordered spot sequence

end
