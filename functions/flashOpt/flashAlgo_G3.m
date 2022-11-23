%% flashAlgo_G3
% Compute the optimum sequence of spot in order to achieve high average dose rate
% * Identify the major axis of the scanned area
% * Set the main axis of the scarves paralell to the major axis
% * Place multiple scarves next to each other, aligned with the major axis
% * The number of paralell scarves is qual to the minor axis length divided by the width of a scarf
%
%% Syntax
% |[spotSequence , Nmaps] = flashAlgo_G3(spot , sigmaAtIso , BDL , DRcritical , NbScarves , plotNb , scanAlgoGW)|
%
%
%% Description
% |[spotSequence , Nmaps] = flashAlgo_G3(spot , sigmaAtIso , BDL , DRcritical , NbScarves , plotNb , scanAlgoGW)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y]
%
% |sigmaAtIso| - _SCLAR_ - Sigma (m) of tha Gaussian beam profile of one PBS spot
%
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
% |DRcritical| -_SCALAR VECTOR_- Vector with the bemalet number of all the beamlets that hit a structure for which there is a dose rate constraint
%
% |NbScarves| -_NTEGER_- Number of scarve to lay paralell to the main ellipse axis. Each scarf will be scanned sequentially
%
% |plotNb| -_SCALAR_- [OTPIONAL]  |plotNb(b)| Figure number where to plot the trajectory of beam b
%
% |scanAlgoGW| -_STRUCTURE_- [OPTIONAL. Only required to use a scanAlgo Gateway to get the timing] Information about the scanAlgo gateway
%    * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%    * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%    * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%    * |scanAlgoGW.spot_id| -_STRING_-  beam supply point as defined in the site config jar site.properties
%
%% Output arguments
%
% |spotSequence| -_SCALAR VECTOR_- Order of the indices of |spot| to sort the spots. |OrderedSpot = spot(spotSequence,:)|
%
%  |Nmaps| -_STRUCTURE_- Topological information on the initial spot sequence
% * |Nmaps.NeighbourghMap| -_SCALAR VECTOR_- |NeighbourghMap(d,i)| d=# of delivered spot; i=# of impacted spot; |NeighbourghMap(d,i)|  = 0 if there is an impact
% * |Nmaps.NeighbourghWeightMap| -_SCALAR VECTOR_- |NeighbourghWeightMap(d,i)| fraction of the dose of spot d (# of delivered spot) that is also delivered at spot i (i=# of impacted spot)
% * |Nmaps.NeighbourghTimeMap| -_SCALAR MATRIX_- |NeighbourghTimeMap(d,i)| Time (ms) required to go from spot d to spot i
%
%% REFERENCE
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [spotSequence , Nmaps] = flashAlgo_G3(spot , sigmaAtIso , BDL , DRcritical , NbScarves , plotNb , scanAlgoGW)

  if(nargin < 6)
    plotNb = [];
  end

  if(nargin < 7)
    scanAlgoGW = [];
  end

  %Add an index to each spot, so we can trace the sorting at the end
  indx = 1:size(spot,1);
  spot = [spot , indx'];
  param = getMachineParam(BDL);
  Nmaps = getTopologicalMaps(spot , BDL , sigmaAtIso , scanAlgoGW); %Get the inital topological map
  [Vx , Vy] = getScanSpeed(spot , BDL, scanAlgoGW); %Get the scan speed along each axis

  %First, optimize trajectory in the structures requiring minimum dose rate
  %------------------------------------------------------------------------
  %The first spot in the sequence is the one further away from isocentre.
  %This will force the scanning to start from the remotest edge of the OAR
  OARspots = spot(DRcritical,1:2);
  [~ , OrderedCoord] = orderALongScarves(OARspots, [Vx , Vy] , NbScarves , plotNb);

  OrderedSpot1 = spot(DRcritical(OrderedCoord),:);
  spotSequence1 = OrderedSpot1(:,3);

  %Then optimize the rest of the spots
  %----------------------------------
  SpotIndices = 1:size(spot,1); %The lsit of all the spots
  SpotIndices(DRcritical) = []; %Remove the spots that are already ordered
  if (~isempty(SpotIndices))
    OrderedSpot2 = scan2Neighbor(spot(SpotIndices,:) , sigmaAtIso , BDL); %Order the remaining spots
    spotSequence2 = OrderedSpot2(:,3); %Get the indice or the newly ordered spots
  else
    %All the spots are in a critical structure
    spotSequence2 = [];
  end

  %Compiled sequence
  spotSequence =  [spotSequence1 ; spotSequence2]; %Collect the indices of all the sorted spots
  OrderedSpot = spot(spotSequence,:); %The ordered spot sequence

end
