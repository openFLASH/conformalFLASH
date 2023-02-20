%% getTopologicalMaps
% Compute the matrices describing the neighbourhood of the PBS spots
%
%% Syntax
% |Nmaps = getTopologicalMaps(spot , BDL , sigmaAtIso , scanAlgoGW)|
%
%
%% Description
% |Nmaps = getTopologicalMaps(spot , BDL , sigmaAtIso , scanAlgoGW)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y,w] Coordinate units: mm
%
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
% |sigmaAtIso| - _SCLAR_ - Sigma (m) of tha Gaussian beam profile of one PBS spot
%
% |scanAlgoGW| -_STRUCTURE_- [OPTIONAL. Force the use of the sscanAlgo gateway] Information about the scanAlgo gateway
%    * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%    * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%    * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%    * |scanAlgoGW.spot_id| -_STRING_- Spot tune ID as defined in the gateway
%
%% Output arguments
% |Nmaps| -_STRUCTURE_-
% * |Nmaps.NeighbourghMap| -_SCALAR VECTOR_- |NeighbourghMap(d,i)| d=# of delivered spot; i=# of impacted spot; |NeighbourghMap(d,i)|  = 0 if there is an impact
% * |Nmaps.NeighbourghWeightMap| -_SCALAR VECTOR_- |NeighbourghWeightMap(d,i)| fraction of the dose of spot d (# of delivered spot) that is also delivered at spot i (i=# of impacted spot)
% * |Nmaps.NeighbourghTimeMap| -_SCALAR MATRIX_- |NeighbourghTimeMap(d,i)| Time (ms) required to go from spot d to spot i
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Nmaps = getTopologicalMaps(spot , BDL , sigmaAtIso , scanAlgoGW)

  if nargin < 4
    scanAlgoGW = struct;
  end
  if isempty(scanAlgoGW)
    scanAlgoGW = struct;
  end

  if min(size(spot))==1
    %There is a single spot. Neighbourghiid is quite reduced
    Nmaps.NeighbourghMap = 1 ;% -_SCALAR VECTOR_- |NeighbourghMap(d,i)| d=# of delivered spot; i=# of impacted spot; |NeighbourghMap(d,i)|  = 0 if there is an impact
    Nmaps.NeighbourghWeightMap = 1 ;% -_SCALAR VECTOR_- |NeighbourghWeightMap(d,i)| fraction of the dose of spot d (# of delivered spot) that is also delivered at spot i (i=# of impacted spot)
    Nmaps.NeighbourghTimeMap = 0;% -_SCALAR MATRIX_- |NeighbourghTimeMap(d,i)| Time (ms) required to go from spot d to spot i
    return
  end

  DistMat = interSpotDistance(spot);
  NeighbourghWeightMap = exp(-(DistMat ./ sigmaAtIso).^2); % dose delivered by Gaussian tail to other locations
  Nmaps.NeighbourghMap = (NeighbourghWeightMap >= 0.001); %Define a neghboug only if it contributes to dose
  Nmaps.NeighbourghWeightMap = NeighbourghWeightMap .* Nmaps.NeighbourghMap; %Set to zero Gaussian tail < 0.1%
  %      NeighbourghWeightMap(d,i) d=# of delivered spot; i=# of impacted spot; NeighbourghWeightMap = 0 if there is an impact

  if isfield(scanAlgoGW , 'scanalgoGateway_IP')
    %Get timing from scnaAlgo
    fprintf('Getting timing from scanAlgo gateway: %s \n',scanAlgoGW.scanalgoGateway_IP)
    Nmaps.NeighbourghTimeMap = getTimeMapScanAlgo(spot , BDL, scanAlgoGW);

  else
    %Use an approximation of the timing and scanning speed
    param = getMachineParam(BDL);
    Nmaps.NeighbourghTimeMap = 1000 .* DistMat ./ param.ScanSpeed;
  end

end

%-------------------------------
function NeighbourghTimeMap = getTimeMapScanAlgo(spot , BDL, scanAlgoGW)

NbSpots = size(spot,1);
list = [];
param = getMachineParam(BDL);

if (NbSpots > 1)
      %There are more than 1 spot
      for stpSize = 1:NbSpots-1 %Step to jump
        for Si = 1:stpSize %index of starting spot
          for stp = Si:stpSize:NbSpots %jump from spot to spot
            list = [list ; spot(stp,1:2)];
          end %for stp
        end %for Si
      end %for stpMax

      %Get the sweep time from scanAlgo
      Tsweep1 = getScanalgoTiming(list , param.MAXenergy , scanAlgoGW , [] , false); %Timing for  motion A -> B
      list = flip(list,1);
      Tsweep2 = getScanalgoTiming(list , param.MAXenergy , scanAlgoGW , [] , false); %Timing for  motion B -> A
      Tsweep2 = Tsweep2(2:end);
      Tsweep2 = flip(Tsweep2,1);
      Tsweep2 = [0 ; Tsweep2];

      idx = 1;
      NeighbourghTimeMap = zeros(NbSpots,NbSpots);

      for stpSize = 1:NbSpots-1 %Step to jump
        for Si = 1:stpSize %index of starting spot
          idx = idx+1; %Skip the delta at the begining of a cycle
          for stp = Si+stpSize:stpSize:NbSpots
            %fprintf('%d : T(%d,%d) = %f ms and %f ms \n',idx,stp-1,stp,Tsweep1(idx), Tsweep2(idx))
            NeighbourghTimeMap(stp-stpSize,stp) = Tsweep1(idx);
            NeighbourghTimeMap(stp,stp-stpSize) = Tsweep2(idx);
            idx = idx+1;
          end %for stp
        end %for Si
      end %for stpMax
else
  %There is only 1 spot
  % NeighbourghTimeMap is a 1x1 matrix
  NeighbourghTimeMap =0;
end

end
