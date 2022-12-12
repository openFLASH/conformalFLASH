%% getScanSpeed
% Compute the magnet scanning speed along the X and Y axis
% If |scanAlgoGW|  is provided, use a call to scanAlgo gateway to get info on spot timing.
% Otherwise, use a simple model for the scan speed provided by getMachineParam.m
%
%% Syntax
% |[Vx , Vy] = getScanSpeed(spot , BDL , scanAlgoGW)|
%
%
%% Description
% |[Vx , Vy] = getScanSpeed(spot , BDL , scanAlgoGW)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y] Coordinate units: mm
%
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
% |scanAlgoGW| -_STRUCTURE_- [OPTIONAL] Information about the scanAlgo gateway
%    * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%    * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%    * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%    * |scanAlgoGW.spot_id| -_STRING_-  beam supply point as defined in the site config jar site.properties
%
%% Output arguments
%
% |Vx| - _SCALAR_ - Speed (m/s) of displacement of the spot along the X axis
%
% |Vy| - _SCALAR_ - Speed (m/s) of displacement of the spot along the Y axis
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Vx , Vy] = getScanSpeed(spot , BDL , scanAlgoGW)

  if nargin < 3
    scanAlgoGW = [];
  end

  param = getMachineParam(BDL);

  if ~isempty(scanAlgoGW)

      %Get timing from scanAlgo
      fprintf('Getting scan speed data from %s \n', scanAlgoGW.scanalgoGateway_IP)
      step = 15; %mm Arbitrary distance between spots

      % Scanning Y axis
      %========================
      spot = [];
      X=0;
      for Y = -10:1:10
          spot = [spot ; [X,Y].*step];
      end
      Nmaps = getTopologicalMaps(spot , BDL , 3 , scanAlgoGW); %Get the sweep time. The spot sigma is irrelevant as we are only interested in spot timing, not beighbourhood
      DistMat  = interSpotDistance(spot); %Get the spot distance

      dT = Nmaps.NeighbourghTimeMap(:);  %mm
      dX =   DistMat(:); %ms
      Vy = dT\dX; % X\Y;
      fprintf('Y Scanning speed D/T : %f m/s \n', Vy)

      % Scanning X axis
      %========================
      spot = [];
      Y=0;

      for X = -8:1:8
        spot = [spot ; [X,Y].*step];
      end

      Nmaps = getTopologicalMaps(spot , BDL , 3, scanAlgoGW); %Get the sweep time
      DistMat  = interSpotDistance(spot); %Get the spot distance
      dT = Nmaps.NeighbourghTimeMap(:) ; %mm
      dX =   DistMat(:); %ms
      Vx = dT\dX; % X\Y
      fprintf('X Scanning speed D/T : %f m/s \n', Vx)

else

    %No scanAlgo available. Use same scan speed in X and Y
    fprintf('Getting default scan speed \n')
    Vx = param.ScanSpeed ./1000; %m/s
    Vy = param.ScanSpeed ./1000;
end

end
