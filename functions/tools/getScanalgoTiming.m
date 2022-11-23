%% getScanalgoTiming
% Get the spot timing from a call to the scanalgo gateway.
%
%% Syntax
% |[Tsweep , TimePerSpot] = getScanalgoTiming(spot , w , Energy , scanAlgoGW)|
%
%
%% Description
% |[Tsweep , TimePerSpot] = getScanalgoTiming(spot , w , Energy , scanAlgoGW)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y] Coordinate units: mm in the IC gantry CS
%
% |Energy| -_SCALAR_- Energy (MeV) of the spot layer
%
% |scanAlgoGW| -_STRUCTURE_- Information about the scanAlgo gateway
%    * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%    * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%    * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%    * |scanAlgoGW.spot_id| -_STRING_-  beam supply point as defined in the site config jar site.properties
%
% |w| -_SCLAR VECTOR_- [OPTIONAL. Used only to compute |TimePerSpot|] |w(j)| the weight of the j-th spot.
%
% |aggregateSpots| -_BOOL_- [OPTIONAL: default = false] If true, the merge the PBS spots having the same X,Y position.
%                       If MU is too large, scanAlgo divides a spot into multiple spots.
%
%
%% Output arguments
%
% |Tsweep| -_SCALAR VECTOR_- |Tsweep(i)| Sweep time (ms) to go from the (i-1)-th spot to the i-th spot. The first index corresponds to the second spot
%                         |numel(Tsweep) = numel(TimePerSpot)-1|
%
% |TimePerSpot| -_SCALAR VECTOR_- |TimePerSpot(s)| Time (ms) required to deliver all the protons in s-th spot of the seqeunce
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Tsweep , TimePerSpot] = getScanalgoTiming(spot , Energy , scanAlgoGW , w , aggregateSpots)

  if nargin < 4
    %No spot weight provided. Use arbitrary spot weight
    Plan{1}.spots.weight = ones(size(spot,1),1);
  else
    %spot weight provided. Use proper weight in scanAlgo
    if isempty(w)
      Plan{1}.spots.weight = ones(size(spot,1),1);
    else
      Plan{1}.spots.weight = w;
    end
  end

  if nargin < 5
    aggregateSpots = false;
  end

  %scanalgo gateway parameters
  nb_paintings = 1;

  %Create a treatment plan with the provided spots
  Plan{1}.spots.nb_paintings = nb_paintings;
  Plan{1}.spots.energy = Energy;
  Plan{1}.spots.xy = spot(:,1:2); %pbs_convert_ScanAlgo takes IEC gantry coordinates

  %Call scanAlgo gateway to obtain the spot timing
  delivery = pbs_convert_ScanAlgo(Plan,{'nb_paintings',nb_paintings,'gateway_IP',scanAlgoGW.scanalgoGateway_IP,'room_id',scanAlgoGW.room_id,'spot_id',scanAlgoGW.spot_id,'snout_id',scanAlgoGW.snout_id,'sortSpot','false'} , false);

  %If the MU of a spot is too large, then scanAlgo will split the spot in multiple spots
  %Let's aggregate the spot back together, based on the spot position
  if aggregateSpots
    delivery = aggregate_PBS_paintings(delivery, false , true);
  end

  %compute the sweep time between spots
  Tsweep = getSweepTime(delivery{1}.spots.time, delivery{1}.spots.duration);

  if nargin < 3
    %No weight were provided. We cannot compute sensible spot delivery time
    TimePerSpot = [];
  else
    TimePerSpot = delivery{1}.spots.duration' ;
  end

end
