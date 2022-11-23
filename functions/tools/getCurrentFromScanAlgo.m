%% getCurrentFromScanAlgo
% Get an estimate of the nozzle current from scanAlgo.
% A set of spot with increasing MU are sent to scanAlgo. A linear regression on the spot duration
% is done in order to estimate the nozzle current
%
%% Syntax
% |Inozzle  = getCurrentFromScanAlgo(scanAlgoGW , BDL)|
%
%
%% Description
% |Inozzle  = getCurrentFromScanAlgo(scanAlgoGW , BDL)| Description
%
%
%% Input arguments
% |scanAlgoGW| -_STRUCTURE_- Information about the scanAlgo gateway
%    * |scanAlgoGW.scanalgoGateway_IP| -_STRING_- IP address, including port, to the scanAlkgo gatewat
%    * |scanAlgoGW.room_id| -_STRING_- Room ID as defined  inthe gateway
%    * |scanAlgoGW.snout_id|  -_STRING_- snout ID as defined in the gateway
%    * |scanAlgoGW.spot_id| -_STRING_-  beam supply point as defined in the site config jar site.properties
%
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
%% Output arguments
%
% |Inozzle| - _SCALAR_ - Average nozzle current (nA) estimated from scanAlgo
%
%
%% Contributors
% Authors : Rudi Labarbe (open.reggui@gmail.com)

function Inozzle  = getCurrentFromScanAlgo(scanAlgoGW , BDL)

  eV = 1.6021766208e-19 ; % J/eV
  param = getMachineParam(BDL);

  %Create a treatment plan with the provided spots
  nb_paintings = 1;
  Plan{1}.spots.nb_paintings = nb_paintings;
  Plan{1}.spots.energy = param.MAXenergy;

  %Create set of spot at arbitrary positions
  spot = [];
  X=0;
  step = 1;
  for Y = -10:1:10
      spot = [spot ; [X,Y].*step];
  end
  Plan{1}.spots.xy = spot; %pbs_convert_ScanAlgo takes IEC gantry coordinates
  Plan{1}.spots.weight = 1:size(spot,1);

  Charge = MU_to_NumProtons(Plan{1}.spots.weight, param.MAXenergy) .* eV; %Cb Charge per spot

  %Call scanAlgo gateway to obtain the spot timing
  delivery = pbs_convert_ScanAlgo(Plan,{'nb_paintings',nb_paintings,'gateway_IP',scanAlgoGW.scanalgoGateway_IP,'room_id',scanAlgoGW.room_id,'spot_id',scanAlgoGW.spot_id,'snout_id',scanAlgoGW.snout_id,'sortSpot','false'} , false);

  %If the MU of a spot is too large, then scanAlgo will split the spot in multiple spots
  %Let's aggregate the spot back together, based on the spot position
  delivery = aggregate_PBS_paintings(delivery, false , true);

  Inozzle  = 1e3 .* 1e9 .*  (delivery{1}.spots.duration \ Charge');   % X\Y; Nozzle current in nA


end
