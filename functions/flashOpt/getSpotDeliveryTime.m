%% getSpotDeliveryTime
% Estimate the time (s) to deliver a PBS spot
%
%% Syntax
% |T = getSpotDeliveryTime(w , Inozzle ,  BDL , Energy)|
%
%
%% Description
% |T = getSpotDeliveryTime(w , Inozzle ,  BDL , Energy)| Description
%
%
%% Input arguments
% |w| -_SCALAR VECTOR_- Weights (MU) to deliver the dose in the PBS spot
%
% |Inozzle| -_SCALAR_- Nozzle average current (A)
%
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
%
%% Output arguments
%
% |T| -_SCALAR VECTOR_- Time (s) to deliver the spot
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function T = getSpotDeliveryTime(w , Inozzle ,  BDL , Energy)

  physicsConstants;
  param = getMachineParam(BDL);

  if nargin < 4
    Charge = MU_to_NumProtons(w, param.MAXenergy) .* eV; %Convert MU into proton charge (Cb) in the spot
  else
    Charge = MU_to_NumProtons(w, Energy) .* eV; %Convert MU into proton charge (Cb) in the spot
  end

  if (round2nanoA(Inozzle) > round2nanoA(param.MAXcurrent .* 1e-6))
    %Requested current is too high. Set limit to high current
    warning(['Requested nozzle current is too high. Throttling current I = ',num2str(param.MAXcurrent),' uA'])
    Inozzle = round2nanoA(param.MAXcurrent .* 1e-6);
  end


  switch param.MachineType
    case 'PROTEUSplus'
      %The P+ has a continuous beam at averate current. The computation of delivery time is straightforward
      T = Charge ./ Inozzle; %Time (s) to deliver each spot. Assume nozzle current at 235MeV
      T(T < param.MinimumSpotDuration) = param.MinimumSpotDuration; %Make sure no spot is delivered in less than the minimum delivery time

    case 'PROTEUSone'
      T = ceil(Charge ./ Inozzle); %P1 delivers an integral number of pulses. So it must be integer
      T = T + 2; %ms

    otherwise
      T = Charge ./ Inozzle; %Time (s) to deliver each spot
  end

end

%===========================
% Round the current to the closest nano amp
%===========================
function A = round2nanoA(A)
    A = round(A .* 1e9) .* 1e-9;
end
