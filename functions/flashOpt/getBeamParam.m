%% getBeamParam
% Description
%
%% Syntax
% |beam=getBeamParam(MACHINE,beam)|
%
%
%% Description
% |beam=getBeamParam(MACHINE,beam)| Description
%
%
%% Input arguments
% |MACHINE| - _STRING_ - Label of the requested beam delivery
%
% |beam| -_STRUCTURE_- Structure to which the beam delivery informatio nwill be added
%
%
%% Output arguments
%
% |beam| -_STRUCTURE_- Updated beam delivery structure
%
%
%% Contributors
% Authors : R.Labarbe (open.reggui@gmail.com)

function beam=getBeamParam(MACHINE,beam,material)

    switch MACHINE
      case 'Pplus'
        %Machine specification for C230
        beam.MINspotTime = 1; %ms Shortest time in which a spot can be delivered
        beam.MAXcurrent = 0.5; %uA Maximum current available in treatment room
        beam.spotSweepTime = 1;% ms
        beam.lineSweepTime = 1; %ms Time to move spot back to begining of next line

        beam.current = 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 16; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,

        beam.R0 = 32.5;%  _SCALAR_ -  Range (cm) -	Do a “set range” at the maximal avilable range (32.5)
        beam.SpotSigma = 0.3; %cm Spot sigma 6-7mm at isocentre
        beam.RangeShifter = 0; %No Conformal Energy Filter. One single layer
        beam.RSweight = 1;


      case 'Pone'
        %Machine specification for S2C2
        beam.MINspotTime = 3; %ms Shortest time in which a spot can be delivered. 3 painting at 1ms => 3ms
        beam.MAXcurrent = 0.04; %uA Maximum current available in treatment room

        beam.spotSweepTime = 0;% ms If the motion is less than 8mm, then the otion is done between 2 pulses
        beam.lineSweepTime = 0; %ms Time to move spot back to begining of next line
        beam.current = 0.034; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. =0.07 uA -	Optics + slits settings were adjusted to maximize the transmission (based on simuations done by Quentin and his team), transmission was around 30-35%
        beam.SpotDuration = 10; %ms Time to deliver one spot. Constraint > 1ms


      case 'Pone-Oxford'
        %Reading experiment
        E = 230; %MeV
        beam.R0 =  material.alpha .* E^material.p;%  _SCALAR_ -  Range (cm) % E = (R./alpha).^(1./p) => alpha . .* E^p = R
        beam.MINspotTime = 1; %ms -	Shoot the map with the BDCU service screens -> only one “burst”, no dose control
        beam.MAXcurrent = 1; %=0.07 uA -	Do not worry about the current
        beam.SpotSigma = 0.3; %cm Spot sigma 6-7mm at isocentre
        beam.RangeShifter = 0; %No Conformal Energy Filter. One single layer
        beam.RSweight = 1;
        beam.MaxOffsetX = 1; %cm
        beam.MaxOffsetY = 1; %cm
        beam.spacing = 1.5; %Spacing between spot, in multiple of spot sigma

        beam.spotSweepTime = 0;% ms If the motion is less than 8mm, then the otion is done between 2 pulses
        beam.lineSweepTime = 0; %ms Time to move spot back to begining of next line
        Nbpulses = 4;
        beam.current =   2 .* 0.0705 .* 3 ./ Nbpulses; %Average beam current in uA in nozzle
        beam.SpotDuration = Nbpulses; %ms Time to deliver one spot. There were 6 pulses per spot

        beam.RidgeFilterFunction = @ridgeFilterVector;

      case 'Pone-Reading'
        %Reading experiment
        E = 230; %MeV
        beam.R0 =  material.alpha .* E^material.p;%  _SCALAR_ -  Range (cm) % E = (R./alpha).^(1./p) => alpha . .* E^p = R
        beam.MINspotTime = 1; %ms -	Shoot the map with the BDCU service screens -> only one “burst”, no dose control
        beam.MAXcurrent = 1; %=0.07 uA -	Do not worry about the current
        beam.SpotSigma = 0.3; %cm Spot sigma 6-7mm at isocentre
        beam.RangeShifter = 0; %No Conformal Energy Filter. One single layer
        beam.RSweight = 1;

        beam.spotSweepTime = 0;% ms If the motion is less than 8mm, then the otion is done between 2 pulses
        beam.lineSweepTime = 0; %ms Time to move spot back to begining of next line
        beam.current = 0.235 .* 0.35; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. =0.07 uA -	Optics + slits settings were adjusted to maximize the transmission (based on simuations done by Quentin and his team), transmission was around 30-35%
        %According to ODE: courant approximatif en nozzle est de 80 pC/pulse => 80e-12Cb * 6pulses / 6e-3s = 0.08uA
        beam.SpotDuration = 6; %ms Time to deliver one spot. There were 6 pulses per spot

      case 'Pone-Reading-large'
        %Reading experiment
        E = 230; %MeV
        beam.R0 =  material.alpha .* E^material.p;%  _SCALAR_ -  Range (cm) % E = (R./alpha).^(1./p) => alpha . .* E^p = R
        beam.MINspotTime = 1; %ms -	Shoot the map with the BDCU service screens -> only one “burst”, no dose control
        beam.MAXcurrent = 1; %=0.07 uA -	Do not worry about the current
        beam.SpotSigma = 0.3; %cm Spot sigma 6-7mm at isocentre
        beam.RangeShifter = 0; %No Conformal Energy Filter. One single layer
        beam.RSweight = 1;

        beam.spotSweepTime = 0;% ms If the motion is less than 8mm, then the otion is done between 2 pulses
        beam.lineSweepTime = 0; %ms Time to move spot back to begining of next line
        beam.current = 0.235 .* 0.35; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. =0.07 uA -	Optics + slits settings were adjusted to maximize the transmission (based on simuations done by Quentin and his team), transmission was around 30-35%
        %According to ODE: courant approximatif en nozzle est de 80 pC/pulse => 80e-12Cb * 6pulses / 6e-3s = 0.08uA
        beam.SpotDuration = 6; %ms Time to deliver one spot. There were 6 pulses per spot

      case 'Pone-Reading-plus'
        %Reading experiment
        E = 230; %MeV
        beam.R0 =  material.alpha .* E^material.p;%  _SCALAR_ -  Range (cm) % E = (R./alpha).^(1./p) => alpha . .* E^p = R
        beam.MINspotTime = 1; %ms -	Shoot the map with the BDCU service screens -> only one “burst”, no dose control
        beam.MAXcurrent = 1; %=0.07 uA -	Do not worry about the current
        beam.SpotSigma = 0.3; %cm Spot sigma 6-7mm at isocentre
        beam.RangeShifter = 0; %No Conformal Energy Filter. One single layer
        beam.RSweight = 1;

        beam.spotSweepTime = 0;% ms If the motion is less than 8mm, then the otion is done between 2 pulses
        beam.lineSweepTime = 0; %ms Time to move spot back to begining of next line
        %ODE: En pratique on sort sans trop de problème 130 pC/pulse des machines => 130 nA moyen et on sait qu’on sait sortir un peu plus (ce qu’on a fait pour le test FLASH) mais alors plus en conditions « cliniques ».
        beam.current = 0.27; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. =0.07 uA -	Optics + slits settings were adjusted to maximize the transmission (based on simuations done by Quentin and his team), transmission was around 30-35%
        beam.SpotDuration = 1; %ms Time to deliver one spot.

      case 'Pone-Reading-small'
        %Reading experiment
        E = 230; %MeV
        beam.R0 =  material.alpha .* E^material.p;%  _SCALAR_ -  Range (cm) % E = (R./alpha).^(1./p) => alpha . .* E^p = R
        beam.MINspotTime = 1; %ms -	Shoot the map with the BDCU service screens -> only one “burst”, no dose control
        beam.MAXcurrent = 1; %=0.07 uA -	Do not worry about the current
        beam.SpotSigma = 0.3; %cm Spot sigma 6-7mm at isocentre
        beam.RangeShifter = 0; %No Conformal Energy Filter. One single layer
        beam.RSweight = 1;

        beam.spotSweepTime = 0;% ms If the motion is less than 8mm, then the otion is done between 2 pulses
        beam.lineSweepTime = 0; %ms Time to move spot back to begining of next line
        beam.current = 0.235 .* 0.3; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. =0.07 uA -	Optics + slits settings were adjusted to maximize the transmission (based on simuations done by Quentin and his team), transmission was around 30-35%
        %According to ODE: courant approximatif en nozzle est de 80 pC/pulse => 80e-12Cb * 6pulses / 6e-3s = 0.08uA
        beam.SpotDuration = 24; %ms Time to deliver one spot. There were 6 pulses per spot

      case 'Pplus-Groningen'
        beam = getGroningenBeamParam(beam);

      case 'C230option1'
        beam = getGroningenBeamParam(beam);
        beam.current = 3 .* 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 16 ./3; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,

      case 'C230option2'
        beam = getGroningenBeamParam(beam);
        beam.current = 3 .* 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 1; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,

      case 'C230option3'
        beam = getGroningenBeamParam(beam);
        beam.current = 3 .* 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,

      case 'C230option4'
        beam = getGroningenBeamParam(beam);
        beam.current = 2 .* 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,

      case 'C230option5'
        beam = getGroningenBeamParam(beam);
        beam.current = 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
      case '0.5cm2'
        beam = getGroningenBeamParam(beam);
        beam.current = 5.*0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 2.5.*1.0135; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
          %Normalise the dose to be equal to that of 3cm2
      case '1cm2'
        beam = getGroningenBeamParam(beam);
        %beam.current = 0.001.*0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        %beam.SpotDuration = 1000.*2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
        beam.current = 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,

      case '1cm2Ip'
        beam = getGroningenBeamParam(beam);
        beam.current = 10.*3.* 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 0.1.*2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
      case '1cm2Fast'
        beam = getGroningenBeamParam(beam);
        beam.current = 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 2.5.*1.0135; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
        beam.spotSweepTime = 0;% ms
        beam.lineSweepTime = 0; %ms Time to move spot back to begining of next line

      case '2cm2'
        beam = getGroningenBeamParam(beam);
        %beam.current = 0.001.*0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        %beam.SpotDuration = 1000.*2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
        beam.current = 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,


      case '3cm2'
        beam = getGroningenBeamParam(beam);
        %beam.current = 0.001.*0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        %beam.SpotDuration = 1000.*2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
        beam.current = 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,

      case 'scanned'
        beam = getGroningenBeamParam(beam);
        %beam.current = 0.001.*0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        %beam.SpotDuration = 1000.*2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
        beam.current = 0.043 * 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 10; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
        beam.MaxOffsetX = 5; %cm
        beam.MaxOffsetY = 5; %cm

      otherwise
        beam = getGroningenBeamParam(beam);
        beam.current = 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
        beam.SpotDuration = 2.5; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,
    end
end




function beam = getGroningenBeamParam(beam)

  %Machine specification for C230
  beam.MINspotTime = 1; %ms Shortest time in which a spot can be delivered
  %beam.MAXcurrent = 0.5; %uA Maximum current available in treatment room
  beam.MAXcurrent = 5; %uA Maximum current available in treatment room
  beam.spotSweepTime = 1;% ms
  beam.lineSweepTime = 1; %ms Time to move spot back to begining of next line

  beam.current = 0.3 .* 0.7; %Average beam current in uA in nozzle for the spot with beam.RSweight =1. Line transmission is 70%
  beam.SpotDuration = 16; %ms Time to deliver one spot. modify the maximum duration of the spots so that the map is “time driven”,

  beam.R0 = 32.5;%  _SCALAR_ -  Range (cm) -	Do a “set range” at the maximal avilable range (32.5)
  beam.SpotSigma = 0.3; %cm Spot sigma 6-7mm at isocentre
  beam.RangeShifter = 0; %No Conformal Energy Filter. One single layer
  beam.RSweight = 1;

  beam.RidgeFilterFunction = @ridgeFilterStep;
end

%================================
% Conformal Energy Filter Definition
%
% H = ridgeFilter() : return the height of the Conformal Energy Filter
% weight = ridgeFilter(Rshift) : Fraction of proton undergoing that range shift
%
% The Conformal Energy Filter is a cone with a base radius r0 and a height H.
% The number of proton with a range shift deltaR is proportional to the perimeter of the circle at height delta-R (coutned from basis)
% The perimeter are normalised to that the integral from basis to top is equal to 1, so as to spread the proton in volume
%
% Input
%
% |Rshift| _-SCALAR-_ Range shift (cm)
%
% Ouput
%
% |weight| _-SCALAR-_ Fraction of proton undergoing that range shift
%================================
function weight = ridgeFilter(Rshift)

%Description of the conical Conformal Energy Filter
r0 = 0.7; %cm Radius of the base of the cylinder spike
H = 10; %cm Height of the Conformal Energy Filter

if (nargin <1)
    % Return the height of the Conformal Energy Filter
    %======================================
    weight = H;
  else
    %Return the weight for this range shift
    %======================================
    %Compute the radius of the cone at the specified thickness of Ridge
    r = r0 - (r0./H) .* Rshift; % Radius of the cone at the height of the range shift

    %The weight is the ratio of the area at height Rshift to the basis area
    weight = r; %Normalise so that the volume of the cone is 1 (the proton are spread in depth without loosing protons)
  end
end

function weight = ridgeFilterVector(Rshift)

  RidgeFilterWeight = [2.5 0.1 0.1  0.5 0.1];
  RidgeFilterShift =  [0 0.5  1   1.5  2]; %cm WET

  if (nargin <1)
      % Return the height of the Conformal Energy Filter
      %======================================
      weight = max(RidgeFilterShift);
    else
      weight = interp1(RidgeFilterShift , RidgeFilterWeight , Rshift , [],0);
  end
end
