%% function
% Compute the proton beam energy after the range shifter
% Use the CSDA approximation (CSDArange.m) to make the range estimates
%
%% Syntax
% |Energy = energyAfterRangeShifter(Beams)|
%
%
%% Description
% |Energy = energyAfterRangeShifter(Beams)| Description
%
%
%% Input arguments
% |material| - _STRING_ - Name of the material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%
% |E| -_SCALAR_- Energy (MeV) of the incoming proton beam
%
% |thick| -_SCALAR_- Thickness (mm) of the range shifter
%
%
%% Output arguments
%
% |Eout| - _SCALAR_ - Energy (MeV) after the range shifter
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Eout = energyAfterRangeShifter(material , Ein , thick)

  [range , Ep] = CSDArange(material); %Range (cm) in the material vs energy

  Rmax = interp1(Ep , range , Ein); %Range (cm) of the incoming particle
  Rmin = Rmax - thick/10;  %Range (cm) of the out-going particle

  Eout = interp1(range, Ep , Rmin); %Energy of the outgoing particle


end
