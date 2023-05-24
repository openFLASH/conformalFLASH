%% WETFromThickness
% Compute the WET of a given thickness of material
%
%% Syntax
% |[RangeShifterWET , EoutRS]= WETFromThickness(RangeShifterMaterial, E, thickness)|
%
%% Description
% |[RangeShifterWET , EoutRS]= WETFromThickness(RangeShifterMaterial, E, thickness)| - Compute the WET of a given thickness of material
%
%% Input arguments
% |material| - _STRING_ - Name of the material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
% |E| -_SCALAR VECTOR_- Energy (MeV) of the incoming proton
% |thickness| -_SCALAR VECTOR_- Thickness (mm) of the range shifting material%
%
%
%% Output arguments
% |RangeShifterWET| -_SCALAR VECTOR_- Water equivalent thickness (cm) of the corrsponding thickness of material
% |EoutRS| -_SCALAR VECTOR_- Energy (MeV) of the proton beam coming out of the corresponding thicknes of material


%% Contributors
% Authors : L. Hotoiu, R. Labarbe (open.reggui@gmail.com)


function [RangeShifterWET , EoutRS]= WETFromThickness(RangeShifterMaterial , E , thickness)

  water = materialDescription('water');
  R_beam = energy2range(E, water.alpha,water.p); %Range (cm) in water of the proton delivered by the machine
  EoutRS = energyAfterRangeShifter(RangeShifterMaterial , E , thickness); %Proton energy out of the range shifter
  R_downRS = energy2range(EoutRS, water.alpha,water.p); %Range (cm) in water of the beam coming out of the range shifter material
  RangeShifterWET = (R_beam - R_downRS); %WET of the range shifter
end
