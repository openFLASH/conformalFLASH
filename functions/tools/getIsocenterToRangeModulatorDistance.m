%% getIsocenterToRangeModulatorDistance
% Compute the distance from isocentre to the base of the CEF
% The CEF is upstream to the range shifter
% |IsocenterToRangeModulatorDistance| is the distance to the upstream surface of the CEM.
%
%% Syntax
% |[IsocenterToRangeModulatorDistance  , ModulatorMountingPosition] = getIsocenterToRangeModulatorDistance(Beams, BDL)|
%
%
%% Description
% |[IsocenterToRangeModulatorDistance  , ModulatorMountingPosition] = getIsocenterToRangeModulatorDistance(Beams, BDL)| Description
%
%
%% Input arguments
% |Beams| - _STRUCT_ - Information about beam structure
%     * |Beam.RangeModulator.CEMThicknessData| -_SCALAR MATRIX_- [OPTIONAL] |CompensatorThicknessData(x,y)| Thickness (mm) of the CEF pixel at position (x;y) in the IEC beam Limiting device CS
%     * |Beams.RSinfo.RangeShifterMaterial| - _STRING_ - Name of the range shifter material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%     * |Beams.RSinfo.RangeShifterWET| -_SCALAR_- Water equivalent thickness (mm) of the range shifter to add to the CEF
%     * |Beams.Iso2Skin| -_SCALAR_- Isocentre (mm) To Skin Distance along the proton beam axis
%     * |Beams.Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
%
%% Output arguments
%
% |IsocenterToRangeModulatorDistance| -_SCALAR_- Isocenter to UPSTERAM edge of CEF (mm).
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function IsocenterToRangeModulatorDistance  = getIsocenterToRangeModulatorDistance(Beams, BDL)

  IsocenterToBlockTrayDistance = getIsocenterToBlockTrayDistance(Beams);
  param = getMachineParam(BDL);
  IsocenterToRangeModulatorDistance = round(IsocenterToBlockTrayDistance + param.snout.CEMOffset,1);


end
