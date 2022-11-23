%% getIsocenterToRangeShifterDistance(Beams
% Compute the distance from isocentre to the DOWNSTREAM edge of range shifter
% The range shifter can be placed either:
% * on the skin surface if no aperture is present
% * on the aperture block if an aperture block is provided
%
% From DICOM standard : Isocenter to Range Shifter Distance (300A,0364)	Isocenter to DOWNSTREAM edge of range shifter (mm) at current control point.
% |IsocenterToRangeShifterDistance| is therefore the distance to the patient side of the range shifter

%
%% Syntax
% |IsocenterToRangeShifterDistance = getIsocenterToRangeShifterDistance(Beams , ScannerDirectory , BDL)|
%
%
%% Description
% |IsocenterToRangeShifterDistance = getIsocenterToRangeShifterDistance(Beams , ScannerDirectory , BDL)| Compute the range shifter to isocentre distance by placing the RS on the aperture block
%
%
%% Input arguments
% |Beams| - _STRUCT_ - Information about beam structure
%     * |Beams.RSinfo.RangeShifterMaterial| - _STRING_ - Name of the range shifter material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%     * |Beams.RSinfo.RangeShifterWET| -_SCALAR_- Water equivalent thickness (mm) of the range shifter to add to the CEF
%     * |Beams.Iso2Skin| -_SCALAR_- Isocentre (mm) To Skin Distance along the proton beam axis
%     * |Beams.Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%
% |ScannerDirectory| - _STRING_ - Name of the folder containing the definition of the CT scanner properties in MCsquare in folder "plugins\openMCsquare\lib\Scanners"
%
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
%% Output arguments
%
% |IsocenterToRangeShifterDistance| -_SCALAR_- Isocenter to DOWNSTEAM edge of range shifter (mm).
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function IsocenterToRangeShifterDistance = getIsocenterToRangeShifterDistance(Beams , ScannerDirectory , BDL)

  %The range shifter is upstream to the aperture block

  if ~isfield(Beams, 'ApertureBlock')
    ApertureBlock = 0;
  else
    ApertureBlock = Beams.ApertureBlock;
  end

  param = getMachineParam(BDL);
  dwstRS2Aper = param.snout.RangeShifterOffset(Beams.RSinfo.NbSlabs) - Beams.RSinfo.RSslabThickness(end); %distance (mm) from upstream aperture surface to downstream RS surface
  %distance        upstream surface or RS as defined in snout design         slab thickness
  %downstream surface
  %of RS to Upstream
  %surface of aperture

  if ApertureBlock
    %There is an aperture. The range shifter touches the aperture
    IsocenterToBlockTrayDistance = getIsocenterToBlockTrayDistance(Beams);
    IsocenterToRangeShifterDistance = IsocenterToBlockTrayDistance + dwstRS2Aper;
    % From DICOM standard : Isocenter to Range Shifter Distance (300A,0364)	Isocenter to DOWNSTREAM edge of range shifter (mm) at current control point.
    % |IsocenterToRangeShifterDistance| is therefore the distance to the patient side of the range shifter
    % The range shifter is upstream to the aperture block and touching the aperture block
  else
    %There is no aperture. The range shifter touches the skin
    IsocenterToRangeShifterDistance = Beams.Iso2Skin + dwstRS2Aper;
  end

end
