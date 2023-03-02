%% getIsocenterToRangeShifterDistance
% Compute the distance from isocentre to the DOWNSTREAM edge of range shifter
% The range shifter can be placed either:
% * on the skin surface if no aperture is present
% * on the aperture block if an aperture block is provided
%
% From DICOM standard : Isocenter to Range Shifter Distance (300A,0364)	Isocenter to DOWNSTREAM edge of range shifter (mm) at current control point.
% |IsocenterToRangeShifterDistance| is therefore the distance to the patient side of the range shifter

%
%% Syntax
% |IsocenterToRangeShifterDistance = getIsocenterToRangeShifterDistance( Beams )|
%
%
%% Description
% |IsocenterToRangeShifterDistance = getIsocenterToRangeShifterDistance( Beams )| Compute the range shifter to isocentre distance by placing the RS on the aperture block
%
%
%% Input arguments
% |Beams| - _STRUCT_ - Information about beam structure
%     * |Beams.RSinfo.SlabOffset| -_SCALAR or SCALAR VECTOR_- Distance (mm) from upstream aperture surface (= snout position) to the upstream surface of the i-th slab of range shifter
%     * |Beams.RSinfo.RSslabThickness| -_SCALAR or SCALAR VECTOR_- |RSslabThickness(s)| Thickness (mm) of the s-th slab of the range shifter
%     * |Beams.Iso2Skin| -_SCALAR_- Isocentre (mm) To Skin Distance along the proton beam axis
%
% |Snout2RSOffset| -_SCALAR_- distance (mm) between the plane defining the snout position and the donwstream surface of the range shifter
%
%% Output arguments
%
% |IsocenterToRangeShifterDistance| -_SCALAR_- Isocenter to DOWNSTEAM edge of range shifter (mm).
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function IsocenterToRangeShifterDistance = getIsocenterToRangeShifterDistance( Beams , Snout2RSOffset)

  %The range shifter is upstream to the aperture block

  if ~isfield(Beams, 'ApertureBlock')
    ApertureBlock = 0;
  else
    ApertureBlock = Beams.ApertureBlock;
  end

  if ApertureBlock
    %There is an aperture. The range shifter touches the aperture
    IsocenterToBlockTrayDistance = getIsocenterToBlockTrayDistance(Beams);
    IsocenterToRangeShifterDistance = IsocenterToBlockTrayDistance + Snout2RSOffset;
    % From DICOM standard : Isocenter to Range Shifter Distance (300A,0364)	Isocenter to DOWNSTREAM edge of range shifter (mm) at current control point.
    % |IsocenterToRangeShifterDistance| is therefore the distance to the patient side of the range shifter
    % The range shifter is upstream to the aperture block and touching the aperture block
  else
    %There is no aperture. The range shifter touches the skin
    IsocenterToRangeShifterDistance = Beams.Iso2Skin + Snout2RSOffset;
  end

  IsocenterToRangeShifterDistance = round(IsocenterToRangeShifterDistance,1);

end
