%% getIsocenterToBlockTrayDistance
% Compute the distance from isocentre to the base of the aperture
% Compute the thickness of the aperture required to stop the proton. The aperture block will stop the proton of the maximum energy layer.
% If there is a range shifter, it is assumed that the range shifter is upstream of the aperture and therefore degrades the proton energy
% Place the aperture upstream to the skin position
%
% The dicom tag (300A,00FB) BlockMountingPosition defines the reference surface for the measurement of |IsocenterToBlockTrayDistance|

%
%% Syntax
% |[IsocenterToBlockTrayDistance , BlockThickness , Zskin , BlockMountingPosition] = getIsocenterToBlockTrayDistance(Beam)|
%
%
%% Description
% |[IsocenterToBlockTrayDistance , BlockThickness , Zskin , BlockMountingPosition] = getIsocenterToBlockTrayDistance(Beam)| Description
%
%
%% Input arguments
% |Beam| -_STRUCTURES_- Information about the beam
%     * |Beam.Iso2Skin| -_SCALAR_- Isocentre (mm) To Skin Distance along the proton beam axis
%     * |Beam.Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%
%
%% Output arguments
%
% |IsocenterToBlockTrayDistance| - _SCALAR_ - Distance (mm) from isocentre to the UPSTREAM base of the aperture
%
% |BlockThickness| - _SCALAR_ - Thickness (mm) of the aperture block
%
% |Zskin| - _SCALAR_ - Z coordinate (mm) in IEC gantry of the skin surface
%
% |BlockMountingPosition| -_STRING_- The dicom tag (300A,00FB) BlockMountingPosition
%                           SOURCE_SIDE is using the downstream face of the block as a reference position for expressing the isocenter to block tray distance.
%                           PATIENT_SIDE is using the upstream face
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [IsocenterToBlockTrayDistance , BlockThickness , Zskin , BlockMountingPosition] = getIsocenterToBlockTrayDistance(Beam)

  %There is an aperture block for this beam. Create it now
  Energy = max([Beam.Layers(:).Energy]); %Incoming beam energy

  %If there is a range shifter, degrade energy in the range shifter
  if (Beam.NumberOfRangeShifters)
    Energy = energyAfterRangeShifter(Beam.RSinfo.RangeShifterMaterial , Energy , sum(Beam.RSinfo.RSslabThickness));
  end

  %Compute the thickness of the brass block
  % Find the depth at which place the aperture
  % It will touch the skin
  Zskin = Beam.Iso2Skin; % Isocentre (mm) To Skin Distance along the proton beam axis
  BlockThickness = getApertureThickness(Energy);

  if isfield(Beam , 'AirGap')
      AirGap = Beam.AirGap;
  else
      AirGap = 0;
  end
  IsocenterToBlockTrayDistance = round(Zskin + BlockThickness + AirGap , 1); %Distance from isocentre to upstream side of of aperture block;
  BlockMountingPosition = 'PATIENT_SIDE';
  %SOURCE_SIDE and PATIENT_SIDE are used to indicate which face of the accessory is used for defining its position
  %SOURCE_SIDE is using the downstream face of the block as a reference position for expressing the isocenter to block tray distance.
  %PATIENT_SIDE is using the upstream face
  %We use the upstream side of the block to define the distance so that it is independent of the block thickness
  % By defining the snout position at the same reference point as the iso to block tray distance,
  % both distances are equal and independent of the block thickness


end
