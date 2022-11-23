%% getIsocentreToSkinDistance
% Compute the distance between the isocentre and the patient skin for the specified beam
% Take into account the radius of the aperture and avoid collision with |Body| contour within
% |Rmax| radius around the beam axis
%
%% Syntax
% |Zskin = getIsocentreToSkinDistance(Beam , Body , Rmax , Spacing , ImagePositionPatient)|
%
%
%% Description
% |Zskin = getIsocentreToSkinDistance(Beam , Body , Rmax , Spacing , ImagePositionPatient)| Description
%
%
%% Input arguments
%
% |Beams| -_STRUCTURES_- Information about the beam
%     * |Beams.GantryAngle| -_SCALAR_- Gantry Angle (deg)
%     * |Beams.PatientSupportAngle| -_SCALAR_- Couch Angle (deg)
%     * |Beams.isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
%     * |Beams.TargetMargin| -_SCALAR VECTOR_- [OPTIONAL]  [x,y,z] Margin around the target to take into account penumbra and setup errors (mm). Used for the definition of the spot positions
%
% |Body| - _SCALAR MATRIX_ - Mask defining the contour of the body. This is used to define the position of hte aperture |PTV(x,y,z)=1| if the voxel is inside the PTV
%
% |Rmax| -_SCALAR_- Maximum radisu of the aperture block to consider to estimate collision with |Body|  contour
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size (|mm|) of the displayed images in GUI
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
%
%% Output arguments
%
% |Zskin| - _SCALAR_ - Distance (mm) from isocentre to skin surface along the axis of the specified beam
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Zskin = getIsocentreToSkinDistance(Beam , Body , Rmax , Spacing , ImagePositionPatient)

  Bdcm = getDICOMcoord(Body, Spacing , ImagePositionPatient , Beam.isocenter); %Coordinate of body voxel in DICOM CS (with origin at isocentre)
  M = matDICOM2IECgantry(Beam.GantryAngle,Beam.PatientSupportAngle);
  Bbev = M * Bdcm; %Coordinate of the voxels of the body in IC gantry
  [~, nearAxis] = find(sum(Bbev(1:2,:).^2 , 1) <= Rmax.^2); %Find the voxels closest to optical axis of beamlet. Their distance to beam aixs is less than SpotSpacing
  Zskin = max(Bbev(3,nearAxis)); %Position on beam axis of the skin + some margin in mm
  %The point at the skin surface has the highest Zg (it is the further away from isocentre and towards the source) and is near the optical axis

end
