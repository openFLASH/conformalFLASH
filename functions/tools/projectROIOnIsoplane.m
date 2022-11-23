%% projectPTVOnIsoplane
% Project a mask from a matrix aligned with the axis of the DICOM CS onto a plane aligned with the IEC gantry CS
%
%% Syntax
% |[bevROI , X, Y , pxlSize] = projectROIOnIsoplane(Beams , ROI , Spacing , ImagePositionPatient)|
%
%
%% Description
% |[bevROI , X, Y , pxlSize] = projectROIOnIsoplane(Beams , ROI , Spacing , ImagePositionPatient)| Description
%
%
%% Input arguments
% * |Beams| -_STRUCTURES_- Information about the beam
%     * |Beams.GantryAngle| -_SCALAR_- Gantry Angle (deg)
%     * |Beams.PatientSupportAngle| -_SCALAR_- Couch Angle (deg)
%     * |Beams.isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
%
% |ROI| - _SCALAR MATRIX_ - Mask defining the position of the ROI |ROI(x,y,z)=1| if the voxel is inside the region of interest
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
% |pxlSize| -_SCALAR VECTOR_- [OPTIONAL: default pxlSize = min(Spacing) ./2] Size (mm) of the square pixels of the |bevROI|
%
%% Output arguments
%
% |bevROI| - SCALAR MATRIX_ -  Binary mask |bevROI(x,y)=1| if the pixel is in the projection of the PTV in the isocentre plane
%
% |X| -_SCALAR VECTOR_- |X(i,j)| X Coordinate in IEC gantry of the pixel at |bevROI(x,y)|
%
% |Y| -_SCALAR VECTOR_- |X(i,j)| Y Coordinate in IEC gantry of the pixel at |bevROI(x,y)|
%
% |pxlSize| -_SCALAR VECTOR_- [x,y] Size (mm) of the pixels of the |bevROI|
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [bevROI , X, Y , pxlSize] = projectROIOnIsoplane(Beams , ROI , Spacing , ImagePositionPatient , pxlSize)

  if nargin < 5
    %Set the default pixel size of the output
    pxlSize = min(Spacing) ./2;
  end

  %Create a description of the projection axis
  beam.gantry_angle = Beams.GantryAngle; %- _SCALAR_ - Gantry angle (in degree) of the treatment beam beam
  if (isfield(Beams,'PatientSupportAngle'))
    beam.table_angle = Beams.PatientSupportAngle;
  else
    beam.table_angle = 0; % - _SCALAR_ - Table top yaw angle (degree) of the treatment beam
  end
  beam.isocenter = Beams.isocenter; %_SCALAR VECTOR_ - |beam.isocenter= [x,y,z]| Cooridnate (in mm) of the isocentre in the CT scan for the treatment beam


  %Define the dimension of the BEV window
  Ind = find(ROI);
  Ind = unique(Ind);
  [x,y,z] = ind2sub(size(ROI),Ind); %Find the voxels inside the body to define the BEV map

  bev_size = getBEVsize(x,y,z , beam);
  bev_size = round(2.5 .* max(bev_size)); %BEV window is equal to the longest diameter of the RT struc + a margin


  %Create a projection grid in the plane of the isocentre
  SAD = 1000; %mm Arbitrary as the rays are paralell. Place the plane of the source outside of CT scan
  handles.spacing = Spacing; %| - _SCALAR VECTOR_ - Pixel size (|mm|) of the displayed images in GUI
  handles.origin = ImagePositionPatient; %| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image

  [bevROI , bev_v , pts_iso] = Project_on_isoplane(handles,ROI,beam,SAD,bev_size , pxlSize);
  bevROI = bevROI'; %The first index is the Ygantry axis. We need to transpose
  bevROI = bevROI > 0;

  res = bev_v(2) - bev_v(1); %Spatial resolution (mm). All points are equidistant
  pxlSize = [res,res];
  [Y,X] = meshgrid(bev_v,bev_v); %First index is Y in meshgrid


end
