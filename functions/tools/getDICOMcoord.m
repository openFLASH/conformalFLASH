%% getDICOMcoord
% Get a 4 vector [x,y,z,1] with the coordinates of the voxels contained in a mask
% Coordinates expressed in the DICOM CS with the origin moved to isocentre.
%
%% Syntax
% |Bdcm = getDICOMcoord(ROI, Spacing , ImagePositionPatient , isocenter)|
%
%
%% Description
% |Bdcm = getDICOMcoord(ROI, Spacing , ImagePositionPatient , isocenter)| Description
%
%
%% Input arguments
% |ROI| -_SCALAR MATRIX_- |A(x,y,z) = 1| if the voxel belong to the ROI
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size (|mm|) of the displayed images in GUI
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
% |isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the DICOM CS
%
%
%% Output arguments
%
% |Bdcm| 4 vector [x,y,z,1] with the coordiantes of the mask voxels expressed in DICOM CS with oriting at isocentre
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Bdcm = getDICOMcoord(ROI, Spacing , ImagePositionPatient , isocenter)

  Ind = find(ROI);

  [x,y,z] = ind2sub(size(ROI),Ind);
  Ind = [];
  Bdcm = PXLindex2DICOM([x,y,z] , Spacing , ImagePositionPatient , isocenter);
  Bdcm = Bdcm';

end
