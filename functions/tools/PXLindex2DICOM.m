%% DICOM2PXLindex
% Convert from physical coordiante (mm) in DICOM CS in to the (x,y,z) pixel index in the 3D image
% If |isocentre| is present, then the origin of the DICOM CS is placed at isocentre
%
%% Syntax
% |P = DICOM2PXLindex(Axyz , Spacing , ImagePositionPatient )|
%
% |P = DICOM2PXLindex(Axyz , Spacing , ImagePositionPatient , isocenter)|
%
%
%% Description
% |Axyz = DICOM2PXLindex(P , Spacing , ImagePositionPatient)| Compute indices when |P| is expressed in a coordinate system with origin at origin of DICOM CS
%
% |Axyz = DICOM2PXLindex(P , Spacing , ImagePositionPatient , isocenter)| Compute indices when |P| is expressed in a coordinate system with origin at isocentre
%
%
%% Input arguments
% |Axyz| -_SCALAR MATRIX_- |Axyz(i,:) = [xi , yi , zi]| The pixel indices along each axis of the CT scan for the i-th pixel
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size (|mm|) of the displayed images in GUI
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
% |isocenter| -_SCALAR VECTOR_- [OPTIONAL] [x,y,z] Coordiantes (mm) of the isocentre in the DICOM CS
%
%
%% Output arguments
%
% |P| -_SCALAR MATRIX_- |P(i,:) = [x,y,z,1]| The coordinate (mm) of the i-th voxel in DICOM IEC
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function P = PXLindex2DICOM(Axyz , Spacing , ImagePositionPatient , isocenter)

  if nargin < 4
    isocenter = [0,0,0];
  end

  %NB Axyz-1 : the indexing in Matlab starts at 1. The pixel at origin is at index (1,1,1) and not (0,0,0)
  P(:,1) = (Axyz(:,1)-1) .* Spacing(1) + ImagePositionPatient(1) - isocenter(1); %Coordiante of the body voxels in DICOM CS
  P(:,2) = (Axyz(:,2)-1) .* Spacing(2) + ImagePositionPatient(2) - isocenter(2);
  P(:,3) = (Axyz(:,3)-1) .* Spacing(3) + ImagePositionPatient(3) - isocenter(3);

  P(:,4) = ones(size(P,1),1);
end
