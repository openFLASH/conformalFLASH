%% DICOM2PXLindex
% Convert from physical coordinate (mm) in DICOM CS in to the (x,y,z) pixel index in the 3D image
%
%% Syntax
% |Axyz = DICOM2PXLindex(P , Spacing , ImagePositionPatient)|
%
% |Axyz = DICOM2PXLindex(P , Spacing , ImagePositionPatient , isInteger)|
%
% |[~ , X, Y , Z]= DICOM2PXLindex([] , Spacing , ImagePositionPatient , isInteger, X , Y , Z)|
%
%% Description
% |Axyz = DICOM2PXLindex(P , Spacing , ImagePositionPatient)| Compute indices when |P| is expressed in a coordinate system with origin at origin of DICOM CS
%
% |[~ , X, Y , Z]= DICOM2PXLindex([] , Spacing , ImagePositionPatient , isInteger, X , Y , Z)| Compute the indices when coodinates are expressed as separate X,Y,Z vecotr. This call is  4 * faster.
%
%% Input arguments
% |P| -_SCALAR MATRIX_- |A(i,:) = [x,y,z]| The coordinate of the i-th voxel in DICOM IEC
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size (|mm|) of the displayed images in GUI
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
% |isInteger| -_BOOLEAN_- [OTPIONAL: default: true] |isInteger=true| indicates that the pixel index is rounded to the closest integer
%
% |X| -_SCALAR VECTOR_- |X(i)| [OPTIONAL. required only if |P| is empty]. The X coordinate (mm) of the i-th voxel in DICOM IEC
%
% |Y| -_SCALAR VECTOR_- |Y(i)| [OPTIONAL. required only if |P| is empty]. The Y coordinate (mm) of the i-th voxel in DICOM IEC
%
% |Z| -_SCALAR VECTOR_- |Z(i)| [OPTIONAL. required only if |P| is empty]. The Z coordinate (mm) of the i-th voxel in DICOM IEC
%
%% Output arguments
%
% |Axyz| -_SCALAR MATRIX_- |Axyz(i,:) = [xi , yi , zi]| The pixel indices along each axis of the CT scan for the i-th pixel.
%
% |X| -_SCALAR VECTOR_- |X(i)|  The X coordinate (pixel index) of the i-th voxel in DICOM IEC
%
% |Y| -_SCALAR VECTOR_- |Y(i)| The Y coordinate (pixel index) of the i-th voxel in DICOM IEC
%
% |Z| -_SCALAR VECTOR_- |Z(i)| The Z coordinate (pixel index) of the i-th voxel in DICOM IEC

%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Axyz , X, Y , Z ]= DICOM2PXLindex(P , Spacing , ImagePositionPatient , isInteger, X , Y , Z )

  if nargin < 4
    isInteger = true;
  end

  if isInteger
    %Round the pixel index to the closest integer
    if ~isempty(P)
        %NB: The computation of Axyz is 4* slower than the computation of [X,Y,Z]
        %NB +1 : the indexing in Matlab starts at 1. The pixel at origin is at index (1,1,1) and not (0,0,0)
        Axyz(:,1) = 1 + round((P(:,1) - ImagePositionPatient(1)) ./ Spacing(1));
        Axyz(:,2) = 1 + round((P(:,2) - ImagePositionPatient(2)) ./ Spacing(2));
        Axyz(:,3) = 1 + round((P(:,3) - ImagePositionPatient(3)) ./ Spacing(3));
        X = [];
        Y = [];
        Z = [];
    else
        %NB: This computation is much faster
        X = 1 + round((X - ImagePositionPatient(1)) ./ Spacing(1));
        Y = 1 + round((Y - ImagePositionPatient(2)) ./ Spacing(2));
        Z = 1 + round((Z - ImagePositionPatient(3)) ./ Spacing(3));
        Axyz = [];
    end
  else
    %Do not round the pixel index to an integer value
    %This is usefull if the index is used in an interpolation.
      if ~isempty(P)
          Axyz(:,1) = 1 + (P(:,1) - ImagePositionPatient(1)) ./ Spacing(1);
          Axyz(:,2) = 1 + (P(:,2) - ImagePositionPatient(2)) ./ Spacing(2);
          Axyz(:,3) = 1 + (P(:,3) - ImagePositionPatient(3)) ./ Spacing(3);
          X = [];
          Y = [];
          Z = [];
      else
        X = 1 + (X - ImagePositionPatient(1)) ./ Spacing(1);
        Y = 1 + (Y - ImagePositionPatient(2)) ./ Spacing(2);
        Z = 1 + (Z - ImagePositionPatient(3)) ./ Spacing(3);
        Axyz = [];
      end
  end

end
