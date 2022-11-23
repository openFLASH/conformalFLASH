%% getCTaxes
% Compute the physical coordinate of the 3 axes of the CT scan
% from the information in handles
%
%% Syntax
% |[Xvec, Yvec , Zvec , X , Y , Z] = getCTaxes(ImagePositionPatient , Spacing , sizeCT , isocenter)|
%
%
%% Description
% |[Xvec, Yvec , Zvec , X , Y , Z] = getCTaxes(ImagePositionPatient , Spacing , sizeCT , isocenter)| Description
%
%
%% Input arguments
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size [x,y,z] (|mm|) of the voxels of the CT scan
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size [x,y,z] (|mm|) of the voxels of the CT scan
%
% |sizeCT| - _SCALAR VECTOR_ - [x,y,z] Number of voxel in the 3 axes of the CT scan
%
% |isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
%
%% Output arguments
%
% |Xvec| -_SCALAR VECTOR_- |Xvec(i)| Coordinate (mm) of the i-th voxel along the axis of the first pixel index
%
% |Yvec| -_SCALAR VECTOR_-  |Yvec(j)| Coordinate (mm) of the j-th voxel along the axis of the second pixel index
%
% |Zvec| -_SCALAR VECTOR_-  |Zvec(k)| Coordinate (mm) of the k-th voxel along the axis of the third pixel index
%
% |X| -_SCALAR MATRIX_-  |X(i,j,k)| X coordinate (mm) of the voxel at pixel position (i,j,k) in the CT scan
%
% |Y| -_SCALAR MATRIX_-  |Y(i,j,k)| X coordinate (mm) of the voxel at pixel position (i,j,k) in the CT scan
%
% |Z| -_SCALAR MATRIX_-  |Z(i,j,k)| X coordinate (mm) of the voxel at pixel position (i,j,k) in the CT scan
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Xvec, Yvec , Zvec , X , Y , Z] = getCTaxes(ImagePositionPatient , Spacing , sizeCT , isocenter)

  Xvec = ImagePositionPatient(1) - isocenter(1) + [0:sizeCT(1)-1]' * Spacing(1);  %Origin of the original CT is moved to isocentre
  Yvec = ImagePositionPatient(2) - isocenter(2) + [0:sizeCT(2)-1]' * Spacing(2);
  Zvec = ImagePositionPatient(3) - isocenter(3) + [0:sizeCT(3)-1]' * Spacing(3);

  [Y , X , Z] = meshgrid( Yvec , Xvec , Zvec);

end
