%% interpolateCT
% Interpolate the CT on the new grid from the CT on the old grid
% Make sure the the indices of the new CT are oriented along the same axes as the old CT
% The coordinate of (Xint, Yint, Zint) of the voxels of the interpolated CT scan are expressed in the smae coordinate than the original CT scan
%
%% Syntax
% |CTintrp = interpolateCT(CT, CTx, CTy, CTz , Xint, Yint, Zint , HUextpl , CTintSize)|
%
%
%% Description
% |CTintrp = interpolateCT(CT, CTx, CTy, CTz , Xint, Yint, Zint , HUextpl , CTintSize)| Description
%
%
%% Input arguments
% |CT| - _SCALAR MATRIX_ - |CT(i,j,k)=HU| gives the intensity (Hounsfield unit) of the voxel at location (CTx(i),CTy(j),CTz(k))
%
% |CTx| -_SCALAR VECTOR_- |CTx(i)| Coordinate (mm) of the i-th voxel along the axis of the first pixel index
%
% |CTy| -_SCALAR VECTOR_-  |CTy(j)| Coordinate (mm) of the j-th voxel along the axis of the second pixel index
%
% |CTz| -_SCALAR VECTOR_-  |CTz(k)| Coordinate (mm) of the k-th voxel along the axis of the third pixel index
%
% |Xint| -_SCALAR VECTOR_- |Xint(u)| Coordinate (mm) of the u-th voxel of the interpolated CT scan
%
% |Yint| -_SCALAR VECTOR_- |Yint(u)| Coordinate (mm) of the u-th voxel of the interpolated CT scan
%
% |Zint| -_SCALAR VECTOR_- |Zint(u)| Coordinate (mm) of the u-th voxel of the interpolated CT scan
%
% |HUextpl| -_SCALAR_- HU value to place in the voxels that are extrapolated (outside the box of the original CT scan)
%
% |CTintSize| -_SCALAR VECTOR_- Number of pixels along each axis of the interpolated CT scan. prod(CTintSize) = numel(Xint)*numel(Yint)*numel(Zint)
%
%
%% Output arguments
%
% |CTintrp| - _SCALAR MATRIX_ - Interpolated CT |CTintrp|. This is a 3D matrix with size |CTintSize|. The orering of the voxels is determined by the ordering of the voxels in Xint, Yint, Zint
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function CTintrp = interpolateCT(CT, CTx, CTy, CTz , Xint, Yint, Zint , HUextpl , CTintSize , interMethod)

  if nargin < 10
    interMethod = 'nearest';
  end

  [Yct,Xct,Zct] = meshgrid(CTy ,CTx, CTz ); %Meshgrid inverse X and Y. And it is inversed again in interp3 for Xct,Yct
  CTintrp = interp3(double(Yct),double(Xct),double(Zct), CT, double(Yint), double(Xint), double(Zint),interMethod,HUextpl); %The Yint, Xint in interp3
  CTintrp = reshape(CTintrp, CTintSize(1) , CTintSize(2), CTintSize(3));


end
