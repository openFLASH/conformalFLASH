%% function
% Convert coordinates from pixel index to map
%
%% Syntax
% |pt = pxl2map(pt , FistCorner , gap)|
%
%
%% Description
% |pt = pxl2map(pxlIdx , FistCorner , gap)| Description
%
%
%% Input arguments
% |pxlIdx| -_SCALAR MATRIX_- |pxlIdx(i,:)= [j,k]| Pixel coordinate (j,k) of the i-th point
%
% |FistCorner| -_SCALAR VECTOR_- Coordinate of top left corner in the physics coordinate system
%
% |gap| -_SCALAR VECTOR_- Pixel size of |mapImage|  in the physics coordinate system
%
%
%% Output arguments
%
% |pt| - _SCLAR MATRIX_ - |pt(i,:) = [x,y]| The physics coordinates of the i-th pixels
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)


function pt = pxl2map(pxlIdx , FistCorner , gap)
  pt = (pxlIdx-1) .* gap + repmat(FistCorner,size(pxlIdx,1),1);
end
