%% elevationMap2mask3D
% Convert a 2D elevation map into a 3D binary mask
% The (x,y) size of the pixels in the mask is the size of the pixels in |ElvMap|
%
%% Syntax
% |Mask = elevationMap2mask3D(ElvMap , Z_cem)|
%
%
%% Description
% |Mask = elevationMap2mask3D(ElvMap , Z_cem)| Description
%
%
%% Input arguments
% |ElvMap| - _SCALAR MATRIX_ - |ElvMap(i,j)| Height (mm) of the column at position (i,j). The height is measured from |Z_cem=0|
%
% |Z_cem|  -_SCALAR VECTOR_- Z coordinate (mm) in CEM obect. |Z_cem=0| at the base of CEM. |Z_cem| increases when moving up in the CEM
%                      |Z_cem| must be a continuously increasing vector
%
%% Output arguments
%
% |Mask| - _SCALAR MATRIX_ - Binary mask |Mask(i,j,k) = 1| if the voxel is contained inside the CEM defined by the elevation map
%               |Mask(i,j,k)| is inside the CEM if |ElvMap(i,j) > Z_cem(k)|
%               The size of |Mask| is [size(ElvMap,1) , size(ElvMap,2) , numel(Z_cem)]
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Mask = elevationMap2mask3D(ElvMap , Z_cem)

  %Check that the Z vector is properly sorted
  if min(diff(Z_cem)) <0
    error('Z_cem must be a continuously increasing vector')
  end

  GridSize = size(ElvMap);
  NbLayers = numel(Z_cem);
  Mask = zeros(GridSize(1) , GridSize(2) , NbLayers); %Make the 3D mask

  %Add the CEM in the mask
  baseIDX = min(find(Z_cem >= 0)); %The base of the spikes are at Z_cem=0

  for Zgidx = baseIDX:NbLayers
    binImg = ElvMap > Z_cem(Zgidx); %Start at 0, so must be strictly > to have the correct delta
    idx = find(binImg) ; %Indices of the pixels within this layer
    [i,j] = ind2sub(GridSize , idx); %2d indices
    Aidx = sub2ind(size(Mask) , i , j , repmat(Zgidx,numel(i),1));
    Mask(Aidx) = 1; %Fill the mask of the spike in the 'CT' of the CEM
  end

end
