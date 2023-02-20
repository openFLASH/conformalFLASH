%% updateROI
% Parse the RT structure and create the sparse matrice that will be used by ComputePijsMatrix
%
%% Syntax
% |res = help_header(im1,im2)|
%
%
%% Description
% |res = help_header(im1,im2)| Description
%
%
%% Input arguments
% |handles| - _STRUCT_ - REGGUI data structure.
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |ROI| - _struct_ - MIROpt structure containing information about all
%           volumes in the RTSTRUCT file. The following data must be present in the structure:
%     * |ROI(i).mask3D.value|- _array_ - |ROI(i).mask3D.value(x,y,z)=1| if the voxel at (x,y,z) is located inside the RT struct
%
%
%% Output arguments
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
% |ROI| - _struct_ - MIROpt structure containing information about all
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [handles , Plan , ROI] = updateROI(handles , Plan , ROI)

  nvoxels = prod(handles.size);
  Plan.OptROIVoxels_nominal = false(nvoxels,1); % initialize to zeros
  Plan.OptROIVoxels_robust = false(nvoxels,1); % initialize to zeros

  %Create the ROI structure
  for i = 1:numel(ROI)

      % Calculate mask for each ROI
      ROI(i).mask3D.value = Get_reggui_data(handles, ROI(i).name,'images');

      % Convert to 1D sparse mask with z inverse format (to be consistent
      % with beamlets format)
      temp = flip(ROI(i).mask3D.value,3);
      ROI(i).mask1D = sparse(logical(double(temp(:))));
      ROI(i).nvoxels = sum(ROI(i).mask1D);
      ROI(i).voxelsID = find(ROI(i).mask1D);

      if (i < numel(ROI)) %Do not add EXTERNAL to the mask
        Plan.OptROIVoxels_nominal = Plan.OptROIVoxels_nominal | ROI(i).mask1D; %Mask defining the voxels of the dose influence matrix to be loaded in the sparse matrix
      end
  end

end
