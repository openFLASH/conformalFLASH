%% updateAllImages
% Update the CT scan of |handles.images| with name |CTname|
% Use |Plan.CTinfo| to add information to the image in |handles|
% Resize all images in handles.images to fit the enlarged CT
% The old smaller images are inserted in the larger new image at position |oldOrigin| of the DICOM CS
%
%% Syntax
% |[handles , Plan] = updateAllImages(handles , Plan , CT , NewOrigin ,  HUvoid  , CTname)|
%
%
%% Description
% |[handles , Plan] = updateAllImages(handles , Plan , CT , NewOrigin ,  HUvoid  , CTname)| Description
%
%
%% Input arguments
% |handles| -_STRUCTURE_- REggui data handle
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are
%
% |CT| - _SCALAR MATRIX_ - CT scan to be added |CT(x,y,z)=HU| gives the intensity (Hounsfield unit) of the voxel at location (x,y,z)
%
% |NewOrigin|  -_SCALAR VECTOR_- |NewOrigin= [x,y,z]| Coordinate (mm) of the first pixel of the new larger image
%
% |HUvoid| -_SCALAR_- Hounsfield unit to use for unoccupied voxels in larger CT scan. If the image is a binary mask |HUvoid| is ignored and the new space is filled with 0
%
% |CTname| -_STRING_- Name of the CT scan in |handles.images| where to save the results
%
%% Output arguments
%
% |handles| -_STRUCTURE_- Updated REggui data handle with larger images
%
% |Plan| - _struct_ - MIROpt structure with updated plan parameters
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [handles , Plan] = updateAllImages(handles , Plan , CT , NewOrigin ,  HUvoid  , CTname)

  NewSize = size(CT);

  %Resize all images in handles to fit new image format

  if (sum(handles.origin - NewOrigin) | sum(handles.origin - NewSize))
    %If the CT size was updated, resize all images in handles
    handles = resizeAllHandles(handles , NewOrigin , NewSize , HUvoid);
  end %if (sum(CTsize - oldSize))

  %Update the CT in the handles
  Plan  = updatePlanCTparam(handles , Plan);
  handles = Set_reggui_data(handles , CTname , CT , Plan.CTinfo , 'images',1); %Update the CT scan with the aperture block in handles

end

%===============================================
% Resize all images contained in handles
% If the image is a grey level image, the new voxels are filled with |HUvoid|
% The the image is a binary mask, the new voxel are filled with 0
%
%% Input arguments
% |handles| -_STRUCTURE_- REggui data handle
% |NewOrigin|  -_SCALAR VECTOR_- |NewOrigin= [x,y,z]| Coordinate (mm) of the first pixel of the new larger image
% |NewSize| -_SCALAR VECTOR_- [x,y,z] number of voxels of the enlarged CT scan
% |HUvoid| -_SCALAR_- Hounsfield unit to use for unoccupied voxels in larger CT scan. If the image is a binary mask |HUvoid| is ignored and the new space is filled with 0
%
% OUTPUT
% |handles| -_STRUCTURE_- Updated REggui data handle with larger images
%     * |handles.images|
%     * |handles.origin|
%     * |handles.size|
%===============================================
function handles = resizeAllHandles(handles , NewOrigin , NewSize , HUvoid)

    %Start at 2 to avoid image 'none'
    for i = 2:numel(handles.images.name)

      if sum(size(handles.images.data{i}) - NewSize)
         values = unique(handles.images.data{i});
         if (min(values)==0 & max(values) == 1)
           %This is a binary mask. Fill with zeros
           handles.images.data{i} = enlargeCT(handles.images.data{i} , handles.origin , NewOrigin , handles.spacing , NewSize , 0);
         else
           %This is a grey level image. Fill with HUvoid
           handles.images.data{i} = enlargeCT(handles.images.data{i} , handles.origin , NewOrigin , handles.spacing , NewSize , HUvoid);
         end
      end

    end

    %Update the handles
    handles.origin = NewOrigin;
    handles.size = NewSize'; %handles.size is a column vector

end
