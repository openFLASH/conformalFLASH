%% loadEmptyStructs
% Load the RT strcut that was used to create the RT plan on the Raystation
% Update the REGGUI handles and the PLan structure and crate the ROI structure
%
%% Syntax
% |[handles , Plan ] = loadEmptyStructs( rtstructFileName, handles, Plan, ExternalROI ,percentile)|
%
%
%% Description
% |[handles , Plan ] = loadEmptyStructs( rtstructFileName, handles, Plan, ExternalROI ,percentile)| Description
%
%
%% Input arguments
%
% |rtstructFileName| -_STRING_- File name and full path to the DICOM RT struct. If not empty, load the file from disk.
%                               If empty, the Rt strcut must be in handles.
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
% |ExternalROI| -_STRING_- Name of the RT struct with the body contour
%
% |percentile| -_SCLAR_- [OPTIONAL. If absent use default of DRaEstimate] The percentile to compute the percentile dose rate
%
%% Output arguments
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
% |ROI| - _struct_ - MIROpt structure containing information about all
%           volumes in the RTSTRUCT file. The following data must be present in the structure:
%     * |ROI(i).mask3D.value|- _array_ - |ROI(i).mask3D.value(x,y,z)=1| if the voxel at (x,y,z) is located inside the RT struct
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [handles , Plan ] = loadEmptyStructs( rtstructFileName, handles, Plan, ExternalROI ,percentile)

  %Load the RT struct
  if ~isempty(rtstructFileName)
    %We received a file name. Load the RT strcut from file
    [handles , structName] = Import_contour(rtstructFileName,{ExternalROI},Plan.CTname,1,handles);
    Plan.ExternalROI = structName{1};
  else
    %We did not receive a file name. The structure is already loaded in handles. Just copy the name for future references
    Plan.ExternalROI = ExternalROI;
  end

  % Create ROI specific properties for all structures
  Plan.DoseGrid.size = handles.size';
  Plan.OptROIVoxels_nominal = createROI(Plan.ExternalROI, handles);
  Plan.OptROIVoxels_robust = []; % initialize to zeros

  %Create the ROI structure
  optFidx = 1;
  Plan.optFunction(optFidx) = createOptFunctionStruct();
  Plan.optFunction(optFidx).ROIname = Plan.ExternalROI;

  if nargin >= 5
    Plan.optFunction(optFidx).Vref = percentile;
  end

end

%------------------------------------
% Create a sttructure defining the dose rate objective function
%------------------------------------
function optFunction = createOptFunctionStruct()

  optFunction.ID = 8; %Select this ID because the optimisation function ID=8 requires the computation of dose rate
  optFunction.Dref  = 0; %Set the threshold to 0, so that all voxels are included in the computation of the dose rate

end

%---------------------------
% add the ROI to the list of ROIs
%-----------------------------
function mask1D = createROI(ROIname, handles)

    ROImask = Get_reggui_data(handles, ROIname);

    % Convert to 1D sparse mask with z inverse format (to be consistent
    % with beamlets format)
    temp = flip(ROImask,3);
    mask1D = sparse(logical(double(temp(:))));

  end
