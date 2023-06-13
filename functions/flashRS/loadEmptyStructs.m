%% loadRSrtStructs
% Load the RT strcut that was used to create the RT plan on the Raystation
% Update the REGGUI handles and the PLan structure and crate the ROI structure
%
%% Syntax
% |[handles , Plan , ROI] = loadRSrtStructs(rtstructFileName , handles, Plan, RTstruct , DoseRate)|
%
%
%% Description
% |[handles , Plan , ROI] = loadRSrtStructs(rtstructFileName , handles, Plan, RTstruct , DoseRate)| Description
%
%
%% Input arguments
%
% |rtstructFileName| -_STRING_- File name and full path to the DICOM RT struct
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
% |RTstruct| -_STRUCTURE_- Information about the RT structs used in the TPS
%   * |RTstruct.selected_ROIs| -_CELL VECTOR_- |RTstruct.selected_ROIs{i}| is a stirng with the name of the i-th RT struct in which the dose rate is to be computed
%   * |RTstruct.ExternalROI| -_STRING_- Name of the RT struct with the body contour
%   * |RTstruct.TargetROI| -_STRING_- Name of the RT struct with the PTV
%
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

function [handles , Plan ] = loadEmptyStructs( rtstructFileName, handles, Plan, RTstruct )


  %Load the RT struct
  strName =  RTstruct.ExternalROI;
  handles = Import_contour(rtstructFileName,{strName},Plan.CTname,1,handles);

  % Create ROI specific properties for all structures
  Plan.ExternalROI = [Plan.CTname '_' remove_bad_chars(RTstruct.ExternalROI)];

  Plan.DoseGrid.size = handles.size';
  Plan.OptROIVoxels_nominal = createROI(Plan.ExternalROI, handles);
  Plan.OptROIVoxels_robust = []; % initialize to zeros

  %Create the ROI structure
  optFidx = 1;
  Plan.optFunction(optFidx) = createOptFunctionStruct();
  Plan.optFunction(optFidx).ROIname = Plan.ExternalROI;



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
