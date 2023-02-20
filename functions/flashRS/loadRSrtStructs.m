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
% |DoseRate| -_STRUCTURE_- Definition of parameters for the computation of dose rates
%   * |DoseRate.Dref| -_SCALAR_- Threshold dose (Gy / FRACTION) in OAR above which the dose rate condition must be respected. Voxels below the threshold dose are not included in DR condition
%   * |DoseRate.DMF| -_SCALAR_- Dose modifying factor of flash for this organ. Not used for dose constraint.
%   * |DoseRate.DR50| -_SCALAR_- Gy/s 50% dose rate of the sigmoid probability function. Not used for dose constraint.
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

function [handles , Plan , ROI , usedROI] = loadRSrtStructs(rtstructFileName, handles, Plan, RTstruct, DoseRate)

if nargin < 5
    DoseRate = [];
end

%Load the RT struct
NbStruct = numel(RTstruct.selected_ROIs);
strName = unique({RTstruct.selected_ROIs{:}, RTstruct.TargetROI, RTstruct.ExternalROI}); %Make sure there is no duplicate of name so that we load the structure only once
handles = Import_contour(rtstructFileName,strName,Plan.CTname,1,handles);

% Create ROI specific properties for all structures
Plan.ExternalROI = [Plan.CTname '_' remove_bad_chars(RTstruct.ExternalROI)];
Plan.TargetROI = [Plan.CTname '_' remove_bad_chars(RTstruct.TargetROI)];
Plan.TargetROI_ID = numel(strName) - 1;
Plan.ExternalROI_ID = numel(strName);

Plan.DoseGrid.size = handles.size;
nvoxels = prod(handles.size);
Plan.OptROIVoxels_nominal = false(nvoxels,1); % initialize to zeros
Plan.OptROIVoxels_robust = false(nvoxels,1); % initialize to zeros
ROItotalMask = zeros(Plan.DoseGrid.size'); %Create a structure which is the sum of all ROI

%Create the ROI structure
for i = 1:numel(strName)
    % Set ROI name and create list of ROIs
    RTstruct.selected_ROIs{i} = [Plan.CTname '_' remove_bad_chars(strName{i})];
    ROImask = Get_reggui_data(handles, RTstruct.selected_ROIs{i});
    ROI(i) = createROI(RTstruct.selected_ROIs{i}, ROImask , handles);

    if (sum(strcmp(strName{i} , RTstruct.DRCritical_ROIs)))
        %This is a dose rate critical structure to include in the spot trajectory
        Plan.optFunction(i) = createOptFunctionStruct(remove_bad_chars(strName{i}) , ROI , DoseRate , true);
      else
        %This is not a dose rate critical structure. Ignore it for the trajectory computation
        Plan.optFunction(i) = createOptFunctionStruct(remove_bad_chars(strName{i}) , ROI , DoseRate , false);
    end

    if(i<numel(strName)) %Do not add EXTERNAL to the mask
        Plan.OptROIVoxels_nominal = Plan.OptROIVoxels_nominal | ROI(i).mask1D;
        ROItotalMask = ROItotalMask  | ROImask;
    end
end

%Add the mask for the total volume in which the dose rate is to be computed
TotalMaskName = [Plan.CTname '_' remove_bad_chars('TotalMask')];
RTstruct.selected_ROIs{end+1} = TotalMaskName;
handles = Set_reggui_data(handles,TotalMaskName,ROItotalMask , [] , 'images');
ROImask = Get_reggui_data(handles, TotalMaskName);
ROI(end+1) = createROI(TotalMaskName, ROImask , handles);
Plan.optFunction(end+1) = createOptFunctionStruct(remove_bad_chars('TotalMask') , ROI , DoseRate , false);

Plan.TargetROI_ID = getROIByName(ROI,  [Plan.CTname '_' remove_bad_chars(RTstruct.TargetROI)]);

%Get the list of structure names
usedROI = RTstruct.selected_ROIs;

end

%------------------------------------
% Create a sttructure defining the dose rate objective function
%------------------------------------
function optFunction = createOptFunctionStruct(ROIName , ROI , DoseRate , isInTrajectory)

  if ~isInTrajectory
    optFunction.ID = 1; %This is not an objective on dose rate. Ignore from trajectory computation
  else
    optFunction.ID = 8; %Select this ID because the optimisation function ID=8 requires the computation of dose rate
  end

  optFunction.ROIindex = getROIByName(ROI, ['ct_' , ROIName]);
  optFunction.ROIname = ROIName;

  if isfield(DoseRate,'Dref')
    optFunction.Dref  = DoseRate.Dref;
  end

  if isfield(DoseRate,'DMF')
    optFunction.DMF  = DoseRate.DMF;
  else
    optFunction.DMF = [];
  end
  if isfield(DoseRate,'DR50')
    optFunction.DR50  = DoseRate.DR50;
  else
    optFunction.DR50 = [];
  end
end


%---------------------------
% add the ROI to thje list of ROIs
%-----------------------------
function ROI = createROI(ROIname, ROImask , handles)

    ROI.name = ROIname;
    % Calculate mask for each ROI
    ROI.mask3D.value = ROImask;

    % Convert to 1D sparse mask with z inverse format (to be consistent
    % with beamlets format)
    temp = flip(ROI.mask3D.value,3);
    ROI.mask1D = sparse(logical(double(temp(:))));

  end
