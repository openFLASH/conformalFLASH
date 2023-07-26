%% loadRSctScan
% Load the CT scan that was used to create the RT plan in hte Raystation
% Update the REGGUI handles and the PLan structure.
%
%% Syntax
% |[handles, Plan] =  loadRSctScan(CTname , handles, Plan)|
%
%
%% Description
% |[handles, Plan] =  loadRSctScan(CTname , handles, Plan)| Description
%
%
%% Input arguments
%
% |CTname| -_STRING_- File name and full path to the CT scan
%                   If empty, the CT scan must be provided in handles
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
%
%% Output arguments
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are stored.
%
% |handles| - _STRUCT_ - REGGUI data structure.
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [handles, Plan] =  loadRSctScan(CTname , handles, Plan)

if isempty(handles)
    handles = struct;
    handles.path = Plan.output_path;
    handles = Initialize_reggui_handles(handles);
end

if ~isempty(CTname)
  %CT scan file name is provided. Load it from disk
  Plan.CTname = 'ct';
  [CTdirectory,CTfileName,EXT] = fileparts(CTname);
  handles = Import_image(CTdirectory,[CTfileName EXT],1,Plan.CTname,handles);
end

[~,info] = Get_reggui_data(handles,Plan.CTname);

%Read the scanner model from the plan in order to identify the correct scanner calibration folder for MCsquare
fprintf('CT scanner calibration file : %s \n', Plan.ScannerDirectory);

Plan.CTinfo = info;
Plan.CTinfo.Spacing = info.Spacing;
Plan.DoseGrid.resolution = Plan.CTinfo.Spacing;
Plan.CTinfo.ImagePositionPatient = info.ImagePositionPatient;

end
