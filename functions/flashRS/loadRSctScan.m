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

CTimageName = 'ct';
[CTdirectory,CTfileName,EXT] = fileparts(CTname);
handles = Import_image(CTdirectory,[CTfileName EXT],1,CTimageName,handles);
[~,info] = Get_reggui_data(handles,CTimageName);

%Read the scanner model from the plan in order to identify the correct scanner calibration folder for MCsquare
if isfield(info.OriginalHeader , 'ManufacturerModelName')
  Plan.ScannerDirectory = remove_bad_chars(info.OriginalHeader.ManufacturerModelName);
else
  warning('(008,1090) ManufacturerModelName undefined in CT scan')
  Plan.ScannerDirectory = 'default';
  fprintf('Using default folder for CT scanner calibration : %s \n', Plan.ScannerDirectory);
end
fprintf('CT scanner calibration file : %s \n', Plan.ScannerDirectory);

%Check that the scanner directory exists.
[pluginPath , MCsqExecPath , BDLpath , MaterialsPath , ScannersPath] = get_MCsquare_folders();
if (~exist(fullfile(ScannersPath, Plan.ScannerDirectory),'dir'))
  %TODO This should be an error and the program should stop here in a well configured system
  warning('(008,1090) ManufacturerModelName references an unknown scanner')
  fprintf('Folder for CT scanner calibration : %s \n', Plan.ScannerDirectory);
  Plan.ScannerDirectory = 'default';
  fprintf('Using default folder for CT scanner calibration : %s \n', Plan.ScannerDirectory);
end

Plan.CTinfo = info;
Plan.CTname = CTimageName;

Plan.CTinfo.Spacing = info.Spacing;
Plan.DoseGrid.resolution = Plan.CTinfo.Spacing;
Plan.CTinfo.ImagePositionPatient = info.ImagePositionPatient;

end
