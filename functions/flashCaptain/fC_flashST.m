%% fC_flashST
% Load a DICOM RT plan for a shoot-through irradiation and compute:
%  * The dose for a shoot throught irradiation
%  * The dose rate
%
%% Syntax
% |[handles, Plan] = fC_flashST(configFile)|
%
%
%% Description
% |[handles, Plan] = fC_flashST(configFile)| Description
%
%
%% Input arguments
%
% |configFile| -_STRING_- Full path and file name of the JSON with the scrip parameters
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

function [handles, Plan] = fC_flashST(configFile)

  % Load the JSON file with the parameters for the computation
  %-----------------------------------------------------------
  config = loadjson(configFile)

  RSplanFileName = config.files.planFileName;
  CTname = config.files.CTname;
  rtstructFileName = config.files.rtstructFileName;


  %Load the data
  %-------------
  handles = struct;
  handles.path = config.files.output_path;
  handles.dataPath = config.files.output_path;
  handles = Initialize_reggui_handles(handles);


  %Load reference plan from TPS
  %----------------------------
  [plan_dir,plan_file] = fileparts2(config.files.planFileName);
  handles= Import_plan(plan_dir,plan_file,'dcm','plan', handles);

  %Update beam information with provided input parameters
  BeamProp.NbScarves = 1; %umber of scarves to paint on the BEV
  BeamProp.FLAGOptimiseSpotOrder = false; %Do not optimise trajectory. Use the one read from logs
  BeamProp.FLAGcheckSpotOrdering = false; %Check that spot ordering in plan matches scanAlgo output

  BeamProp = copyFields(config.BeamProp , BeamProp);
  BeamProp.CEFDoseGrid =  num2cell(BeamProp.CEFDoseGrid);


  %Load plan from TPS and create a MIROPT |PLan| structure with the monolayer plan
  CEMprop = struct;
  [~, Plan] = flashSTLoadAndCompute(RSplanFileName, CTname , rtstructFileName , config.files.output_path , BeamProp , config.RTstruct );

end
