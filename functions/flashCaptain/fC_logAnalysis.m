%% fC_logAnalysis
% Load the irradiation logs and compare the logs to the reference treatment plan:
%   * Compare the distance between spots and the Monitoring units of two PBS treatment plans
%   * Compute the dose map using records from the logs
%
%% Syntax
% |handles = fC_logAnalysis(configFile)|
%
%
%% Description
% |handles = fC_logAnalysis(configFile)| Description
%
%
%% Input arguments
% |configFile| -_STRING_- Full path and file name of the JSON with the scrip parameters
%
%
%% Output arguments
%
% |handles| -_STRUCT_- REGGUI structure
%   * |handles.plans| -CELL VECTOR of STRUCT- The treatment plan from TPS and the plan reconstructed from logs
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function   [handles, Plan] = fC_logAnalysis(configFile)


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

  %Load irradiation logs
  %---------------------
  [record_dir,record_file] = fileparts2(config.files.RecordName);
  %handles = Import_tx_records(record_dir,record_file,'iba','record',handles,'plan',config.files.AggregatePaintings);
  handles = Import_tx_records(record_dir,record_file,'iba','logs',handles,'plan',config.files.AggregatePaintings);


  %Compare plan and treatment record
  %Save in JSON file the spot position and MU differences
  %-------------------------------------------------------
  %Compare_plans('record' , 'plan' , handles , 0 , fullfile(config.files.output_path,config.files.ResultTXT) );

  %Compute the dose map
  %--------------------
  %logID = find(strcmp(handles.plans.name , 'record')); %identify the logs in handles
  logID = find(strcmp(handles.plans.name , 'logs')); %identify the logs in handles

  BeamProp.NbScarves = 1; %umber of scarves to paint on the BEV
  BeamProp.FLAGOptimiseSpotOrder = false; %Do not optimise trajectory. Use the one read from logs
  BeamProp.FLAGcheckSpotOrdering = false; %Check that spot ordering in plan matches scanAlgo output

  BeamProp = copyFields(config.BeamProp , BeamProp);
  BeamProp.CEFDoseGrid =  num2cell(BeamProp.CEFDoseGrid);
  CEMprop.makeSTL = false;

  %Load plan from TPS and create a MIROPT |PLan| structure with the monolayer plan
  [~, Plan] = flashLoadAndCompute(RSplanFileName, CTname , rtstructFileName , config.files.output_path , BeamProp , config.RTstruct.ExternalROI , CEMprop , [] , handles.plans.data{logID}{1});

end
