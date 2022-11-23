%% fC_ClinicalIndicators
% Compute, for the provided dose or doserate maps:
%   * The clinical indicators defined in the JSON |config.files.indicators|. If no indicator are provided, then nothing is computed.
%   * The DVH for all the beams and all the structures present in the RT strcut file
%   * The gamma index between the test map and a reference map for for all provided beams
%
% Save the results in the output path as JSON files (the clinical indicators are saved at the format defined in save_Indicators.m)
%
%% Syntax
% |[indicators,handles] = fC_ClinicalIndicators(configFile)|
%
%
%% Description
% |[indicators,handles] = fC_ClinicalIndicators(configFile)| Description
%
%
%% Input arguments
% |configFile| -_STRING_- Full path and file name of the JSON with the scrip parameters. See at the bottom of this .m file for a description of the format of the JSON file.
%
%
%% Output arguments
%
% |handles| -_STRUCT_- REGGUI structure
%   * |handles.images| -_STRUCTURE_- The dose map and the RT structs to usde to compute clinical indicators
%   * |handles.indicators| -_STRUCTURE_- The results of the clinical indicators computations. See Compute_indicators.m for structure format
%   * |handles.dvhs| -_STRUCTURE_- The dose volume histogram on the dose map and dose rate maps for all structures used in clnical indicators
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function handles = fC_ClinicalIndicators(configFile)

  % Load the JSON file with the parameters for the computation
  %-----------------------------------------------------------
  config = loadjson(configFile)

  %Load the data
  handles = struct;
  handles.path = config.files.output_path;
  handles = Initialize_reggui_handles(handles);

  %Load CT scan
  CTimageName = 'ct';
  [CTdirectory,CTfileName,EXT] = fileparts(config.files.CTname);
  handles = Import_image(CTdirectory,[CTfileName EXT],1,CTimageName,handles);

  %Load all the RT structs
  [handles , contour_names] = Import_contour(config.files.rtstructFileName , 'all' , CTimageName , 1 , handles);

  %Loop for every JSON defining clinical indicators and make the computations
  for idx = 1:numel(config.files.indicators)
    handles  = processMaps(handles , config.files.TESTMap{idx} , config.files.indicators{idx} , CTimageName ,  config.files.output_path , config.files.REFMap{idx} , config.GammaIndex , contour_names);
  end

  %Save the DVH in a JSON file
  dvh_data = handles.dvhs;
  savejson('' , dvh_data , fullfile(config.files.output_path,'dvh_data.json') );

end


%----------------------------
% Compute all the indicators for the specififed dose (rate) maps
%
% INPUT
% |handles| -_STRUCT_- REGGUI structure
%
% |doseMap| -_CELL VECTOR of STRINGS_- |doseMap{b}| Full path and file name of the dose map for b-th beam
%
% |indicatorFilenName| -_STRING_- Full path and filename to the JSON file containing the clinical indicaotr. The file has the format defined in save_Indicators.m
%
% |CTimageName| -_STRING_- Name of the CT image in |handles.images|
%
% |REFdosemap| -_CELL VECTOR of STRINGS_- |doseMap{b}| Full path and file name of the REFERENCE dose map for b-th beam (for gamma index computation)
%
% |GammaIndex| -_STRUCTURE_- Parameter for the gamma index computation
%     * |GammaIndex.DistanceTolerance| _SCALAR_ Distance tolerance in mm (Default = 1mm)
%     * |GammaIndex.DoseTolerance| _SCALAR_ Dose tolerance in % (Default = 1%)
%
% |contour_names| -_CELL VECTOR of STRING_- List of the RT struct name in |handles|
%
% OUTPUT
% |handles| -_STRUCT_- REGGUI structure
%   * |handles.images| : The dose map used to compute clinical indicators
%   * |handles.indicators| -_STRUCTURE_- The results of the clinical indicators computations. See Compute_indicators.m for structure format
%   * |handles.dvhs| -_STRUCTURE_- The dose volume histogram on the dose map and dose rate maps for all structures used in clnical indicators
%
%----------------------------
function handles = processMaps(handles , doseMap , indicatorFilenName , CTimageName , output_path , REFdosemap , GammaIndex, contour_names)

  %Load maps for each beam
  DoseImageName = [];
  TotNbBeams = numel(doseMap); %Total number of dose map provided

  for idx = 1:numel(doseMap)
    [Dodirectory , DofileName , EXT] = fileparts(doseMap{idx});
    fprintf('File for beam %d : %s \n' , idx , [DofileName EXT]);
    [handles , DoseImageName{end+1} ] = Import_image(Dodirectory, [DofileName EXT] , 1 , [DofileName num2str(idx)] , handles);
  end
  MapName = DofileName; %The filename of the indicator results will be the one of the last beam

  if numel(DoseImageName) == 1
    %There is only one beam provided
    %Virtually create a second beam by copying the first beam so that Compute_indicators.m will assume that it receives dose map per beam
    DoseImageName{end+1} = DoseImageName{end};
  end

  if exist(indicatorFilenName) == 2
      %The file defining the indicator is a valid file name.
      % Compute the indicators on dose
      indicator_name = ['indicators_results_' MapName];
      [indicators , handles ] = makeIndicatorComputation(handles , indicator_name , CTimageName , {DoseImageName{:}} , indicatorFilenName );

      resultFileName = fullfile(output_path , [indicator_name , '.json']);
      save_Indicators(indicators,resultFileName,'json');

    else
      fprintf('Did not receive file to clinical indicator. \n')
  end

  % compute DVH on all RT structures
  show_dvh = false; %silent computation
  handles = DVH_computation(handles , show_dvh , {DoseImageName{1:TotNbBeams}} , contour_names);


  %Load reference dose maps
  REFDoseImageName = [];
  for idx = 1:numel(REFdosemap)
    [Dodirectory , DofileName , EXT] = fileparts(REFdosemap{idx});
    [handles , REFDoseImageName{end+1}] = Import_image(Dodirectory,[DofileName EXT] , 1 , ['REF' DofileName num2str(idx)] , handles);
  end

  if numel(REFDoseImageName) == 1
    %There is only one bem in this plan
    %Virtually create a second beam by copying the first beam so that Compute_indicators.m will assume that it receives dose map per beam
    REFDoseImageName{end+1} = REFDoseImageName{end};
  end

  %Compute the gamma index in a mask that is the uniuon of all RT structures
  options(1) = GammaIndex.DistanceTolerance; % _SCALAR_ Distance tolerance in mm (Default = 1mm)
  options(2) = GammaIndex.DoseTolerance; % _SCALAR_ Dose tolerance in % (Default = 1%)
  handles = Addition(contour_names, 'mask_AllStruct'  , handles); %Compute gamma index only inside all the loaded RT structures

  for idx = 1:TotNbBeams

        fprintf('Computing gamma index on map %s \n',DoseImageName{idx})
        GammaIndexName = ['GammaIndex_results_' DoseImageName{idx}];

        [handles , ~ , passing_rate , average_gamma] = Gamma_index(DoseImageName{idx} , REFDoseImageName{idx} , 'mask_AllStruct' , GammaIndexName , handles , options);

        fprintf('Passing rate  : %f \n', passing_rate(1))
        fprintf('Average gamma : %f \n', average_gamma(1))

        GammaIndexResults = struct;
        GammaIndexResults.TestMap = DoseImageName{idx};
        GammaIndexResults.ReferenceMap = REFDoseImageName{idx};
        GammaIndexResults.passing_rate = passing_rate(1); %The first element is gamma for all computing region. The second element is whithin the mask
        GammaIndexResults.average_gamma = average_gamma(1); %The first element is gamma for all computing region. The second element is whithin the mask

        %Save gamma index to disk
        handles = Export_image(GammaIndexName,fullfile(output_path, GammaIndexName) , 'dcm' , handles);
        savejson('' , GammaIndexResults , fullfile(output_path , [GammaIndexName , '.json']) );
  end

end


%--------------------------------
% Compute the clinical indicators
%  * Load the JSON file with  definition of clinical indicators
%  * Make the computation of the requested indicators
%
% INPUT
%
% |handles| -_STRUCT_- REGGUI structure
%   * |handles.images| : The dose map and the RT structs to usde to compute clinical indicators
%
% |indicator_name| -_STING_- Name of the indicator set to use to store in |handles.indicators|
%
% |CTimageName| -_STRING_- Name of the CT image in |handles.images|
%
% |DoseImageName| -_CELL VECTOR of STRING_- |DoseImageName{b}| is the name in |handles.images|  of the dose map for the b-th beam
%
% |indicatorFilenName| -_STRING_- Full path and filename to the JSON file containing the clinical indicaotr. The file has the format defined in save_Indicators.m
%
%
% OUTPUT
% |indicators| -_STRUCTURE_- The results of the clinical indicators computations. See Compute_indicators.m for structure format
%
% |handles| -_STRUCT_- REGGUI structure
%   * |handles.images| -_STRUCTURE_- The dose map and the RT structs to usde to compute clinical indicators
%   * |handles.indicators| -_STRUCTURE_- The results of the clinical indicators computations. See Compute_indicators.m for structure format
%
%--------------------------------

function [indicators , handles ] = makeIndicatorComputation(handles , indicator_name , CTimageName , DoseImageName , indicatorFilenName )

  [indicators_dir,indicators_file , ext] = fileparts(indicatorFilenName);
  handles = Import_indicators(indicators_dir,[indicators_file , ext],'json',indicator_name,handles);

  indIdx = find(strcmp(handles.indicators.name , indicator_name));
  selected_ROIs = {};
  for idx = 1:numel(handles.indicators.data{indIdx})
      selected_ROIs{end+1} = [handles.indicators.data{indIdx}{idx}.struct];
  end
  selected_ROIs = unique(selected_ROIs); %Make sure there is no duplicate of structure names

  %Define the nema of the structures on which to compute clinical indicators
  strName = {selected_ROIs{:}};
  for idx = 1:numel(handles.indicators.data{2})
      strName{idx} = [CTimageName '_' strName{idx}];
  end
  strName = {strName{:}};


 %Compute the clinical indicators
 %------------------------------
 [indicators,handles] = Compute_indicators(handles,indicator_name,strName,DoseImageName);

end




%=============================================================================
% Description of the format of the input JSON file
%=============================================================================
% {
% "files": {
%   "planFileName" : "path & file name of the DICOM RT ion plan",
%   "rtstructFileName" : "path & file name of the RT struct file",
%   "CTname" : "path & file name of one slice ofthe CT scan",
%   "output_path" : "Path where to save the results",
%         #Several computation of indicator, gamma index and DVH can be done sequentially on different dose and dose rate map
%         #For example, one on the dose maps and one on the dose rate map
%
%   "indicators" : [
%                 "path & file name to the indicator JSON file for the 1st compuation set" ,
%                 "path & file name to the indicator JSON file for the 2nd compuation set",
%                 "" #It is possible to skip the ocmputation of the indicator in one set by using an empty string. Only gamma index and DVH will be computed
%                 ],
%         #For each item of "indicators", there must be a corresponding item defining the test map (dose or dose rate).
%         #If there are several beams in one set, then the file name of each beam is to be provided
%   "TESTMap" : [
%                 [ #1st computation set
%                   "path & file name to dose (rate) map of the 1st beam for the 1st computation set"
%                 ],
%                 [ #2nd computation set
%                   "path & file name to dose (rate) map of the 1st beam for the 2nd computation set",
%                   "path & file name to dose (rate) map of the 2nd beam for the 2nd computation set"
%                 ],
%                 [ #3rd computation set
%                   "path & file name to dose (rate) map of the 1st beam for the 3rd computation set"
%                 ]
%               ],
%         #For each item of "indicators", there must be a corresponding item defining the reference map (dose or dose rate).
%         #If there are several beams in one set, then the file name of each beam is to be provided
%   "REFMap" :  [
%                 [ #1st computation set
%                   "path & file name to reference dose (rate) map of the 1st beam for the 1st computation set"
%                 ],
%                 [ #2nd computation set
%                   "path & file name to reference dose (rate) map of the 1st beam for the 2nd computation set",
%                   "path & file name to reference dose (rate) map of the 2nd beam for the 2nd computation set"
%                 ],
%                 [ #3rd computation set
%                   "path & file name to reference dose (rate) map of the 1st beam for the 3rd computation set"
%                 ]
%               ]
%   },
%
%  # PArameter of the gamma index computation
% "GammaIndex": {
%   "DistanceTolerance" : 3,
%   "DoseTolerance" : 3
%   }
%
% }
