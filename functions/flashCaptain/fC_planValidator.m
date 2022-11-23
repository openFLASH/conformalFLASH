%% fC_planValidator
% Load a DICOM RT plan and compute:
%  * The dose through the CEM
%  * The dose rate in the |RTstruct.selected_ROIs|
%
%% Syntax
% |[handles, Plan] = fC_planValidator(configFile)|
%
%
%% Description
% |[handles, Plan] = fC_planValidator(configFile)| Description
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

function [handles, Plan] = fC_planValidator(configFile)

% Load the YAML file with the parameters for the computation
%-----------------------------------------------------------
config = loadjson(configFile)


RSplanFileName = config.files.planFileName;
CTname = config.files.CTname;
rtstructFileName = config.files.rtstructFileName;
output_path = config.files.output_path;

BeamProp.NbScarves = 1; %umber of scarves to paint on the BEV
BeamProp.CEFDoseGrid = {1, 1, 1}; % Size (mm) of final dose scoring grid. Compute the final dose through CEF on a different grid than the high-res
BeamProp.FLAGOptimiseSpotOrder = false;
BeamProp.FLAGcheckSpotOrdering = false;

RTstruct = config.RTstruct;
DoseRate  = config.DoseRate;

CEMprop.makeSTL = true;

%Process the RayStation plan
%---------------------------
[handles, Plan] = flashLoadAndCompute(RSplanFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct, DoseRate , CEMprop);

end
