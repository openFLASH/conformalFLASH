%---------------------------------------------------------------------------
% Process the Machine log file to compute the dose and dose rate
% Use the TWIGA scripts in REGGUI.
%
% Change the configuration parameters in the two JSON file:
% 1) scripts\twiga\config_twiga_iba_log_converter.json
% 2) scripts\twiga\config_twiga_mc_computation.json
% to describe thje location of the data on your computer and to define the MCsqaure parameters
%---------------------------------------------------------------------------

%Folder where the configuiration JSON files are located.
% Change it to match the path on your computer
dataPath = 'D:\programs\github\openFLASH\conformalFLASH\scripts\twiga';


% STEP 1: Convert the logs from ZIP to CSV
% Call the TWIGA function that converts the machine logs (ZIP format)
% into a text CSV file
%----------------------------------
callReggui(fullfile(dataPath,'config_twiga_iba_log_converter.json'));


%STEP 2 : Process the logs
%  Call the Twiga function that process the logs and compute the dose and dose rate
% using MC square
%---------------------
callReggui(fullfile(dataPath,'config_twiga_mc_computation.json'));
