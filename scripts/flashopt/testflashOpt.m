clear
close all

dataFolder = 'D:\programs\openREGGUI\REGGUI_userdata\water_data\50_50_50';
planFileName = 'D:\programs\openREGGUI\REGGUI_userdata\water_data\50_50_50\Outputs\Outputs_beam1\Plan_CEM.dcm'
rtstructFileName = 'D:\programs\openREGGUI\REGGUI_userdata\water_data\water15\rtstructs_50_50_50.dcm';
CTname = 'D:\programs\openREGGUI\REGGUI_userdata\water_data\water15\water15_0001.dcm';
output_path = 'D:\programs\openREGGUI\REGGUI_userdata\water_data\50_50_50\reloaded';

RTstruct.ExternalROI = 'BODY'; %name for external ROI - the body contour

BeamProp.CEFDoseGrid = {1, 1, 1}; % Size (mm) of final dose scoring grid. Compute the final dose through CEF on a different grid than the high-res
BeamProp.protonsHighResDose = 1e4; %TODO Number of protons in the dose in high resolution CT
BeamProp.BDL = 'D:\programs\openREGGUI\flash\openMCsquare\lib\BDL\BDL_default_UN1_G0_Al_RangeShifter_tilted.txt'; %Identify the BDL file name from the treatment machine name in the plan
BeamProp.ScannerDirectory = 'D:\programs\openREGGUI\REGGUI\plugins\openMCsquare\lib\Scanners\default';
BeamProp.MCsqExecPath = 'D:\programs\openREGGUI\REGGUI\plugins\openMCsquare\lib';

CEMprop.makeSTL = false;


[handles, Plan] = flashLoadAndCompute(planFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct, CEMprop);
