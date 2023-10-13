clear
close all

dataFolder = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\D-58';
planFileName = 'D:/programs/openREGGUI/REGGUI_userdata/raystation/D-58/FP-D58.dcm';
rtstructFileName = fullfile(dataFolder , 'RS1.2.752.243.1.1.20230220152044639.1200.30532.dcm');
CTname = fullfile('D:\programs\openREGGUI\REGGUI_userdata\raystation\D-58\reggui_CT' , 'reggui_CT_0001.dcm');
output_path = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\D-58_output';

RTstruct.ExternalROI = 'WaterCube'; %name for external ROI - the body contour

BeamProp.CEFDoseGrid = {1, 1, 1}; % Size (mm) of final dose scoring grid. Compute the final dose through CEF on a different grid than the high-res
BeamProp.protonsHighResDose = 1e6; %TODO Number of protons in the dose in high resolution CT
BeamProp.BDL = 'D:\programs\openREGGUI\flash\openMCsquare\lib\BDL\BDL_default_UN1_G0_Al_RangeShifter_tilted.txt'; %Identify the BDL file name from the treatment machine name in the plan
BeamProp.ScannerDirectory = 'D:\programs\openREGGUI\REGGUI\plugins\openMCsquare\lib\Scanners\default';
BeamProp.MCsqExecPath = 'D:\programs\openREGGUI\REGGUI\plugins\openMCsquare\lib';

CEMprop.makeSTL = false;


[handles, Plan] = flashLoadAndCompute(planFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct.ExternalROI, CEMprop);
