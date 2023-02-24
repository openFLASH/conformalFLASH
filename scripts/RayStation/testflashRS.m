clear
close all

dataFolder = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\Rombus-Plan-11B-flash';
planFileName = fullfile(dataFolder , 'Rhombus_600.dcm');
rtstructFileName = fullfile(dataFolder , 'Rombus_Plan_11B_flash_rtstruct.dcm');
CTname = fullfile('D:\programs\openREGGUI\REGGUI_userdata\raystation\Spiral_WaterPhantom\reggui_CT2' , 'reggui_CT2_0001.dcm');
output_path = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\Rhombus_out';

RTstruct.selected_ROIs = {'Shallow Sphere' }; %Name of the RT structs for which the dose rate is to be computed
RTstruct.DRCritical_ROIs = {'Shallow Sphere'}; %The OAR structure to include in the trajectory optimisation
RTstruct.ExternalROI = 'body'; %name for external ROI - the body contour
RTstruct.TargetROI = 'rombus'; %name for target ROI

DoseRate.Dref = 2; % Dose (Gy / FRACTION) in OAR above which the dose rate condition must be respected

BeamProp.NbScarves = 1; %umber of scarves to paint on the BEV
BeamProp.CEFDoseGrid = {1, 1, 1}; % Size (mm) of final dose scoring grid. Compute the final dose through CEF on a different grid than the high-res

CEMprop.makeSTL = true;


[handles, Plan] = flashLoadAndCompute(planFileName, CTname , rtstructFileName , output_path , BeamProp , RTstruct, DoseRate , CEMprop);
