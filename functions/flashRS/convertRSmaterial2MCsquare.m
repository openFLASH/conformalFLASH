%----------------------------------------
% Convert the material file exported from the Raystation into
% the format required by MCsquare
%----------------------------------------

clear
close all

materialsPath = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\Output_scanner_calib\Materials' %Output folder where to save the material files
scannerPath = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\Output_scanner_calib\Scanners\Oncentra MasterPlan' %Output folder where to save the scanner calibration files
RSMaterialFile = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\scanner\material.txt' %File from Rasytation with the material composition
HU2densityconversion = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\scanner\HU-density.txt' %HU to density table as defiend in the Raystation.


RS_param = read_RS_param(RSMaterialFile) %REad the material composition stored in the Raystation
modelHU2rho = readRS_HU2rho(HU2densityconversion) %REad the table HU to density

%Create the MAterial files and the scanner calibration file
%from the Raystation export
exportRSConversion2MC2(modelHU2rho, RS_param, materialsPath, scannerPath)
