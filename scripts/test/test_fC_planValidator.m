clear
close all

JSONfileName = 'D:\programs\openREGGUI\flashTPS\scripts\test\fC_planValidator.json'

[handles, Plan] = fC_planValidator(JSONfileName);

JSONfileName = 'D:\programs\openREGGUI\flashTPS\scripts\test\fC_planValidator_indicators.json'
handles = fC_ClinicalIndicators(JSONfileName);
