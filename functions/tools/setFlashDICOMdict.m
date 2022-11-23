%% setFlashDICOMdict
% Define the DICOM dictionary containing the private DICOM tags defined for FLASH
%
%% Syntax
% |dictionary = setFlashDICOMdict()|
%
%
%% Description
% |dictionary = setFlashDICOMdict()| Description
%
%
%% Input arguments
% None
%
%
%% Output arguments
%
% |dictionary| - _STRING_ - Return the full path to the file with the private DICOM dictionary
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function dictionary = setFlashDICOMdict()

  currentDict = dicomdict('get'); %GEt the current DICOM dictionary

  %Define the special FLASH DICOM dictionary
  regguiPath = fileparts(which('reggui'));
  dictionary = fullfile(regguiPath ,'plugins','openMIROpt','functions','io','dicom-dict.txt');

  %If required, change the dictionary
  if ~strcmp(currentDict,dictionary)
    fprintf('Setting DICOM dictionary : %s n',dictionary)
    dicomdict('set',dictionary)
  end


end
