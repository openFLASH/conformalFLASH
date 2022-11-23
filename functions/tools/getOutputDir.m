%% getOutputDir
% Get the full path to the output directory for the beam |b|
% Ifthe folder does not exist, create it.
%
%% Syntax
% |res = help_header(im1,im2)|
%
%
%% Description
% |res = help_header(im1,im2)| Description
%
%
%% Input arguments
% |rootFolder| - _STRING_ - Full path to the root of the working folder
%
% |b| -_INTEGER_- Number of the beam in the plan
%
%
%% Output arguments
%
% |path2beamResults| - _STRING_ - Full path to the folder where the results for beam |b| shall be saved
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function path2beamResults = getOutputDir(rootFolder , b)

  path2beamResults = fullfile(rootFolder , 'Outputs' , ['Outputs_beam' , num2str(b)]);
  if (~exist(fullfile(path2beamResults),'dir'))
    %The folder to save the CT does not exist. Create it
    mkdir (fullfile(path2beamResults))
  end

end
