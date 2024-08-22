%% getOptimumSpotSpacing
% Compute the spot spacing that will give the smallest dose variation in the largest field size
%
%% Syntax
% |[spotSpacing , spotSigma] = getOptimumSpotSpacing(handles , Plan , b , maxE )|
%
%
%% Description
% |[spotSpacing , spotSigma] = getOptimumSpotSpacing(handles , Plan , b , maxE )| Description
%
%
%% Input arguments
% |handles| - _struct_ - REGGUI data structure containing relevant information about the CT and RT DICOM plan.
%       It is used to compute the PBS spot at isocentre. The following data must be present:
%   * |handles.images| - _STRUCTURE_ - CT image data and RT struct data
%   * |handles.plans|  - _STRUCTURE_ - Plan data. New plan will be added to the structure
%
% |Plan| - _STRUCT_ - MIROpt structure where all the plan parameters are stored.
%
% |b| -_INTEGER_- Beam number
%
% |maxE| -_SCALAR_- Maximum energy (MeV) delivered by the treatment machine
%
%% Output arguments
%
% |spotSpacing| -_SCALAR_- Optimum distance (mm) between the spots
%
% |spotSigma| - _SCALAR_ - Sigma (mm) of the lateral Gaussian profile at the depth of the Bragg peak for the spot on the cen,tral axis
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [spotSpacing , spotSigma] = getOptimumSpotSpacing(handles , Plan , b , maxE )

  %Create a monolayer and single spot plan in order to estiumate the spot sigma
  Plan.Beams = Plan.Beams(b); %Keep one single beam only
  Plan.Beams.Layers.Energy = maxE; %One single energy, the maximum one
  Plan.Beams.Layers.SpotWeights = 1; %One single spot with a weight of 1
  Plan.Beams.Layers.SpotPositions = [0,0]; % The single spot is at isocentre
  Plan.fractions = 1; % use fraction = 1 for pseudo-plan. there is one single beam
  Plan.output_path = fullfile(Plan.output_path,'SingleSpot');
  if (~exist(Plan.output_path))
    mkdir (Plan.output_path)
  end

  %Compute the Pij matrix of  this single spot
  fprintf('Computating dose map for central spot \n')
  singleSpot = true;
  [~,Plan] = ComputePijsMatrix(Plan, handles, Plan.CTname, singleSpot); %Compute
  fprintf('Computation completed \n')
  OptConfig.BeamletsMatrixPrecision = 'd';
  Plan = ReadPijsMatrix(Plan, OptConfig, handles); %Load from disk
  fprintf('Dose map for central spot loaded \n')

  %Estimate the spot sigma at the Bragg peak of the maximum energy
  fprintf('Computing spot sigma at isocentre \n')
  [spotSigma , dist]= getSpotSigma(Plan , 1 , [0,0] , 'BP'); %mm Sigma of the Gaussian PBS spot near isocentre
  fprintf('Beam %d : sigma at BP = %.2f mm \n',b,spotSigma);
  fprintf('Beam %d : BP to isocentre distance  = %.2f mm\n',b,dist);

  % Compute the optimum spot spacing for the spot sigma
  [slope , intercp] = getLinearParam();
  spotSpacing = slope .* spotSigma + intercp;

end

%====================================================
% Get the parameters of the linear fit between the spot sigma
% and the optimumu  PBS spot spacing
%
% spacing = slope * spot sigma + intercep
%
% OUTPUT
% |slope| -_SCALAR_- Slope
% intercp| -_SCALAR_- Origin of the line (mm)
%====================================================
function [slope , intercp] = getLinearParam()

  %Parameters were obtained with the function distGauss1D.m
  slope = 1.487059;
  intercp = -0.088235;
end
