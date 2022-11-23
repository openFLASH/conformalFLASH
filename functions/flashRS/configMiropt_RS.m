%% configMiropt_RS
% Define the default parameters to load and process a Raysearch treatment plan
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
% |im1| - _STRING_ -  Name
%
%
%% Output arguments
%
% |res| - _STRUCTURE_ -  Description
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Plan = configMiropt_RS(BeamProp, CEMprop, output_path)

    %Add the default MIROPT parameters
    %---------------------------------
    Plan.YAML.Plan = []; %Do not update any parameter in case of reload of partial results
    % if isfield(BeamProp,'version')
    %   Plan.version  = BeamProp.version;
    % end

    %MCsqaure properties
    %-------------------
    Plan.BeamletsBy = 'MCsquare'; %Algorithm used to compute the beamlets: Monte Carlo ('MCsquare') or Pencil Beam ('FoCa')
    Plan.protonsBeamlet = 1e5; % Number of protons per beamlet. Default 5e4
    Plan.protonsFullDose = 5e7; % Number of protons in the full target. Default 1e7
    Plan.protonsHighResDose = 5e4; %TODO Number of protons in the dose in high resolution CT

    Plan.ComputeCEMdose = true;  % TRUE = compute the dose distribution through the hedgehog, range shifter and aperture
    Plan.SaveDoseBeamlets = 'dcm'; % save the dose of each beamlet in the reference frame of the CT with aperture: dcm (DICOM format) , sparse (sparse matrix) , false (not saved)
    Plan.SaveHighResCT = true; %Do not save the high resolution CT for each beamlet in the reference frame of the beamlet
    Plan.SaveHighResDoseMap = true; % Do not save the dose map at CEFDoseGrid resolution in the reference frame of the beamlet
    Plan.CEFDoseGrid  = {1, 1, 1}; % Size (mm) of final dose scoring grid. Compute the final dose through CEF on a different grid than the high-res
    % if isfield(BeamProp , 'CEFDoseGrid')
    %   Plan.CEFDoseGrid = BeamProp.CEFDoseGrid;
    % end

    %Beam properties
    %----------------
    if ~isempty(BeamProp)
      Plan = copyFields(BeamProp , Plan);
    end

    Plan.output_path = output_path;
    Plan.showGraph = true;
    if isfield(BeamProp , 'NbScarves')
      Plan.Extras.NbScarves = BeamProp.NbScarves;
    end
    Plan.Extras.addPrinterErrors = false; %Add printing error for CEM
    %Plan.Scattering =  'Moliere'; %Scattering model used for CEM design. It relates the CEM material thickness to the spot sigma. POssible options are 'SAM', 'Moliere' and 'User'. At fault to define an algorithm, Moliere is used by default.

    %Aperture properties
    % if isfield(BeamProp , 'ApertureMargin')
    %   Plan.Extras.ApertureMargin = cell2mat(BeamProp.ApertureMargin);
    % end
    % Plan.Extras.ApertureBlock = 1; % 0 = do not use an aperture for this beam 1 = use an aperture for the beam

    %CEM properties
    %--------------
    % Spike = struct;
    % Spike = copyFields(CEMprop , Spike);
    % Plan.Spike = Spike;
    Plan.makeSTL = CEMprop.makeSTL;
    Plan.RidgeFilter = true;
    Plan.exportCEFinCT = false;

    %Count the number of robust scenario
    % There is only one nominal scenario for which we compute the dose rate
    Plan.SystSetUpError = 0;
    Plan.RangeError = 0;
    Plan.rr_nominal = 1;
    Plan.rs_nominal = 1;
    Plan.RandSetUpError = {[0,0,0]};
    Plan.Opt4Dmode = 0;
    Plan = countScenarios( Plan );
end
