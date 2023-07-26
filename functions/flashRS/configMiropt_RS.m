%% configMiropt_RS
% Define the default parameters to load and process a Raysearch treatment plan
%
%% Syntax
% |Plan = configMiropt_RS(BeamProp, CEMprop, output_path)|
%
%
%% Description
% |Plan = configMiropt_RS(BeamProp, CEMprop, output_path)| Description
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
    Plan.SaveHighResDoseMap = true; % Do not save the dose map at CEFDoseGrid resolution in the IEC gantry CS
    Plan.SaveDoseBeamlets = true; % Do not save beamlet dose maps to save time
    Plan.SaveHighResCT = true; %Do not save the high resolution CT for each beamlet in the reference frame of the beamlet

    if ~isfield(Plan , 'MCsqExecPath')
      %No path to MCsqaure executable was provided. Add it now
      [~ , Plan.MCsqExecPath ] = get_MCsquare_folders();
    end

    %Beam properties
    %----------------
    BeamProp.FLAGOptimiseSpotOrder = false; %Do not optimise trajectory. Use the one read from logs
    if ~isempty(BeamProp)
      Plan = copyFields(BeamProp , Plan);
    end

    Plan.output_path = output_path;
    Plan.showGraph = true;
    if isfield(BeamProp , 'NbScarves')
      Plan.Extras.NbScarves = BeamProp.NbScarves;
    end
    Plan.Extras.addPrinterErrors = false; %Add printing error for CEM

    %CEM properties
    %--------------
    if isfield(CEMprop , 'makeSTL')
      Plan.makeSTL = CEMprop.makeSTL;
    else
      Plan.makeSTL = false;
    end
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
