%% getMachineParam
% Get the specification of the treatment machine
% The specs are read from the DBL file
% Somme additional specs are defined within this function.
%
%% Syntax
% |param = getMachineParam(BDL)|
%
%
%% Description
% |param = getMachineParam(BDL)| Description
%
%
%% Input arguments
% |BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
%
%% Output arguments
%
% |param| - _STRUCTURE_ - Proton machine parameters
%   * |param.MachineName| - _STRING_ - Name of the treatment machine
%   * |param.MachineType| - _STRING_ - Description of the treatment machine (PROTEUSone , PROTEUSplus)
%   * |param.MAXcurrent| -_SCALAR_- Maximum average beam current (uA) in nozzle
%   * |param.MAXenergy| -_SCALAR_- Beam energy (MeV) of deepest layer. Defines the beam line transmisssion for MAXcurrent in nozzle
%   * |param.ChargeFraction| -_SCALAR VECTOR_- |param.ChargeFraction(i)| Maximum fraction of the total charge delivered in the i-th burst
%   * |param.PulsePeriod| -_SCALAR_- Period (s) between pulses in a burst
%   * |param.ScanSpeed| -_SCALAR_- Scanning speed (mm/s) to move from one spot location to another spot location
%   * |param.snout| -_STRUCTURE_- Description of the FLASH accessory holder. See getParamSnout.m
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function param = getMachineParam(BDL)

  persistent BDLdata; %Store the BDLdata data in a persistent variable. Later function calls are faster because there will be no need to read again from table

  if ~isstruct(BDLdata)
    BDLdata = struct;
  end
  if ~isfield(BDLdata , 'BDL')
    BDLdata.BDL = '';
  end

  if ~strcmp(BDLdata.BDL , BDL)
      %This is the first time the function is called with this BDL name
      %Load BDL from disk
      [~ , sadX , sadY] = get_sad(BDL);
      param.VDSA = [sadX , sadY];  %mm source to axis distance for the 2 scanning magnets
      [param.MachineType , param.MachineName] = getMachineFromBDL(BDL);

      [regguiRoot,~,~] = fileparts(which('reggui.m'));
      BDLpath = fullfile(regguiRoot, 'plugins','openMCsquare', 'lib', 'BDL', BDL);
      BDLdata.data = load_BDL(BDLpath);
      param.MAXenergy = max(BDLdata.data.NominalEnergy);  %MeV Beam energy of deepest layer. Defines the beam line transmisssion for MAXcurrent in nozzle

      BDLdata.param = param;
      BDLdata.BDL = BDL;

    else
      %The BDL was previously loaded from disk. Reuse it
      param = BDLdata.param;
    end %if ~isfield

  %Ifthe machine type is missing, break here
  if isempty(param.MachineType)
    return
  end

  switch param.MachineType
    case 'PROTEUSplus'
        param.MAXcurrent = 1; %uA
        param.PulsePeriod = 0; %s
        param.ScanSpeed = 8000; %Approximate scanning speed (mm/s) to move from one spot location to another spot location
        param.MinimumSpotDuration = 5 .* 2.5e-4; %'ms) Minimum time required to deliver a spot with sufficient dose accuracy
            % For P+, collect 5 cycles of 250us (2.5e-4s) at least, when registering the charge and to have a good measure
            % The delivery of a PBS spot requires at least 5 cycle. ScanAlgo reduce the spot current to be sure to deliver the spot in 5 cycles

        %Parameters of the FLASH accessory holder
        param.snout = getParamSnout('flash-UN80');

    case 'PROTEUSone'
        param.MAXcurrent = 1; %uA
        %The total spot charge is delivered in several bursts. One burst is composed of several pulses. The period of one pulse is defined in |param.PulsePeriod|
        param.PulsePeriod = 1e-3; %s Period between pulses in a burst
        param.ScanSpeed = 8000; %Approximate scanning speed (mm/s) to move from one spot location to another spot location

        %Parameters of the FLASH accessory holder
        param.snout = getParamSnout('flash-UN80');

    otherwise
      MachineType
      error('Unknown machine type')
  end

end
