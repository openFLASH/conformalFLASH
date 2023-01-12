%% function
% Identify the MCsquare beam data library file for the specified machine name
%
%% Syntax
% |BDLname = BDLname4machine(MachineName)|
%
%
%% Description
% |BDLname = BDLname4machine(MachineName)| Description
%
%
%% Input arguments
% |MachineName| - _STRING_ - Name of the treatment machine as defined in the BDL
%
%
%% Output arguments
%
% |BDLname| Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function BDLname = BDLname4machine(MachineName)

  %Get the folder containing the BDLs
  [pluginPath , MCsqExecPath , BDLpath , MaterialsPath , ScannersPath] = get_MCsquare_folders();
  clist = dir(BDLpath); %List of all the BDL files

  NbIdentified = 0; %Number of BDL containing the machine name |MachineName|

  for fidx = 3:numel(clist) %skip . and ..
    BDL = clist(fidx).name; %Test this BDL file
    param = getMachineParam(BDL);
    if strcmp(param.MachineName,MachineName)
      %The |MachineName| is contained in this file
      BDLname = BDL;
      NbIdentified = NbIdentified +1;
    end
  end

  %Check that there is one and only one BDL containing this machine name
  if NbIdentified ~=1
    fprintf('Machine name = %s \n',MachineName);
    fprintf('Nb of BDL with that machine name = %d \n',NbIdentified);
    error('There should be one and only one BDL with this machine name');
  end

end
