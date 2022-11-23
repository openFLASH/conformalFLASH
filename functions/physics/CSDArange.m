%% CSDArange
% Compute the proton range using the continuous slowing down approximation (CSDA).
% The mass stopping power is extracted from the G4_Stop_Pow.dat file in the MCsquare library
%
%% Syntax
% |[range , Ep] = CSDArange(material , Ep)|
%
% |[range , Ep] = CSDArange(material)|
%
%
%% Description
% |[range , Ep] = CSDArange(material , Ep)| Description
%
%
%% Input arguments
% |material| - _STRING_ - Name of the material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%
% |E| -_SCALAR VECTOR_- [OPTIONAL. Default = E=[]] Energy (MeV) of the incoming proton beam
%
%
%% Output arguments
%
% |range| -_SCALAR VECTOR_- |range(i)| Range (cm) in the material of a proton beam with energy E(i)
%
% |Ep| -_SCALAR VECTOR_- |Ep(i)| The i-th energy (MeV) of the proton beam at which the range is computed. This is |E| if it is provided as input
%                   % If E=[], then the range is computed at all enrgy defined in G4_Stop_Pow.dat file
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [range , Ep] = CSDArange(material , Ep)

  if nargin < 2
    Ep = [];
  end

  [Name , MaterialsDirectory] = getMaterialNameinDB(material);

  Density = Import_MatProperty(fullfile(MaterialsDirectory, Name, 'Material_Properties.dat'),'Density');
  [E, SP] = Import_SP_data(fullfile(MaterialsDirectory, Name, 'G4_Stop_Pow.dat')); %mass stopping power as a function of energy

  %Remove the stopping power equal to zero to avoid infinite 1/SP
  idxNull = find(SP==0);
  E(idxNull) = [];
  SP(idxNull)= [];

  dE_dx = Density .* SP; %MeV/cm linear stopping power
  CSDAr = cumtrapz(E,1./dE_dx); %proton range (cm) in the material for energy E  using the CSDA approximation

  if ~isempty(Ep)
    range = interp1(E,CSDAr,Ep);
  else
    %No energy was provided. Compute range for all energies reported in the G4 database
    Ep = E;
    range = CSDAr;
  end


end

%===================================================
% Read the Geant 4 stopping power calibration curve
% and interpolate for therequested energies
% Author: K. Souris (Modifications: R. Labarbe)
%===================================================
function [E, SP] = Import_SP_data(FileName)

  fid=fopen(FileName,'r');
  if (fid < 0)
      error(['Unable to open file ' FileName])
  end

  data = fscanf(fid, '%f %f', [2 inf]);
  fclose(fid);

  E = data(1,:);
  SP = data(2,:);

end

%=====================================
% Get the correct material name to retrieve Geant4 database
%=====================================
function [Name , MaterialsDirectory] = getMaterialNameinDB(material)

  MaterialsDirectory = fullfile(getPluginPath('openMCsquare'), 'lib', 'Materials');
  Material_List_File = fullfile(MaterialsDirectory, 'list.dat');

  % Import the list of materials
  fid = fopen(Material_List_File);
  Material_List_File = textscan(fid, '%d %s', 'delimiter', '\n', 'MultipleDelimsAsOne',1);
  fclose(fid);
  Material_index = Material_List_File{1};
  Material_name = Material_List_File{2};

  %Check that the requested material is in the list
  checkList = strfind(Material_name,material);
  found =0;
  for idx = 1:numel(checkList)
    if ~isempty(checkList{idx})
      found = idx;
      break; %We have a match
    end
  end
  if (idx==0)
    material
    error('Material not found')
  end

  Name = strsplit(Material_name{found}, {' ', '\t', '#'});
  Name = Name{1};

end

%================================
% Read material properties from material file
% Author: K. Souris (Modifications: R. Labarbe)
%================================
function value = Import_MatProperty(FileName , property)

    fid=fopen(FileName,'r');
    if (fid < 0)
        error(['Unable to open file ' FileName])
    end

    while ~feof(fid)
        tmp = fgetl(fid);
        key = strsplit(tmp, {' ', '\t'});
        id = find(strcmp(key, property));
        if (not(isempty(id)))
            value = key{id+1};
            value = str2num(value);
        end
    end

    fclose(fid);
end
