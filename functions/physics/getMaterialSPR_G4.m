%% getMaterialSPR_G4
% Retrieve the material stopping pwoer from the MCsqaure material library
%
%% Syntax
% |[SPR, SP, Density, RelElecDensity] = getMaterialSPR_G4(material , E)|
%
%
%% Description
% |[SPR, SP, Density, RelElecDensity] = getMaterialSPR_G4(material , E)| Description
%
%
%% Input arguments
% |material| - _STRING_ - Name of the material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%
% |E| -_SCALAR_- Energy (MeV) of the incoming proton beam
%
%% Output arguments
%
% |Density| -_SCALAR_- Density (g/cm3) of the material
%
% |SPR(i)| -_SCALAR_- Relative stopping power (relative to water) at energy E(i) SPR(i) = Density * SP(i) / Water_SP;
%
% |SP(i)| -_SCALAR_- Mass stopping powers  (MeV cm2/g) at energy E(i)
%
% |RelElecDensity| -_SCALAR_- Electron density relative to water
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com) from the code of K. Souris

function [SPR, SP, Density, RelElecDensity] = getMaterialSPR_G4(material , E)

  MaterialsDirectory = fullfile(getPluginPath('openMCsquare'), 'lib', 'Materials');
  Material_List_File = fullfile(MaterialsDirectory, 'list.dat');

  % Import the list of materials
  fid = fopen(Material_List_File);
  Material_List_File = textscan(fid, '%d %s', 'delimiter', '\n', 'MultipleDelimsAsOne',1);
  fclose(fid);
  Material_index = Material_List_File{1};
  Material_name = Material_List_File{2};

  %Check thatthe requested material is in the list
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

  % Import water electron density and stopping powers (SP) at 100 MeV
  Water_SP = Import_SP_data(fullfile(MaterialsDirectory, 'Water', 'G4_Stop_Pow.dat') , E);
  Water_ElecDensity = Import_MatProperty(fullfile(MaterialsDirectory, 'Water', 'Material_Properties.dat'),'Electron_Density');

  %Read material properties
  ElecDensity = Import_MatProperty(fullfile(MaterialsDirectory, Name, 'Material_Properties.dat'),'Electron_Density');
  Density = Import_MatProperty(fullfile(MaterialsDirectory, Name, 'Material_Properties.dat'),'Density');

  % Compute SPR
  SP = Import_SP_data(fullfile(MaterialsDirectory, Name, 'G4_Stop_Pow.dat'), E);
  SPR = Density .* SP ./ Water_SP;
  RelElecDensity = Density * ElecDensity ./ Water_ElecDensity;


end

%===================================================
% Read the Geant 4 stopping power calibration curve
% and interpolate for therequested energies
% Author: K. Souris (Modifications: R. Labarbe)
%===================================================
function SP = Import_SP_data(FileName , E)

  fid=fopen(FileName,'r');
  if (fid < 0)
      error(['Unable to open file ' FileName])
  end

  data = fscanf(fid, '%f %f', [2 inf]);
  fclose(fid);
  SP = interp1(data(1,:), data(2,:), E, 'spline', 'extrap');

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
