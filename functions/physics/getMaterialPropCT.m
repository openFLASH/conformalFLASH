%% getMaterialPropCT
% Retrieve the stopping power, relative electron density, density and Hounsfield unit for a given material.
% The information is retrieved from the MCsquare database for the specified CT scanner
%
%% Syntax
% |[HU, Density, SPR, SP, RelElecDensity] =  getMaterialPropCT(material)|
%
%
%% Description
% |[HU, Density, SPR, SP, RelElecDensity] =  getMaterialPropCT(material)| Description
%
%
%% Input arguments
% |material| - _STRING_ - Name of the material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%
% |ScannerDirectory| - _STRING_ - Name of the folder containing the definition of the CT scanner properties in MCsquare in folder "plugins\openMCsquare\lib\Scanners"
%
%% Output arguments
%
% |HU| -_SCALAR_- Hounsfield unit of the material in the scanner
%
% |Density| -_SCALAR_- Density (g/cm3) of the material
%
% |SPR| -_SCALAR_- Relative stopping power (relative to water) at 100 MeV SPR(i) = Density(i) * SP(i) / Water_SP;
%
% |SP| -_SCALAR_- Mass stopping powers  (MeV cm2/g) at 100 MeV
%
% |RelElecDensity| -_SCALAR_- Electron density relative to water
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [HU, Density, SPR, SP, RelElecDensity] =  getMaterialPropCT(material, ScannerDirectory)


  HUmaterialPath = fullfile(ScannerDirectory,'HU_Material_Conversion.txt');

  [HUm, materialID , descr] = MC2_import_scanner_file(HUmaterialPath); %Read link between HU and material ID
  [HUct, Densitiesct, SPRct, SPct, RelElecDensityct] = Compute_SPR_data(ScannerDirectory); %Read link between HU and physical properties

  checkList = strfind(descr,material);
  found =0;
  for idx = 1:numel(checkList)
    if ~isempty(checkList{idx})
      found = idx;
      break; %We have a match
    end
  end

  if found
    %We have found the requested material
    HU = HUm(found);
    Density = interp1(HUct, Densitiesct , HUm(found), 'linear','extrap') ;
    SPR = interp1(HUct, SPRct , HUm(found), 'linear','extrap') ;
    SP = interp1(HUct, SPct, HUm(found), 'linear','extrap') ;
    RelElecDensity = interp1(HUct, RelElecDensityct , HUm(found) , 'linear','extrap') ;
  else
    %The material is not in the list
    material
    ScannerDirectory
    error('The requested material is not available in ScannerDirectory')
  end

end
