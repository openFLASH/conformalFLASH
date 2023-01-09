%% readRS_HU2rho
% Read the text file containing the HU to density calibration of the Raystation.
% This is a text file with 2 columns:
%  * First column is the HU
%  * Second column is the density (g/cm3)
%
%% Syntax
% |modelHU2rho = readRS_HU2rho(filename)|
%
%
%% Description
% |modelHU2rho = readRS_HU2rho(filename)| Description
%
%
%% Input arguments
% |filename| - _STRING_ - File name (and path) of the text file containing the HU to density calibration exported from the Raystation
%
%
%% Output arguments
%
% |modelHU2rho| - _STRUCTURE_ -  conversion from HU to density
% * |modelHU2rho.Density| - _CELL MATRIX_ - Definition of the mass density of the tissues
% * ------|Density{i,1}| - _SCALAR_ - Minimum Hounsfield unit for the tissue |i|. The tissue has HU between Material{i,1}< HU <=Material{i+1,1}
% * ------|Density{i,2}| - _SCALAR_ - Mass density (g/cm3) of the tissue |i|.
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function modelHU2rho = readRS_HU2rho(filename)


  fid = fopen(filename,'r');

  if(fid==-1)
      error(['Cannot open file ',filename]);
  end

  tline = fgetl(fid);
  txt = [];
  while ischar(tline)
      txt = [txt,';',tline];
      tline = fgetl(fid);
  end
  txt = strrep(txt,'\t',' ');% removes tab
  for i=1:30
      txt = strrep(txt,'  ',' ');% removes multiple spaces
  end

  fclose(fid);

  %remove the first ';' and the [] character at the begining of the string
  txt = txt(3:end);
  txt = [ 'D = [' txt '];' ];
  eval(txt);

  modelHU2rho.Density = D;

end
