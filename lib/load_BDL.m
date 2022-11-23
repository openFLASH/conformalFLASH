%% load_BDL
% Load the beam data library from the MCsquare data folder.
% The BDL provides the PBS spot parameters at the exit of the nozzle for the Courant-Snyder formula.
% Code of the function from reference [1].
% Information about the Beam Data Library from [2].
%
%% Syntax
% |BDL = load_BDL(BDL_file)|
%
%
%% Description
% |BDL = load_BDL(BDL_file)| Description
%
%
%% Input arguments
% |BDL_file| -_STRING_- path to the BDL file
%
%
%% Output arguments
%
% |BDL| - _STRUCTURE_ - The infroamtion from the beam data library
%   * |BDL.NominalEnergy| -_SCALAR VECTOR_- the nominal energy (in MeV) as specified in the treatment plan.
%   * |BDL.MeanEnergy(E)| -_SCALAR VECTOR_- the actual mean energy (in MeV) that should be simulated to reproduce the proton range measured in water.
%   * |BDL.EnergySpread(E)| -_SCALAR VECTOR_- the standard deviation of the Gaussian distribution used to model the energy spectrum at the exit of the nozzle (in % of the nominal energy).
%   * |BDL.ProtonsMU(E)| -_SCALAR VECTOR_- the number of protons delivered per MU at this specific nominal energy.
%   * |BDL.Weight1(E)| -_SCALAR VECTOR_- he weight of the first Gaussian for the double Gaussian optical model. The sum of Weight1 and Weight2 should be 1.0.
%   * |BDL.SpotSize1x(E)| -_SCALAR VECTOR_- the in air spot size (in mm) along the x direction at the nozzle exit (for the first 2D Gaussian).
%   * |BDL.Divergence1x(E)| -_SCALAR VECTOR_- the spot divergence (in rad) along the x direction (for the first 2D Gaussian).
%   * |BDL.Correlation1x(E)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the x direction (for the first 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |BDL.SpotSize1y(E)| -_SCALAR VECTOR_- the in air spot size (in mm) along the y direction at the nozzle exit (for the first 2D Gaussian).
%   * |BDL.Divergence1y(E)| -_SCALAR VECTOR_- the spot divergence (in rad) along the y direction (for the first 2D Gaussian).
%   * |BDL.Correlation1y(E)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the y direction (for the first 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |BDL.Weight2(E)| -_SCALAR VECTOR_- the weight of the second Gaussian for the double Gaussian optical model. The sum of Weight1 and Weight2 should be 1.0.
%   * |BDL.SpotSize2x(E)| -_SCALAR VECTOR_- the in air spot size (in mm) along the x direction at the nozzle exit (for the second 2D Gaussian).
%   * |BDL.Divergence2x(E)| -_SCALAR VECTOR_- the spot divergence (in rad) along the x direction (for the second 2D Gaussian).
%   * |BDL.Correlation2x(E)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the x direction (for the second 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |BDL.SpotSize2y(E)| -_SCALAR VECTOR_- the in air spot size (in mm) along the y direction at the nozzle exit (for the second 2D Gaussian).
%   * |BDL.Divergence2y(E)| -_SCALAR VECTOR_- the spot divergence (in rad) along the y direction (for the second 2D Gaussian).
%   * |BDL.Correlation2y(E)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the x direction (for the second 2D Gaussian). Its value must be between -0.99 and 0.99.
%
%% REFERENCE
% [1] https://gitlab.com/openmcsquare/commissioning/-/blob/master/Tools/load_BDL.m
% [2] http://www.openmcsquare.org/documentation_commissioning.html
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function BDL = load_BDL(BDL_file)

fid = fopen(BDL_file);
tags = {};

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if(not(isempty(strfind(tline,'NominalEnergy'))))
        tags = get_values(tline,0);
        for i=1:length(tags)
            BDL.(tags{i}) = [];
        end
    elseif(not(isempty(tags)))
        values = get_values(tline,1);
        for i=1:length(tags)
            BDL.(tags{i})(end+1) = values(i);
        end
    end
end

fclose(fid);

end

function values = get_values(tline,num)
tline = regexprep(tline,'\t',' ');
for i=1:20
    tline = strrep(tline,'  ',' ');
end
temp = strsplit(tline,' ');
if(num)
    for j=1:length(temp)
        values(j) = str2double(temp{j});
    end
else
    values = temp;
end
end
