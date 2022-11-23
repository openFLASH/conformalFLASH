%% getSAmodelParam
% Get the polynomial coefficient of the semi-analytical model
% The coefficient are used to predict the spot lateral sigma and range straglging
% when the proton beam goes through a specified tickness of CEF
%
%% Syntax
% |[pS , pE ] = getSAmodelParam()|
%
%
%% Description
% |[pS , pE ] = getSAmodelParam()| Description
%
%
%% Input arguments
% |im1| - _STRING_ -  Name
%
%
%% Output arguments
%
% |pS| - _SCALAR VECTOR_ - Polynomial coefficients to predict the sigma of the lateral profile of the spot
%
% |pE| - _SCALAR VECTOR_ - Polynomial coefficients to predict the sigma of the range straggling
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [pS , pE ] = getSAmodelParam()

%Polynomial coefficient optimised on 3 spikes
 % pS = [0.332711 , -16.877937 , 2356.661926 ] ./ 1000; %Polynomial coefficient for spot sigma vs thickness
 % pE = [-29.556142 , -102.338602  ,5.158769 ] ./ 1000; %Polynomial coefficient for range vs thickness

%Polynomial coefficient optimised on 6 spikes
% pS = [0.454753 , -18.181926  , 2697.250498 ] ./ 1000; %Polynomial coefficient for spot sigma vs thickness
% pE = [5.011806 , -147.685855 , 4.522610 ] ./ 1000; %Polynomial coefficient for range vs thickness

%Polynomial coefficient optimised on 6 spikes
%Range shifter is B4C. Range shifter position defined by downstream surface
pS = [0.353219 , -20.121663  , 2917.469822 ] ./ 1000; %Polynomial coefficient for spot sigma vs thickness
pE = [3.997292 , -118.637974 , 5.279108 ] ./ 1000; %Polynomial coefficient for spot sigma vs thickness
end
