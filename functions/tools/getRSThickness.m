%% getRSThickness
% Compute the thickness of material required to reduce proton energy
% from |Ein| to |Eout|
%
%% Syntax
% |RSthick = getRSThickness(Ein , Eout , material)|
%
%% Description
% |RSthick = getRSThickness(Ein , Eout , material)| - Compute the thickness of material required to reduce proton energy
% from |Ein| to |Eout|
%
%% Input arguments
% |Ein| -_SCALAR_- Energy (MeV) of the incoming proton beam
% |Ein| -_SCALAR_- Energy (MeV) of the outgoing proton beam
% |material| - _STRING_ - Name of the material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%
%
%% Output arguments
% |RSthick| -_SCALAR_ Thickness (mm) of the mateiral required to reduce the proton energy

%% Contributors
% Authors : L. Hotoiu, R. Labarbe (open.reggui@gmail.com)

function  RSthick = getRSThickness(Ein , Eout , material)

    [range , Ep] = CSDArange(material); %Range (cm) in the RS material vs energy
    R = interp1(Ep , range , [Eout , Ein]); %Range (cm) in the material of the incoming particle
    RSthick = (R(2) - R(1)) .* 10;
end
