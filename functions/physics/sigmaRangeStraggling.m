%% sigmaRangeStraggling
% Standard deviation on the distribution of range of monoenergetic proton beam
% due to range stranggling after traversing some material.
%The statistical process of proton energy loss in material worsens the beam distribution
%along the beam path. This effect, usually called straggling increase the beam spread in
%energy.
% The distribution in proton range, called range straggling, is a Gaussian with a width |sigmaR|.
% This standard deviation is computed using equation 2.16 in [1]
%
%% Syntax
% |res = help_header(im1,im2)|
%
%
%% Description
% |res = help_header(im1,im2)| Description
%
%
%% Input arguments
% |R| -_SCALAR_- Range (cm) of of the incoming proton beam
%
% |mater| - _STRING_ -  Description of the material
%
%
%% Output arguments
%
% |sigmaR| - _SCALAR_ - Standard deviation (cm) of the range after traversing the material
%
%% REFERENCE
% [1] Boon, S. N. (1998). Dosimetry and quality control of scanning proton beams. Groningen. Retrieved from https://www.rug.nl/research/portal/en/publications/dosimetry-and-quality-control-of-scanning-proton-beams(5539d204-4c09-47cd-9cc0-b7b56c51e082).html
% [2] ICRU. (1993). Ranges for Protons and Alpha Particles. International Commission on Radiation Units and Measurements. Report, 49(May).
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function sigmaR = sigmaRangeStraggling(R , mater)

physicsConstants;
material = materialDescription(mater);

p = material.p;
alpha = material.alpha;
Z = material.Z;
A = material.A;
w = material.f; %mass fraction

mc2 = m_e * (c^2); %Rest mass of electron (g)
K = 2 .*pi .* re.^2 .* mc2 .* Na; %Eq 2.2 in [1]. Unit J m2 / g
K = K ./ eV .* 1e-6 .* 1e4; %MeV cm2 /g
mc2 = m_e ./ MeV * (c^2); %electron mass (MeV)

Z_A_eff = sum(w .* Z ./ A); %Eq 2.22 in [2]

%equation 2.16 in [1]
sigmaR = sqrt(K .* 2 .* mc2 .* Z_A_eff  .* (p.^2 .* alpha.^(2./p) ./ (3-2./p)) .* R.^(3-2./p));

end
