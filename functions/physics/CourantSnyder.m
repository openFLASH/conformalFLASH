%% CourantSnyder
% Compute the spot sigma, beam divergence and correlation coefficient of the PBS spot using the Courant Snyder equation
% The Courant Snyder function is eq (1) in [2]
%
% The same equation is valid for X and Y for the two bi-normal functions.
%
%
%% Syntax
% |[sxZ , stZ , rhoZ] = CourantSnyder(sigX0 , sigTh0 , rho0 , Z)|
%
%
%% Description
% |[sxZ , stZ , rhoZ] = CourantSnyder(sigX0 , sigTh0 , rho0 , Z)| Description
%
%
%% Input arguments
% |sigX0| -_SCALAR_- Spot sigma (mm) at the isocenter
%
% |sigTh0| -_SCALAR_- The beam divergence (radian) at isocenter
%
% |rho0| -_SCALAR_- The correlation coefficient between |sigX0| and |sigTh0| at isocenter
%
% |Z| -_SCALAR VECTOR_- The Z IEC gantry position (mm) of the plane in which the properties PBS spot are evaluated
%
%
%% Output arguments
%
% |sigX0| -_SCALAR VECTOR_- Spot sigma (mm) in a plane at Z(i) IEC gantry
%
% |sigTh0| -_SCALAR VECTOR_- The beam divergence (radian) in a plane at Z(i) IEC gantry
%
% |rhoZ| -_SCALAR VECTOR_- The correlation coefficient between |sigX0| and |sigTh0| in a plane at Z(i) IEC gantry
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)
%
% REFERENCE
% [1] http://www.openmcsquare.org/documentation_commissioning.html
% [2] Huang, S. et al. Validation and clinical implementation of an accurate Monte Carlo code for pencil beam scanning proton therapy. J Appl Clin Med Phys. 19, 558–572 (2018).

function [sxZ , stZ , rhoZ] = CourantSnyder(sigX0 , sigTh0 , rho0 , Z)


  sxZ = zeros(1,numel(Z));
  stZ  = zeros(1,numel(Z));
  rhoZ = zeros(1,numel(Z));

  %Vectorize the computation for all distances Z
  for idx = 1:numel(Z)
    [sxZ(idx) , stZ(idx) , rhoZ(idx)] = computeBeamParamAtZ(sigX0 , sigTh0 , rho0 , Z(idx));
  end

end

%----------------------------------
%compute the pram for one single Z
%----------------------------------
function [sxZ , stZ , rhoZ] =computeBeamParamAtZ(sigX0 , sigTh0 , rho0 , Z)

  Zm = [1 , -Z ; 0 , 1];
  Sm = [sigX0.^2 , rho0.*sigX0.*sigTh0 ; rho0.*sigX0 .*sigTh0 , sigTh0.^2];

  Sz = Zm * Sm * Zm'; %Equation (1) in [2]

  sxZ = sqrt(Sz(1,1));
  stZ = sigTh0; %divergence increases slightly with propagation in air due to muti coulomb scattering, [2] approximates it as a constant in air between the nozzle exit and the phantom surface
  rhoZ  = Sz(2,1) ./ (sxZ .* stZ); %positive for a defocusing beam and negative for a focusing beam.
end

%Developping the matrix multiplication gives the following eaquation for spot sigma in plane Z:
% sigX = sqrt(sigX0.^2 - 2 .* rho0 .* sigX0 .* sigTh0 .* z + sigTh0.^2 .* z.^2 );
% Eqaution from [1]
