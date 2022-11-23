%% getSpotSigmaMoliere
% Get the sigma of the lateral spread of the beam due to scattering, as described by the Moliere theory [3,4]
% The sigma is computed at theexit of the scatterer of thickness |l|
% and at a |distance| downstream from the scatterer
%
%% Syntax
% |xe = getSpotSigmaMoliere(T,l,distance, matl)|
%
%
%% Description
% |xe = getSpotSigmaMoliere(T,l,distance, matl)| Description
%
%
%% Input arguments
% |T| -_SCALAR_- energy (MeV)
%
% |l| -_SCALAR VECTOR_- |l(i)| Thickness (mm) of the i-th object
%
% |distance| -_SCALAR_- Distance (mm) from the scatterer at which the sigma of the spot is computed
%
% |matl| - _STRING_ -  Description of the material. See |materialDescription| for more information
%
%
%% Output arguments
%
% |xe(i)| -_SCALAR VECTOR_- sigma (mm) of the lateral spread of the beam at |distance| at the top of the step of height |l(i)|
%
% |x0(i)| - _SCALAR_ - sigma (mm) of the lateral spread of the beam at the top of the step of height |l(i)|
%
% |theta(i)| -_SCALAR_- (radian) Characteristic multiple angle scattering angle at the top of the step of height |l(i)|. Eq 7.11 in [3].
%
% REFERENCE
% [3] Bethe. (1953). Moliere’s theory of multiple scatering. Physical Review, 89(6), 1256–1266.
% [4] https://gray.mgh.harvard.edu/attachments/article/212/pbs.pdf
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [xe , x0 , thtM] = getSpotSigmaMoliere(T,l,distance, matl)

    material = materialDescription(matl);
    x0 = zeros(1,numel(l));
    xe = zeros(1,numel(l));

    for i = 1:numel(l)
      material.t = l(i) / 1000; % -_SCALAR_- Thickness (m) of material.
      material = getMaterialXc(material,T);

      %----calculate x0
      thtM(i)  = getMoliereAngle(material,T); %scattering angle in radian
      x0(i) = thtM(i) .* l(i); %l is in mm so X0 is also in mm
      xe(i) = sqrt((thtM(i) .* distance).^2 + x0(i).^2); %distance in mm
    end
end
