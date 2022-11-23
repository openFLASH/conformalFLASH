%% getMoliereAngle
% Compute the characteristic multiple scattering angle and the reduced target thickness for the Moliere equation [3]
% The incident particle is a proton.
%
%% Syntax
% |[thtM , B] = getMoliereAngle(material,scatterIndex,r1)|
%
%
%% Description
% |[thtM , B] = getMoliereAngle(material,scatterIndex,r1)| Description
%
%
%% Input arguments
% |material| - _STRUCTURE_ - Description of the double scattering hardware
%  * |material| -_SCTRUCTURE_- description of the scaterrer
%    * |material.Xc| -_SCALAR_- (radian) Characteristic single scattering angle of the scatterer if the scatterer has uniform thickness
%    * |material.Xc| -_FUNCTION POINTER_-  @Xc(r1) gives the characteristic single scattering angle of the  scatterer at distance r2 (m) from optical axis if the scatterer has non uniform thickness
%    * |material.Z| -_SCLAR_- Atomic number of the element of scaterer
%
%  |E| -_SCALAR_- (MeV) Kinetic energy of the proton incident on the scatterer
%
% |r1| -_SCALAR VECTOR_- [OPTIONAL. Only needed when |material.Xc| is a _FUNCTION POINTER_]. The distance (m) at which the scattering angle is required
%
%% Output arguments
%
% |thtM| -_SCALAR_- (radian) Characteristic multiple angle scattering angle. Eq 7.11 in [3]. If |r1| is a vector, then it is also a vector of same length
%                   If the angular distribution is approximated by a Gaussian, this is the sigma of the Gaussian
%
% |B| -_SCALAR_- Reduced target thickness (dimensionless). Eq 7.10 in [3]. If |r1| is a vector, then it is also a vector of same length
%
% REFERENCE
% [3] Bethe. (1953). Moliere’s theory of multiple scatering. Physical Review, 89(6), 1256–1266.
% [4] https://gray.mgh.harvard.edu/attachments/article/212/pbs.pdf
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [thtM , B] = getMoliereAngle(material,E,r1)
  physicsConstants;
  z = 1; %Incident particle is a proton, charge = +1
  e2hc = 1./137.036; %Fine sturcture constant [4]
  % Note that ref [4] is not using SI units in the equation. See p10

  if (~isfield(material,'Xc2'))
    %We did not receive a Xc2. Let's compute it from Xc
    if (isa(material.Xc,'function_handle') )
      %This is a function. Return the characterisitc angle at each r1
      Xc2 = material.Xc(r1).^2; %single scattering angle
    else
      %This is a constant characteristic angle
      Xc2 = material.Xc.^2; %single scattering angle
    end
  else
    %We received the Xc2
    Xc2 = material.Xc2;
  end

  if (~isfield(material,'Xa2'))
    %We did not receive a Xa2. Let's compute it from Xa
    c1 = (e2hc.*z.*material.Z).^2; %eq 7.7 in [4]
    c2 = ((1./0.885) .* e2hc .* ((m_e .* c.^2) ./ MeV) .* material.Z.^(1/3)).^2; %eq 7.8 in [4]
    alph2 = c1 ./ T2beta2(E.*MeV); % Born parameter eq 7.6 in [4]

    X02 = c2 ./ (T2pc2(E .* MeV)./MeV.^2); %eq 7.5 in [4]. (pc)2 must be in MeV2
    Xa2 = X02 .* (1.13 + 3.76 .* alph2); % screening angle eq 7.5 in [4]
  else
    %We received the Xa2
    Xa2 = material.Xa2;
  end

  b = log(Xc2 ./ (1.167.*Xa2)); %eq 7.9 in [4]
  B = getB(b); % resolve eq 7.10 in [4] to get reduced target thickness

  thtM = sqrt(Xc2 .* B ./2); %Characteristic multiple angle scattering angle
  %It should be noted that standard Moliere theory [3] speaks instead of Xc * sqrt(B) which is the rms value of the space angle thethM.
  % Gottchalk [4] introduced the 1/sqrt(2) to put thetM on the same footing as the standard deviation in a Gaussian.

end





%===================
% Compute the reduced target thickness B
% Eq 7.10 in [2]
%===================
function B = getB(b)

  Bsc = 1:0.01:30;
  bsc = Bsc-log(Bsc);
  B = interp1(bsc,Bsc,b);
end
