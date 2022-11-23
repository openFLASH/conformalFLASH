%% getMaterialXc
% Compute the characteristic single scattering angle Xc and the Screening angle Xa for a compound material
% The particle beam is proton (z=1).
% Based on equation 7.20 and 7.21 of [4]
%
%% Syntax
% |material = getMaterialXc(material)|
%
%
%% Description
% |material = getMaterialXc(material)| Description
%
%
%% Input arguments
% |material| - _STRUCTURE_ - Description of the scatterer
%    * |material.Z| -_SCLAR VECTOR_- Z(i) Atomic number of the i-th element of scaterer
%    * |material.A| -_SCLAR VECTOR_- A(i) Atomic mass of the i-th element of  scaterer
%    * |material.f| -_SCLAR VECTOR_- f(i) Mass fraction of the i-th element of  scaterer
%    * |material.rho| -_SCALAR_- Density g/cm3
%    * |material.t| -_SCALAR_- Thickness (m) of material
%    * |material.alpha| -_SCALAR_- Factor of the  stopping power vs energy curve (Bragg Kelman equation). [5]
%    * |material.p| -_SCALAR_- Exponent of the Bragg Kelman equation.  [5]
%
%  |E| -_SCALAR_- (MeV) Kinetic energy of the proton incident on the scatterer
%
%
%% Output arguments
%
% |material| - _STRUCTURE_ - Description of the double scattering hardware
%    * |material.Xc2| -_SCALAR_- (radian^2) Square of the Characteristic single scattering angle of the scatterer
%    * |material.Xa2| -_SCALAR_- (radian^2) Square of the Screening angle of the scatterer
%
% REFERENCE
% [1] Damien Prieels (1995)
% [2] http://fisica.ciens.ucv.ve/~svincenz/TISPISGIMR.pdf
% [3] Bethe. (1953). Moliere’s theory of multiple scatering. Physical Review, 89(6), 1256–1266.
% [4] https://gray.mgh.harvard.edu/attachments/article/212/pbs.pdf
% [5] Bortfeld, T., & Introduction, I. (1997). An analytical approximation of the Bragg curve for therapeutic proton beams, 2024–2033.
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function material = getMaterialXc(material,E)

  persistent m_e  c MeV Na; %Make these variables persistent between function calls to reduce time in calls to |physicsConstants|

  if ~isscalar(m_e)
    %The presistent variables have not yet been defined. Load them from physicsConstants
    physicsConstants;
  end

  z = 1; %Incident particle is a proton, charge = +1
  e2hc = 1./137.036; %Fine sturcture constant [4]
  hc =197.327 .* 1e-15 .* 1e2; %MeV.fm See [4] p 10 ==> MeV.cm

  c1 = (e2hc.*z.*material.Z).^2; %eq 7.7 in [4]
  c2 = ((1./0.885) .* e2hc .* ((m_e .* c.^2) ./ MeV) .* material.Z.^(1/3)).^2; %eq 7.8 in [4]
  c3 = 4 .* pi .* Na .* e2hc.^2 .* hc.^2 .* z.^2 .* material.Z.^2 .* material.f ./ material.A;

  thickness = material.t.*100; %cm Convert from m to cm
  if(thickness < 0)
    error('material.t must be positive')
  end

  %Compute Xc Characteristic single scattering angle
  %----------------------------------------------
  Xc2 = integral(@(x) pv_m2(x,material.alpha,material.p,E,material.rho),0,thickness); %convert thickness in cm
  Xc2 = sum(c3) .* Xc2; %Eq 7.20 in [4]
  material.Xc2=Xc2;

  %Compute Xa Screening angle of the scatterer
  %----------------------------------------------
  lnXa2m = sum(c3 .* integral(@(x) lnXa_pv_m2(x,material.alpha,material.p,E,c1,c2,material.rho),0,thickness,'ArrayValued',true)) ./ Xc2;
  Xa2 = exp(lnXa2m);
  material.Xa2=Xa2;

end



% |x| : depth in the scatterer (cm)
% |alpha| - _SCALAR_ - (cm MeV^(-p)) Factor of the  range vs energy curve.
% |p| - _SCALAR_ - Exponent of the Bragg Kelman equation. No unit.
% |E| - _SCALAR_ - Kinetic energy (MeV) of the incoming particle
% |rho| -_SCALAR_- Material density g/cm3
function f = pv_m2(x, alpha,p, E,rho)
  mc2 = 938.27; %MeV for proton mass From [4]
  T = EvsDepth(x, alpha,p, E); %Kinetic energy (MeV) remaining at depth X in the material
  pv = T .* ( T + 2.* mc2) ./ (T + mc2); %eq 3.18 in [4]
  f = rho./ pv.^2; %'t' in Eq 7.20 in [4] is a distance in g/cm2. The function pv_m2 requires distance in cm, so multiply by density
end

% Compute the ln(Xa^2) at depth x inside material
%  * |x| : depth in the scatterer (cm)
%  * |alpha| - _SCALAR_ - (cm MeV^(-p)) Factor of the  range vs energy curve.
%  * |p| - _SCALAR_ - Exponent of the Bragg Kelman equation. No unit.
%  * |E| - _SCALAR_ - Kinetic energy (MeV) of the incoming particle
function f = lnXa2(x, alpha,p, E,c1,c2)
  physicsConstants;
  T = EvsDepth(x, alpha,p, E); %Energy (MeV) remaining at depth X in the material
  alph2 = c1 ./ T2beta2(T.*MeV); % Born materialeter eq 7.6 in [4]
  X02 = c2 ./ (T2pc2(T .* MeV)./MeV.^2); %eq 7.5 in [4]. (pc)2 must be in MeV2
  f = log(X02 .* (1.13 + 3.76 .* alph2)); % screening angle eq 7.5 in [4]

end

% compute integrant of eq 7.18 of [4]
%  * |x| : depth in the scatterer (cm)
%  * |alpha| - _SCALAR_ - (cm MeV^(-p)) Factor of the  range vs energy curve.
%  * |p| - _SCALAR_ - Exponent of the Bragg Kelman equation. No unit.
%  * |E| - _SCALAR_ - Kinetic energy (MeV) of the incoming particle
function f = lnXa_pv_m2(x, alpha,p, E,c1,c2,rho)
  ln_Xa2 = lnXa2(x, alpha,p, E,c1,c2);
  pv_m_2 = pv_m2(x, alpha,p, E,rho);
  f = pv_m_2 .* ln_Xa2;
end
