%% materialDescription
% Parameters for the Bragg Kelman equation [1] for different materials.
%
%
%% Syntax
% |material = materialDescription(name)|
%
%
%% Description
% |material = materialDescription(name)| Description
%
%
%% Input arguments
%
%  |name| - _STRING_ -  Description of the material
%
%% Output arguments
% |material| -_STRUCTURE_- Structure describing the properties of the material:
%    * |material.Name| -_STRING_- Name of the scatterer
%    * |material.Z| -_SCLAR VECTOR_- Z(i) Atomic number of the i-th element of scaterer |scatterIndex|
%    * |material.A| -_SCLAR VECTOR_- A(i) Atomic mass of the i-th element of  scaterer |scatterIndex|
%    * |material.f| -_SCLAR VECTOR_- f(i) Mass fraction of the i-th element of  scaterer |scatterIndex|
%    * |material.rho| -_SCALAR_- Density g/cm3
%    * |material.alpha| -_SCALAR_- Factor of the  stopping power vs energy curve.
%    * |material.p| -_SCALAR_- Exponent of the Bragg Kelman equation.
%    * |material.I| -_SCALAR_- Mean excitation energy computed by the Bragg additivity rule
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)
%
% REFERENCE
% [1] Bortfeld, T., & Introduction, I. (1997). An analytical approximation of the Bragg curve for therapeutic proton beams, 2024–2033.
% [2] Newhauser, W. D., & Zhang, R. (2015). The physics of proton therapy. Physics in Medicine and Biology, 60(8), R155–R209. https://doi.org/10.1088/0031-9155/60/8/R155
% [5] Gottschalk, B. (2009). On the scattering power of radiotherapy protons, 1–33. Retrieved from https://arxiv.org/pdf/0908.1413.pdf
% [13]	https://www.amazon.fr/ELEGOO-R%C3%A9sine-Impression-Liquides-Photopolym%C3%A8re/dp/B07Z96X1YX/ref=redir_mobile_desktop?ie=UTF8&aaxitk=hHAn4qWms-cA913ttsrZOw&hsa_cr_id=2103865140702&pd_rd_r=fc66fbf0-b5ba-49b3-a1da-d64e067da729&pd_rd_w=QeSt2&pd_rd_wg=bR4Qo&ref_=sbx_be_s_sparkle_mcd_asin_1_title
% [14]	https://www.elegoo.com/product/elegoo-3d-rapid-resin-lcd-uv-curing-resin-405nm-standard-photopolymer-resin-for-lcd-3d-printing-1000gram-grey/
% [16]	https://www.chemicalbook.com/ProductChemicalPropertiesCB32131152_EN.htm
% [18]	https://www.chemicalbook.com/CASEN_13048-33-4.htm
% [6] Newhauser, W. D., & Zhang, R. (2015). The physics of proton therapy. Physics in Medicine and Biology, 60(8), R155–R209. https://doi.org/10.1088/0031-9155/60/8/R155

function material = materialDescription(name)

physicsConstants;

switch name
case 'water'
    material.Name = 'water';
    f = [66.667,  33.333]; %H,O Molar fraction
    material.Z = [1,8];
    material.A = A(material.Z);
    material.f = f .* material.A ./ sum(f.*material.A); %Mass fraction
    material.rho=1;
    material.alpha = 2.633e-3;
    material.p=1.735;
    material.Xs = 46.88 ./ material.rho; %Scattering length (cm) for WATER Table 2 in [5]
    material.k = 0.012; % - _SCALAR_ -  material-dependent factor for computing the range stragling sigma. Value for water in [1]
    material.m = 0.935; % - _SCALAR_ -  experimentally determined exponent for computing the range stragling sigma. Value for water in [1]

  case 'PMMA'
      material.Name = 'PMMA';
      f = [53.333,  33.333 , 13.333]; %Molar fraction H 53.333, C 33.333, O 13.333
      material.Z = [1,6,8];
      material.A = A(material.Z);
      material.f = f .* material.A ./ sum(f.*material.A); %Mass fraction
      material.rho=1.185;
      material.alpha = 2.271e-3;
      material.p=1.735;

    case 'lexan'
        material.Name = 'lexan';
        f = [42.424, 48.485, 9.091]; %Molar fraction H 42.424, C 48.485, O 9.091
        material.Z = [1,6,8];
        material.A = A(material.Z);
        material.f = f .* material.A ./ sum(f.*material.A); %Mass fraction
        material.rho=1.2;
        material.alpha = 2.310e-3;
        material.p=1.735;
        material.Xs = 55.05 ./ material.rho; %Scattering length (cm) for WATER Table 2 in [5]

  case 'Pb'
    material.Name = 'Pb';
    material.Z = 82;
    material.A = 327.502;
    material.f = 1; %Mass fraction
    material.rho=11.322; %g/cm3
    material.alpha = 6.505e-4;
    material.p=1.676;

  case 'Cu'
    %Computed by fitting Bragg Kelman equation to the the CSDA range from https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html
    material.Name = 'Cu';
    material.Z = 29;
    material.A = 63.546;
    material.f = 1; %Mass fraction
    material.rho=8.96; %g/cm3
    material.alpha = 0.000510 ;
    material.p=1.707528 ;

    case 'carbon'
      material.Name = 'C';
      material.Z = Z(6);
      material.A = A(6);
      material.f = 1; %Mass fraction
      material.rho=2; %g/cm3

    case 'Al'
      material.Name = 'Al';
      material.Z = Z(13);
      material.A = A(13);
      material.f = 1; %Mass fraction
      material.rho=2.6989; %g/cm3
      material.alpha = 0.001543;
      material.p = 1.700463;
      %Computed by fitting Bragg Kelman equation to the stopping power equation (2) in [6]
      %The stopping power were estimated using function |stoppingpower|

  case 'brass'
    material.Name = 'brass';
    material.Z = [29 , 30 , 82];
    material.A = A(material.Z);
    material.f = [57 , 40.5 , 2.5] ./100;
    material.rho= 8.4 ; %g/cm3
    material.alpha = 0.000616 ;
    material.p= 1.688013  ;
    %Computed by fitting Bragg Kelman equation to the stopping power equation (2) in [6]
    %The stopping power were estimated using function |stoppingpower|

  case 'B4C' %Boron carbide
    material.Name = 'B4C';
    material.Z = [5 , 6];
    material.A = [A(5) , A(6)];
    f = [4 , 1]; %Molar fraction
    material.f = f .* material.A ./ sum(f.*material.A);
    material.rho= 2.52 ; %g/cm3
    material.alpha = 0.001420 ;
    material.p= 1.713110  ;
    %Computed by fitting Bragg Kelman equation to the stopping power equation (2) in [6]
    %The stopping power were estimated using function |stoppingpower|

  case 'UVcuredResin'
  % References [13,14,16,18]
  material.Name = 'UVcuredResin';
  %              H    C      N     O      F     Cl
  material.Z = [ 1   ,6     ,7    ,8     ,9    ,17];
  material.f = [6.81  66.37  2.43  14.97  3.29  6.14]; %Mass fraction
  material.f = material.f ./ sum(material.f); %normalise to 1
  material.A = A(material.Z);
  material.rho=1.195; %density

  case 'Tusk'
  %
  material.Name = 'Tusk';
  %              H    C      O      F     S     Sb
  material.Z = [ 1   ,6     ,8     ,9    ,16   ,51];
  material.f = [9.55 ,66.68 ,19.65 ,1.62 ,0.77 ,1.73]; %Mass fraction
  material.f = material.f ./ sum(material.f); %normalise to 1
  material.A = A(material.Z);
  material.rho=1.19; %density
  material.alpha = 0.002491 ;
  material.p= 1.716063   ;
  %Computed by fitting Bragg Kelman equation to the stopping power equation (2) in [6]
  %The stopping power were estimated using function |stoppingpower|

  otherwise
    error('Unknown material')
end

ZA = sum(material.f .* material.Z ./ material.A); %Effective Z
material.I = exp(sum( (material.f .* material.Z ./ material.A) .* log(Iz(material.Z)) ) ./ ZA);  % Mean excitation energy (J) by Bragg additivity. Unit: J (divide by eV to get unit in eV)


end
