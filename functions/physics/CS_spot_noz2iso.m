%% CS_spot_noz2iso
% Convert the Courant Snyder parameters from the nozzle exit plane to isocenter plane
% The MCsquare BDL specifies the Courant Snyder parameter in the plane of the nozzle exit
% The function CourantSnyder.m requires the parameters to be speciifed in the isocenter plane
%
%% Syntax
% |SpotDataISO = CS_spot_noz2iso(SpotDataNOZL)|
%
%
%% Description
% |SpotDataISO = CS_spot_noz2iso(SpotDataNOZL)| Description
%
%
%% Input arguments
% |SpotDataNOZL| - _STRUCTURE_ - The infroamtion about the PBS spot at the specified energy in the NOZZLE EXIT PLANE
%   * |SpotData.NominalEnergy| -_SCALAR VECTOR_- the nominal energy (in MeV) of the i-th energy.
%   * |SpotData.iso2Nozzle| -_SCALAR_- Distance (mm) from isocentre to nozzle exit. This is the plane in which the spot has a sigma|SpotData.SpotSizeNx|
%   * |SpotData.SpotSize1x(i)| -_SCALAR VECTOR_- the in air spot size (in mm) along the x direction at the nozzle exit (for the first 2D Gaussian).
%   * |SpotData.Divergence1x(i)| -_SCALAR VECTOR_- the spot divergence (in rad) along the x direction (for the first 2D Gaussian).
%   * |SpotData.Correlation1x(i)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the x direction (for the first 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |SpotData.SpotSize1y(i)| -_SCALAR VECTOR_- the in air spot size (in mm) along the y direction at the nozzle exit (for the first 2D Gaussian).
%   * |SpotData.Divergence1y(i)| -_SCALAR VECTOR_- the spot divergence (in rad) along the y direction (for the first 2D Gaussian).
%   * |SpotData.Correlation1y(i)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the y direction (for the first 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |SpotData.SpotSize2x(i)| -_SCALAR VECTOR_- the in air spot size (in mm) along the x direction at the nozzle exit (for the second 2D Gaussian).
%   * |SpotData.Divergence2x(i)| -_SCALAR VECTOR_- the spot divergence (in rad) along the x direction (for the second 2D Gaussian).
%   * |SpotData.Correlation2x(i)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the x direction (for the second 2D Gaussian). Its value must be between -0.99 and 0.99.
%   * |SpotData.SpotSize2y(i)| -_SCALAR VECTOR_- the in air spot size (in mm) along the y direction at the nozzle exit (for the second 2D Gaussian).
%   * |SpotData.Divergence2y(i)| -_SCALAR VECTOR_- the spot divergence (in rad) along the y direction (for the second 2D Gaussian).
%   * |SpotData.Correlation2y(i)| -_SCALAR VECTOR_- the correlation between the spot size and the divergence along the x direction (for the second 2D Gaussian). Its value must be between -0.99 and 0.99.
%
%
%% Output arguments
% |SpotDataISO| - _STRUCTURE_ - The information about the PBS spot at the specified energy in the ISOCENTER PLANE
%
%
% REFERENCE
% [1] http://www.openmcsquare.org/documentation_commissioning.html
% [2] 1. Huang, S. et al. Validation and clinical implementation of an accurate Monte Carlo code for pencil beam scanning proton therapy. J Appl Clin Med Phys. 19, 558–572 (2018).

%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function SpotDataISO = CS_spot_noz2iso(SpotDataNOZL)


  SpotDataISO = SpotDataNOZL;
  z = SpotDataISO.iso2Nozzle; %Z IEC gantry position of the nozzle exit

  label1 = {'SpotSize' , 'Divergence', 'Correlation'};
  label2 = {'x' , 'y'};


  for GausID = 1:2
    for AxeLabl  = 1:numel(label2)
          sigX0 = getfield(SpotDataNOZL , [label1{1} num2str(GausID) label2{AxeLabl}] ); %at isocenter plane
          sigTh0 = getfield(SpotDataNOZL , [label1{2} num2str(GausID) label2{AxeLabl}] ); %at isocenter plane
          rho0 = getfield(SpotDataNOZL , [label1{3} num2str(GausID) label2{AxeLabl}] ); %at isocenter plane

          [sxZ , stZ , rhoZ] = backPropagate(sigX0 , sigTh0 , rho0 , z);

          SpotDataISO = setfield(SpotDataISO , [label1{1} num2str(GausID) label2{AxeLabl}] , sxZ); %spot sigma at nozzle exit. Just apply Courant Snyder
          SpotDataISO = setfield(SpotDataISO , [label1{2} num2str(GausID) label2{AxeLabl}] , stZ); %divergence is the same at all planes. See [1]
          SpotDataISO = setfield(SpotDataISO , [label1{3} num2str(GausID) label2{AxeLabl}] , rhoZ); %Correlation at nozzle exit. See [1]
    end
  end

end

%Inverse the courant Snyder equation (1) in [2] in order to retrieve the value of the beam parameter
%in the isocenter plane
%
% INPUT
% |sigX0| -_SCALAR_- Spot sigma (mm) in a plane at Z IEC gantry
%
% |sigTh0| -_SCALAR_- The beam divergence (radian) in a plane at Z IEC gantry
%
% |rhoZ| -_SCALAR_- The correlation coefficient between |sigX0| and |sigTh0| in a plane at Z IEC gantry
%
% |z| -_SCALAR_- The Z IEC gantry position of the plane in which the beam parameter are specified
%
%OUTPUT
% |sigX0| -_SCALAR_- Spot sigma (mm) at the isocenter
%
% |sigTh0| -_SCALAR_- The beam divergence (radian) at isocenter
%
% |rho0| -_SCALAR_- The correlation coefficient between |sigX0| and |sigTh0| at isocenter

function [sigX0 , sigTh0 , rho0] = backPropagate(  sxZ , stZ , rhoZ, Z)

  Zm = [1 , -Z ; 0 , 1];
  Sz = [sxZ.^2 , rhoZ.*sxZ.*stZ ; rhoZ.*sxZ.*stZ , stZ.^2];

  Sm = inv(Zm) * Sz * inv(Zm') ; %Inversion of Equation (1) in [2]

  sigX0 = sqrt(Sm(1,1));
  sigTh0 = stZ; %divergence increases slightly with propagation in air due to MCS, [2] approximates it as a constant in air between the nozzle exit and the phantom surface
  rho0  = Sm(2,1) ./ (sigX0 .* sigTh0);

end
