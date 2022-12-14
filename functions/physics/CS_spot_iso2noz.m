%% CS_spot_iso2noz
% Convert the Courant Snyder parameters from isocenter plane to the nozzle exit plane
% The MCsquare BDL specifies the Courant Snyder parameter in the plane of the nozzle exit
% The function CourantSnyder.m requires the parameters to be speciifed in the isocenter plane
%
%% Syntax
% |SpotDataNOZL = CS_spot_iso2noz(SpotDataISO)|
%
%
%% Description
% |SpotDataNOZL = CS_spot_iso2noz(SpotDataISO)| Description
%
%
%% Input arguments
% |SpotDataISO| - _STRUCTURE_ - The infroamtion about the PBS spot at the specified energy in the ISOCENTER PLANE
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
% |SpotDataNOZL| - _STRUCTURE_ - The information about the PBS spot at the specified energy in the NOZZLE EXIT PLANE
%
%
% REFERENCE
% [1] http://www.openmcsquare.org/documentation_commissioning.html
% [2] 1. Huang, S. et al. Validation and clinical implementation of an accurate Monte Carlo code for pencil beam scanning proton therapy. J Appl Clin Med Phys. 19, 558–572 (2018).

%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function SpotDataNOZL = CS_spot_iso2noz(SpotDataISO)


  SpotDataNOZL = SpotDataISO;
  z = SpotDataISO.iso2Nozzle; %Z IEC gantry position of the nozzle exit

  label1 = {'SpotSize' , 'Divergence', 'Correlation'};
  label2 = {'x' , 'y'};


  for GausID = 1:2
    for AxeLabl  = 1:numel(label2)
          sigX0 = getfield(SpotDataISO , [label1{1} num2str(GausID) label2{AxeLabl}] ); %at isocenter plane
          sigTh0 = getfield(SpotDataISO , [label1{2} num2str(GausID) label2{AxeLabl}] ); %at isocenter plane
          rho0 = getfield(SpotDataISO , [label1{3} num2str(GausID) label2{AxeLabl}] ); %at isocenter plane

          [sxZ , stZ , rhoZ] = CourantSnyder(sigX0 , sigTh0 , rho0 , z); %Apply eq (1) in [2]

          SpotDataNOZL = setfield(SpotDataNOZL , [label1{1} num2str(GausID) label2{AxeLabl}] , sxZ); %spot sigma at nozzle exit.
          SpotDataNOZL = setfield(SpotDataNOZL , [label1{2} num2str(GausID) label2{AxeLabl}] , stZ); %divergence
          SpotDataNOZL = setfield(SpotDataNOZL , [label1{3} num2str(GausID) label2{AxeLabl}] , rhoZ); %Correlation at nozzle exit
    end
  end

end
