%% moliere
% Relative number of particle scattered between angle |theta(i)| and |theta(i)+ d_theta|.
% The number of particle is already integrated over d(phi), so that it is the number of particles is given in a shell of thickness dr at distance r.
% Computed according to Moliere theory [2].
% The function is normalised so that integral_0^inf d(theta) integral_0^(2pi)d(phi) f = 1
% At the first call, the function tabulate the F0,F1 and F2 function so that the later calls are interpolated and hence faster
%
%% Syntax
% |[f,F0,F1,F2]=moliere(theta,B,force)|
%
%
%% Description
% |[f,F0,F1,F2]=moliere(theta,B,force)| Description
%
%
%% Input arguments
% |theta| - _SCALAR VECTOR_ - (rad) Angle variable in Moliere's theory
%
% |thetaM| -_SCALAR_- (radian) Characteristic multiple angle scattering angle. Eq 7.11 in [1]. If |r1| is a vector, then it is also a vector of same length
%
% |B| -_SCALAR_- Reduced target thickness (dimensionless). Eq 7.10 in [1]. If |r1| is a vector, then it is also a vector of same length
%
% |force| -_BOOLEAN_- [OPTIONAL, default = false]. Force the computation using Bessel function. Do not use interpolations
%
%% Output arguments
%
% |f| - _SCALAR VECTOR_ -   |f(i)| is the relative number of particle scattered between angle |theta(i)| and |theta(i)+ d_theta|.
%                           This is f(theta) in Moliere's theory
%
% |F0| - _SCALAR VECTOR_ -  First term in Moliere formula: |F(i)|  for angle |theta(i)|
%
% |F1| - _SCALAR VECTOR_ -  Second term in Moliere formula: |F(i)|  for angle |theta(i)|
%
% |F2| - _SCALAR VECTOR_ -  Third term in Moliere formula: |F(i)|  for angle |theta(i)|

%
%% REFERENCE
% [1] https://gray.mgh.harvard.edu/attachments/article/212/pbs.pdf
% [2] Bethe. (1953). Moliere’s theory of multiple scatering. Physical Review, 89(6), 1256–1266.
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)



function [f,F0,F1,F2]=moliere(theta,thetaM,B,force)

  persistent moliereData; %Store the tabulated data in a persistent variable. Later function calls are faster because interpolated from this table
  if nargin <4
    force = false;
  end

  if force
    %do not use the fast interpolation
    %recompute the integral of Bessel function
    [f,F0,F1,F2] = computeMoliere(theta./thetaM,thetaM,B);
    return;
  end

  if(~isfield(moliereData,'theta'))
      fprintf('Computing global variable ....\n')
      %No pre-computation available. Tabulate the function so that we can interpolate in the future
      moliereData = computedigitalMoliere;
      fprintf('done \n')
    end

  %Interpolate the value of Moliere from the tabulated function
  F0 = interp1(moliereData.theta,moliereData.F0,theta ./ thetaM,[],0); %default interpolator. Extrapolate with 0
  F1 = interp1(moliereData.theta,moliereData.F1,theta ./ thetaM,[],0);
  F2 = interp1(moliereData.theta,moliereData.F2,theta ./ thetaM,[],0);

  f = addMoliere(F0,F1,F2,thetaM,B); %Multiply by theta to get the probability density

end

%========================

function f = integrant(y,theta,n)
    f = y .* besselj(0,theta .* y) .* exp(-y.^2 ./ 4) .* ((y.^2 ./4) .* log(y.^2 ./ 4)).^n;
    % Note: in [1], there is a missing - sign in the exponential
end

%========================
% Make a tabulation of the function so that in the future
% we can interpolate the funciton in order to make the computation fast
%========================
function moliereData = computedigitalMoliere()
    thetaM = 38.038; %From table 7.1 in [1] for Lead atthickness =  0.1 * range
    B = 12.877;
    %theta = 1e-5:1e-5:10; %Reduced angle
    %The tabulated points are the same as in table II of [2]
    theta = 0:0.2:4;
    theta = [theta , 4.5:0.5:6];
    theta = [theta , 7:10];
    moliereData.theta=theta;
    [~,moliereData.F0,moliereData.F1,moliereData.F2] = computeMoliere(moliereData.theta,thetaM,B);

end


%====================
% The computation of the integral of Bessel functions
%====================
function [f,F0,F1,F2] = computeMoliere(thetaR,thetaM,B)
  fprintf('Computing Moliere \n')
    %See table II in [2] for tabulated value of these functions
    F0 = 2 .* exp(-thetaR.^2); %eq 7.15 in [1]
    F1 =        integral(@(x) integrant(x, thetaR,1),0,inf,'ArrayValued',true); %1! = 1
    F2 = 0.5 .* integral(@(x) integrant(x, thetaR,2),0,inf,'ArrayValued',true); % 2! = 2
    f = addMoliere(F0,F1,F2,thetaM,B);
end

%=====================================================
function f = addMoliere(F0,F1,F2,thetaM,B)
  %f = (F0 + F1 ./ B + F2 ./ B.^2) ./ (2.*pi.*thetaM.^2);
  f = (F0 + F1 ./ B + F2 ./ B.^2) ./ (sqrt(pi).*thetaM); % Normalised for d(theta)
  %TODO Check normalisation
  % There is no factor (1/2) because we do not use sqrt(2) in the defintion of theta_m
  % The normalisation factor of [1] does not work for the exp(-t^2) function. This is a 1D Gaussian, not 2D
  % The function is normalised so that integral_0^inf d(thetah) integral_0^(2pi)d(phi) f = 1
end
