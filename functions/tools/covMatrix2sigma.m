%% covMatrix2sigma
% Compute the spot sigma from the covariance matrix of the Gaussian lateral dose distribution
%
%% Syntax
% |spotSigma = covMatrix2sigma(p)|
%
%
%% Description
% |spotSigma = covMatrix2sigma(p)| Description
%
%
%% Input arguments
% |p| - _SCALAR VECTOR_ -  [~ , ~ , sigmaX , sigmaY , correlation , ~]
%
%
%% Output arguments
%
% |spotSigma| - _SCALAR_ - lateral sigma (mm) of the PBS spot along the main axis of the ellipse
%
% |sx| -_SCALAR_- Standard deviation (mm) along the X axis
%
% |sy| -_SCALAR_- Standard deviation alongthe X axis
%
% |r| -_SCALAR_- Correlation between X and Y
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [spotSigma , sx , sy , r ] = covMatrix2sigma(p)

  r = p(5);
  sx = p(3);
  sy = p(4);
  S = [sx.^2 r.*sx.*sy ; r.*sx.*sy sy.^2];
  [V, D] = eig(S); %Eigenvector and eigen value of the matrix
  ax = sqrt(diag(D));
  spotSigma = max(ax);


end
