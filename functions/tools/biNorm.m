%% biNorm
% Bi normal Gaussian.
% Compute the value of the bi-normal distribution at point |Pts(i,:)=[x,y]| for a Gaussian centered at point |mu|
% with different sigma along the X and Y axis (i.e. elliptical Gaussian) and a correlation |r| between the 2 axis (tilt of the ellispe)
%
% |mu| can be the same value for all the points |Pts(i,:)| or a different value can be provided for each point.
%
%% Syntax
% |G = biNorm(Pts , A,  mu , sx , sy , r , ang)|
%
%
%% Description
% |G = biNorm(Pts , A,  mu , sx , sy , r , ang)| Description
%
%
%% Input arguments
% |Pts| -_SCALAR MATRIX_- |Pts(i,:)= [x,y]| Coordinates (mm) of the i-th point at which to compute the Gaussian
%
% |A| -_SCALAR_- Total area under the Gaussian
%
% |mu| -_SCALAR VECTOR_- Option 1: |numel(mu)=2| : |mu = [x,y]| Coordinate (mm) of the centre of the Gaussian
%                        Option 2: |numel(mu)=size(Pts,1)| : |mu(i,:) = [x,y]| is the position of the centre of the Gaussian when computing the Gaussian value for the point  |Pts(i,:)|
%
% |sx| -_SCALAR_- Standard deviation (mm) along the X axis
%
% |sy| -_SCALAR_- Standard deviation (mm) along the Y axis
%
% |r| -_SCALAR_- Correlation between X and Y
%
% |ang| -_SCALAR_- Rotation angle (rad) around Ziec axis to go from IEC gantry to the ellipses axes
%
%% Output arguments
%
% |G| - _SCALAR VECTOR_ - |G(i)| Value of the bi-normal Gaussian at the i-th coordinate
%
%
%% REFERENCE
% [1] https://en.wikipedia.org/wiki/Multivariate_normal_distribution
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function G = biNorm(Pts , A,  mu , sx , sy , r , ang)

  if nargin < 7
    ang = 0;
  end

  NbPts = size(Pts,1);
  if numel(mu)==2
    X = Pts - repmat(mu,NbPts,1);
  else
    X = Pts - mu;
  end

  R = rot(-rad2deg(ang) , [0,0,0]); %Rotation matrix around Z
  M = R(1:2,1:2);
  X = M * X';
  X = X';

  Xsx = X(:,1) ./ sx;
  Ysy = X(:,2) ./ sy;
  u = (Xsx.^2 - 2.* r .* Xsx .* Ysy + Ysy.^2) ./ (2.*(1-r.^2));
  G = A .* exp(-u) ./ (2 .* pi .* sx .* sy .* sqrt(1-r.^2));
end
