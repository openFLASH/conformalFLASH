%% fitEllipse
% Compute the ellipse that has the same second-moments as the region containing the PBS spots.
%
%% Syntax
% |stats = fitEllipse(spot)|
%
%
%% Description
% |stats = fitEllipse(spot)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y]
%
%
%% Output arguments
%
% |stats| - _STRUCTURE_ -  Description
%   * |stats.MinorAxisLength| -_SCALAR_- Length (in mm) of the minor axis of the ellipse that has the same normalized second central moments as the region.
%   * |stats.MajorAxisLength| -_SCALAR_- 	Length (in mm) of the major axis of the ellipse that has the same normalized second central moments as the region
%   * |stats.Orientation| -_SCALAR_- 	Angle (degree) between the x-axis and the major axis of the ellipse that has the same second-moments as the region
%   * |stats.stats.Centroid| -_SCALAR_- Center of mass of the region
%
%% REFERENCES
% [1] https://stackoverflow.com/questions/1532168/what-are-the-second-moments-of-a-region
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)


function stats = fitEllipse(spot)

mu = mean(spot,1);
X_minus_mu = spot - repmat(mu, size(spot,1), 1);
Sigma = (X_minus_mu' * X_minus_mu) / size(spot,1); %Covariance matrix = second order moments

[V, D] = eig(Sigma); %Eigenvector and eigen value of the matrix
ax = sqrt(diag(D)) .* 4; %Half width of ellipse is 2 * sigma. The full length is 4 * sigma

stats.MinorAxisLength = min(ax);
[stats.MajorAxisLength , idxMaj]= max(ax);
Vec = V(idxMaj , :);
stats.Orientation =  rad2deg(atan(Vec(2) ./ Vec(1)));
stats.Centroid = mu;

% t = -50:50;
% x = mu(1) + t;
% y = t .* Vec(2) ./ Vec(1) + mu(2);
% plot(x,y,'-r')
% pause


end
