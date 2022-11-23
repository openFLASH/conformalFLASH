%% fitDoubleGauss
% Fit a bi-normal gaussian to the 2D dose distribution
%
%% Syntax
% |[p,yf] = fitDoubleGauss(DoseMap , PtsG)|
%
%
%% Description
% |[p,yf] = fitDoubleGauss(DoseMap , PtsG)| Description
%
%
%% Input arguments
% |DoseMap| - _SCALAR MAP_ - |DoseMap(x,y)| Dose delivered at voxel (x,y)
%
% |PtsG| -_SCLAAR MATRIX_- |PtsG(i,:)=[x,y]| The coordinate (mm) of the i-th pixel in |DoseMap|
%
%% Output arguments
%
% |p| - _SCALAR VECTOR_ -  [Xmean , Ymean , sigmaX , sigmaY , correlation , DoseMAx]
%
% |yf| -_SCALAR MATRIX_- |yf(x,y)| Predicted dose at voxel (x,y)
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [p,yf] = fitDoubleGauss(DoseMap , PtsG)

  %Fit a bi-normal gaussian to the dose distribution
  options = optimset('Display','iter');
  dmax = max(DoseMap,[],'all');
  idx = find(DoseMap == dmax);
  sx = 3;
  sy=3;
  r=0;

  pin = [PtsG(idx(1),1) , PtsG(idx(1),2) , sx , sy , r , dmax .* (2 .* pi .* sx .* sy .* sqrt(1-r.^2))]; %Initial guess of the Gaussian shape
  p = fminsearch(@objective,double(pin),options,PtsG(:,1:2),DoseMap); %The optimum parameters
  yf = biNorm(PtsG(:,1:2) , p(6),  [p(1) , p(2)] , p(3) , p(4) , p(5));

end

%---------------------
% Gaussian fit
%---------------------
function gof = objective(p , x , y )
  A = p(6);
  r = p(5);
  sx = p(3);
  sy = p(4);
  mu = [p(1) , p(2)];

  gof = sum((y - biNorm(x , A,  mu , sx , sy , r)).^2);
end
