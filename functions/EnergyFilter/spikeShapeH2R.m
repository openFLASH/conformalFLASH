%% spikeShapeH2R
% Function defining the shape of the edges of the spike
%
%% Syntax
% |res = help_header(im1,im2)|
%
%
%% Description
% |res = help_header(im1,im2)| Description
%
%
%% Input arguments
% |im1| - _STRING_ -  Name
%
%
%% Output arguments
%
% |res| - _STRUCTURE_ -  Description
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function a = spikeShapeH2R(p, h, apothem)
  %r = 0:0.1:apothem;
  r = p(2)+0.1:0.1:apothem;
  r= [0 , r];
  f = spikeShapeR2H(p,r);

figure(600)
hold on
plot(r,f,':b')

  a = interp1(f,r,h,'linear','extrap'); %Coordinate of centre of step of height h_step
  a(a>apothem) = apothem; %Clip the result in the range [0,apothem]
  a(a<0) = 0;
end
