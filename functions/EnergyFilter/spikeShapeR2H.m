%% spikeShapeR2H
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

function h = spikeShapeR2H(p,r)
  %h = p(1) + p(2) .* exp((r./p(4))-p(3)) + p(5) .* exp(((r./p(7))-p(6)).^2) ;
  h = p(1)   + (p(3) .* (r-p(2)) + p(4) .* (r-p(2)).^2 + p(5) .* (r-p(2)).^3) .* (r >= p(2)) ;

end
