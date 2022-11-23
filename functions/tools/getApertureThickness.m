%% getApertureThickness
%Compute the thickness of the brass block in order to stop a beam with the given energy
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
% |E| - _SCALAR_ - Energy (MeV) of the beam to be stopped by the aperture block
%
%
%% Output arguments
%
% |BlockThickness| - _SCALAR_ - Thickness (mm) of the aperture block
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function BlockThickness = getApertureThickness(E)

  material = materialDescription('Cu');
  %BlockThickness = energy2range(E, material.alpha,material.p) .*10 + 10; %Maximum range of the incoming particle (mm) + safety margin
  BlockThickness = 60; %mm NB: The block thickness will be always be set to the maximum thickness (60mm) so that there are no error in manufacturing

end
