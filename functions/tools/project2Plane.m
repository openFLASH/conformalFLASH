%% project2Plane
% compute the magnification factor to project the coordinate from the isocenter planned
% to a plane located at some distance from isocenter
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
% |VDSA| -_SCALAR VECTOR_- [sadX , sadY]  mm source to axis distance for the 2 scanning magnets
%
% |DisoPlane| -_SCALAR_- distance (mm) from isocenter to projection plane. Positive if plane is between iosnceter and source. Negative if distance is beyond isocenter
%
%
%% Output arguments
%
% |mag| - _SCALAR VECTOR_ - Magnification factor [Mx,My] to project from the isocentre plane to the projection plane
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function mag = project2Plane(VDSA , DisoPlane)

  mag(1) = (VDSA(1) - DisoPlane) ./ VDSA(1); %magnification for projection from isocenter plan to hedgehog plane
  mag(2) = (VDSA(2) - DisoPlane) ./ VDSA(2);

end
