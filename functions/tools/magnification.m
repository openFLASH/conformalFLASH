%% magnification
% Compute the magnification factor when projecting from the isocentre plane to a plane located at distance |DisoPlane|
% from isocentre.
% The projection is done from the proton source, defined in the beam data library.
% Take into account the anamorphic projection due to different SADx and SADy.
%
%% Syntax
% |mag = magnification(DisoPlane , BDL_file)|
%
%
%% Description
% |mag = magnification(DisoPlane , BDL_file)| Description
%
%
%% Input arguments
%
% |DisoPlane| -_SCALAR_- Distance (mm) from isocentre to the projection plane. Positive if projeciton plane is toawards the source. Negative otherwise
%
% |BDL_file| -_STRING_- File name of the beam data library. Used to retreive the SADx and SADy
%
%% Output arguments
%
% |mag| - _SCALAR VECTOR_ - Magnification factor [Mx,My] to project from the isocentre plane to the projection plane
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function mag = magnification(DisoPlane , BDL_file)

  [~ , sad_X , sad_Y] = get_sad(BDL_file);
  mag = project2Plane([sad_X , sad_Y] , DisoPlane);

end
