%% function
% Find the minimum and maximum coordinates of the brass block of aperture
% Find the radius of the circle circumbscribing the square aperture
% The dimensions are projected in the isocentre plane
%
%% Syntax
% |[minSize , maxSize , Rmax , BlockBorder] = findApertureBlockSize(BlockData)|
%
%
%% Description
% |[minSize , maxSize , Rmax , BlockBorder] = findApertureBlockSize(BlockData)| Description
%
%
%% Input arguments
% |ApertureBlockData|  -_SCALAR MATRIX_- |ApertureBlockData(i,:)=[x,y]|  Coordinates (mm) of the i-th point defining the contour of the aperture block projected onto the machine isocentric plane in the IEC BEAM LIMITING DEVICE coordinate system.
%
%
%% Output arguments
%
% |minSize| - _SCALAR VECTOR_ - [x,y] Extreme apex (mm) in the minus direction of the IEC gantry CS
%
% |maxSize| - _SCALAR VECTOR_ - [x,y] Extreme apex (mm) in the plus direction of the IEC gantry CS
%
% |Rmax| -_SCALAR_- Radius (mm) of the circle circumbscribing the suqare aperture block
%
% |BlockBorder| -_SCALAR_- Size (mm) of the border to add around the aperture hole
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [minSize , maxSize , Rmax , BlockBorder] = findApertureBlockSize(BlockData)

  BlockBorder = 50; %mm Width of brass to add outside of the aperture hole

  minSize = Inf;
  maxSize = -Inf;
  Rmax = 0;

  %Define the tranversal size of the range eshifter slabs
  for cntIdx = 1:numel(BlockData)
    minSize = min(min(BlockData{cntIdx},[],1)      , minSize);
    maxSize = max(max(BlockData{cntIdx},[],1)      , maxSize);
    Rmax    = max(max(sum(BlockData{cntIdx}.^2,2)) , Rmax   );
  end

  minSize = minSize - BlockBorder;
  maxSize = maxSize + BlockBorder;

  Rmax = sqrt(Rmax)  .* sqrt(2) + BlockBorder;

end
