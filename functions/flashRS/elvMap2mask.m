%% elvMap2mask
% Convert the elevation map of the CEM form a linear vector of number into a 3D maxk and into a 2D elevation map.
% REinterpolate the elevtation map into a high resolution grid with pixel size |[pixelSize,pixelSize,pixelSize]|
%
%% Syntax
% |[CEM3Dmask , CEMThicknessData] = elvMap2mask(ElvMap , nrPixelsX , nrPixelsY , pixelSize)|
%
%
%% Description
% |[CEM3Dmask , CEMThicknessData] = elvMap2mask(ElvMap , nrPixelsX , nrPixelsY , pixelSize)| Description
%
%
%% Input arguments
% |ElvMap| - _SCALAR VECTOR_ - |ElvMap(i)| Height (mm) of the i-th pixel on the elevation map of the CEM
%
% |nrPixelsX| - _SCALAR_ - Number of pixels along the X-IEC axis in |ElvMap|
%
% |nrPixelsY| - _SCALAR_ - Number of pixels along the Y-IEC axis in |ElvMap|
%
% |pixelSize| -_SCALAR VECTOR_- [x,y,z] size (mm) of the pixels in |ElvMap|
%
%
%% Output arguments
%
% |CEM3Dmask| -_SCALAR MATRIX_- |CEM3Dmask{b}(x,y,z)| 3D mask of the CEF. |CEM3Dmask(x,y,z)=1| if the voxel at location (x,y,z) belongs to the CEF.
%                                     The index z increases in the direction of increasing Zg. x and y indices are aligned with the IEC gantry CS
%
% |CEMThicknessData| -_SCALAR MATRIX_- |CompensatorThicknessData(x,y)| Thickness (mm) of the CEF pixel at position (x;y) in the IEC beam Limiting device CS
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [CEM3Dmask , CEMThicknessData ] = elvMap2mask(ElvMap , nrPixelsX , nrPixelsY , pixelSize)

    % DICOM order the pixels row by row. The row are paralell to the X axis (IEC gantry).
    % The first row is at +Y. The last row is at -Y
    % Within one row, the pixels are ordered from -x to +X
    %
    %                         ^ +Y
    % ---1----------------->   |
    % ---2----------------->   | -----> +X
    % ---3----------------->   |
    % ---4----------------->   |
    %
    % The elemnets are therefore ordered [A(:,N) , A(:,N-1) , A(:,N-2), ....]
    % where the first index is X (row) and the second index is Y (column)
    %
    % On the other hand, in Matlab, "reshape" assumes that the elements are ordered [A(:,1) , A(:,2) , ....]
    % We therefore need to flip the second dimension of the matrix

    ElvMap = reshape(ElvMap , nrPixelsX, nrPixelsY);
    ElvMap = flip(ElvMap , 2); %In DICOM, the row start at +Y and run towrds -Y. We need to flip the second dimension

    CEMThicknessData = ElvMap; % mm Elevation map
    maxEl = max(ElvMap,[],'all'); %Maximum height of the elevation map

    %Create an interpolation grid at the resolution of pixelSize
    X = 1:nrPixelsX;
    Y = 1:nrPixelsY;
    X = X .* pixelSize(1);
    Y = Y .* pixelSize(2);

    [~ , ~ , VertDist] = meshgrid( Y , X , 1:pixelSize(3):maxEl); %meshgrid inversion the 1st and second index
    ElvMap3D = repmat(ElvMap , 1 , 1 , size(VertDist,3)); %Create a 3D map of the vertical distances
    CEM3Dmask = ((VertDist <= ElvMap3D) .* (ElvMap3D > 0)); %Convert the 3D elevationation map into a 3D binary mask

end

%Ordering of the element between matrices and vector
% A = [1 2 3 ; 10 20 30 ; 100 200 300]
% A =
%
%      1     2     3
%     10    20    30
%    100   200   300
%
% >> B = A(:)
% B =
%
%      1
%     10
%    100
%      2
%     20
%    200
%      3
%     30
%    300
%
% >> C = reshape(B, 3 , 3)
% C =
%
%      1     2     3
%     10    20    30
%    100   200   300
