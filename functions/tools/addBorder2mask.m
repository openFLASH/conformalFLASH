%% addBorder2mask
% Add a border of size |BlockBorder| around the mask of an aperture hole
%
%% Syntax
% |Apmask = addBorder2mask(Apmask , BlockBorder , pxlSize)|
%
%
%% Description
% |Apmask = addBorder2mask(Apmask , BlockBorder , pxlSize)| Description
%
%
%% Input arguments
% |Apmask| - _SCALAR MATRIX_ - Mask of the aperture hole |Apmask(x,y)=1| if the pixel is in the aperture hole
%
% |BlockBorder| -_SCALAR_- Width (mm) of the border to add around the hole
%
% |pxlSize| -_SCALAR VECTOR_- [x,y] Pixel size (mm)
%
%% Output arguments
%
% |Apmask| - _SCALAR MATRIX_ - Mask of the aperture block |Apmask(x,y)=1| if the pixel is in aperture block
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Apmask = addBorder2mask(Apmask , BlockBorder , pxlSize)

    AirHole = ~Apmask; %Set zero where the hole of the block is located

    sApmask = size(Apmask);
    Aidx = find(Apmask);
    [Xidx,Yidx] = ind2sub(sApmask , Aidx);

    minX = min(Xidx);
    maxX = max(Xidx);
    minY = min(Yidx);
    maxY = max(Yidx);
    Borderpxl =  round(BlockBorder ./ pxlSize);

    Apmask(minX-Borderpxl(1):maxX+Borderpxl(1) , minY-Borderpxl(2):maxY+Borderpxl(2)) = 1;

    Apmask = Apmask .* AirHole; %Make sure the hole still has the same size

end
