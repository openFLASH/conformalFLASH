%% interpolateCT
% Interpolate the CT on the new grid from the CT on the old grid
% Make sure the the indices of the new CT are oriented along the same axes as the old CT
% The coordinate of (Xint, Yint, Zint) of the voxels of the interpolated CT scan are expressed in the smae coordinate than the original CT scan
%
%% Syntax
% |CTintrp = interpolateCT(CT, CTx, CTy, CTz , Xint, Yint, Zint , HUextpl , CTintSize , CTintSpacing)|
%
%
%% Description
% |CTintrp = interpolateCT(CT, CTx, CTy, CTz , Xint, Yint, Zint , HUextpl , CTintSize , CTintSpacing)| Description
%
%
%% Input arguments
% |CT| - _SCALAR MATRIX_ - |CT(i,j,k)=HU| gives the intensity (Hounsfield unit) of the voxel at location (CTx(i),CTy(j),CTz(k))
%
% |CTx| -_SCALAR VECTOR_- |CTx(i)| Coordinate (mm) of the i-th voxel along the axis of the first pixel index
%
% |CTy| -_SCALAR VECTOR_-  |CTy(j)| Coordinate (mm) of the j-th voxel along the axis of the second pixel index
%
% |CTz| -_SCALAR VECTOR_-  |CTz(k)| Coordinate (mm) of the k-th voxel along the axis of the third pixel index
%
% |Xint| -_SCALAR VECTOR_- |Xint(u)| Coordinate (mm) of the u-th voxel of the interpolated CT scan
%
% |Yint| -_SCALAR VECTOR_- |Yint(u)| Coordinate (mm) of the u-th voxel of the interpolated CT scan
%
% |Zint| -_SCALAR VECTOR_- |Zint(u)| Coordinate (mm) of the u-th voxel of the interpolated CT scan
%
% |HUextpl| -_SCALAR_- HU value to place in the voxels that are extrapolated (outside the box of the original CT scan)
%
% |CTintSize| -_SCALAR VECTOR_- Number of pixels along each axis of the interpolated CT scan. prod(CTintSize) = numel(Xint)*numel(Yint)*numel(Zint)
%
% |CTintSpacing| -_SCALAR VECTOR_- Pixel size (mm , (x,y,z) ) of the interpolated CT scan
%
%% Output arguments
%
% |CTintrp| - _SCALAR MATRIX_ - Interpolated CT |CTintrp|. This is a 3D matrix with size |CTintSize|. The orering of the voxels is determined by the ordering of the voxels in Xint, Yint, Zint
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function CTintrp = interpolateCT(CT, CTx, CTy, CTz , Xint, Yint, Zint , HUextpl , CTintSize , CTintSpacing , interMethod)

  if nargin < 11
    interMethod = 'nearest';
  end

  %If the interpolated CT has larger voxels than the original CT scan:
  % Rebin the voxels of the high resolution CT scan to a pixel size close to the target resolution
  % The larger voxel will be the average of the intensities of all the small voxels contained in the box
  %the rebinning allows a reduction of the statistical noise from the Monte Carlo simulations
  % resCTorg    = [getResol(CTx) , getResol(CTy) , getResol(CTz) ]; % resolution of original CT scan
  % fac = max(resCTorg ./ CTintSpacing,[],'all'); %magnification factor between the 2 CT scans


  % if fac < 1
  %
  %   fprintf('Gaussian kernel \n')

    % The interpolated CT scan has a lower resolution.
    % Apply a Gaussian filter in order to remove the high frequency noise
    % The CT keeps the same size but the noise is filtered out
    % CT = GaussianFilter(fac , 'X' , CT);
    % CT = GaussianFilter(fac , 'Y' , CT);
    % CT = GaussianFilter(fac , 'Z' , CT);

    % CT = imresize3(CT,fac,'box'); %'box' compute the average of the voxel intensity
    % resCTorg = resCTorg ./ fac; %resolution of the re-binned CT scan
    % sCTorg = size(CT); %Number of voxels in the rebinned CT
    % CTx = CTx(1) + resCTorg(1) .* (0:sCTorg(1)-1); %Definition of the axes of the re-binned CT scan
    % CTy = CTy(1) + resCTorg(2) .* (0:sCTorg(2)-1);
    % CTz = CTz(1) + resCTorg(3) .* (0:sCTorg(3)-1);
  % end

  [Yct,Xct,Zct] = meshgrid(CTy ,CTx, CTz ); %Meshgrid inverse X and Y. And it is inversed again in interp3 for Xct,Yct
  CTintrp = interp3(double(Yct),double(Xct),double(Zct), CT, double(Yint), double(Xint), double(Zint),interMethod,HUextpl); %The Yint, Xint in interp3
  CTintrp = reshape(CTintrp, CTintSize(1) , CTintSize(2), CTintSize(3));


end


%==============================================
% Find the smallest gap between 2 values in a list of value
% Round to 2 digits
%==============================================
function res = getResol(X)
  delta = diff(unique(round(X,2))); %List of all delta between the values, rounded to 2 digits
  delta(delta==0) = []; %remove the zeros
  res = min(delta); %The smallest gap between é values of the series X

end


%========================================================
% Apply a 1D Gaussian filter to filter out the high frequency noise
% out of the image before resampling it
%
% INPUT
% |downfactor| -_SCALAR_- Ratio lowRes spacing / HiRes spacing
% |CoorAxis| -_STRING_- Name of the axis on which the filtering is to be applied 'X' , 'Y' or 'Z'
% |myImage| -_SCALAR MATRIX_- |myImage(x,y,z)| High resolution image on which the filter is to be applied
%
% OUTPUT
% |myImage| -_SCALAR MATRIX_- |myImage(x,y,z)| Filtered high resolution image
%========================================================
function myImage = GaussianFilter(downfactor , CoorAxis , myImage)

    sigma = downfactor*.4;
    fsz = round(sigma * 5);
    fsz = fsz + (1-mod(fsz,2));
    filter = gaussian_kernel(fsz, sigma);
    filter = filter/sum(filter);

    switch CoorAxis
      case 'X'
          paddingSize = [length(filter) 0 0];
      case 'Y'
          paddingSize = [0 length(filter) 0];
      case 'Z'
          paddingSize = [0 0 length(filter)];
    end


    myImage = padarray(myImage, paddingSize, 'replicate');
    myImage = conv3f(myImage, single(filter));

    switch CoorAxis
      case 'X'
          myImage = myImage(length(filter)+1:end-length(filter), :, :);
      case 'Y'
          myImage = myImage(:,length(filter)+1:end-length(filter), :);
      case 'Z'
          myImage = myImage(:,:,length(filter)+1:end-length(filter));
    end

end
