%% convScatter
% Convolve the fluence map |Fluence(x,y,T)| with the Gaussian kernel
% representing the lateral scatter due to the material thickness |thickness(T)|.
% Use Moliere theory to compute the sigma of the lateral scatter.
%
% The output is |Fluence(x,y,T)| where the index T refers to the |thickness(T)| of the material.
% After traversing |thickness|, the proton beam travels |distance|, which causes further increase of the spot sigma, due to beam divergence
% The output map is computed at this |distance| away from the base of the scatterer.
%
% By default, all the fluence map are treated separatly for each scatterer thickness.
% If a |maskLayer(x,y,T)| is provided, then the scattering between thickness layer is taken into account: some of the fluence can be transfered into a different map.
% The proton scattered into air are assumed to continue their journey with within the same energy layer |Fluence(x,y,T)|
% Only the proton scattered into a TALLER step need to be tracked.
% |thickness| is sorted in ascending order: the fluence of the proton scatter outside of the |maskLayer(x,y,T)| into the |maskLayer(x,y,T+k)| are allocated to |Fluence(x,y,T+k)|
%
%
%% Syntax
% |[Fluence , sigmas] = convScatter(Fluence, sigmas, maskLayer)|
%
%
%% Description
% |[Fluence , sigmas] = convScatter(Fluence, sigmas, maskLayer)| Compute fluence map
%
%
%% Input arguments
% |Fluence(x,y,T)| -_SCALAR MATRIX_-  |Fluence(x,y,T)| Incident proton fluence at posiiton (x,y) for proton that will go through |thickness(T)| of material
%
% |sigmas| -_SCALAR VECTOR_- [OPTIONAL. If absent, the sigma are computd using Moliere angle] |sigmas(T)|  Sigma (mm) of the Gaussian kernel for lateral scatter for T-th thickness of scatterer
%
% |maskLayer(x,y,T)| -_SCALAR MATRIX_- [OPTIONAL : default, treat each fluence map separately] |maskLayer(x,y)=1| indicate the cross-section of the step corresponding to energy layer E
%
%
%% Output arguments
%
% |Fluence(x,y,E)| -_SCALAR MATRIX_-  |Fluence(x,y,E)| Proton fluence at |Distance| for posiiton (x,y) for proton of energy E
%
% |sigmas| -_SCALAR VECTOR_- |sigmas(E)|  Sigma (mm) of the Gaussian kernel for lateral scatter for E-th energy layer
%
%% REFERENCE
% [1] https://blogs.mathworks.com/steve/2017/01/23/gaussian-filtering-with-imgaussfilt/
% [2] https://www.crisluengo.net/archives/150/
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Fluence , sigmas] = convScatter(Fluence, sigmas, maskLayer)

  %Get lateral spread sigma just after spike and at |distance|
  if nargin < 2
    sigmas = [];
  end


  if nargin < 3
    maskLayer = []; %We did not receive any mask to scatter protons into lower energy layer
  end


  N_thickn = size(Fluence,3);
  smask = [size(Fluence,1) , size(Fluence,2)];

  %Convolution of the scattering kernel with the fluence map
  %The layers are ordered with indices increasing from outside ring towards inside ring
  %As the scattered protons are re-assigned to the lower energy layer, we shall process from the highest energy layer
  %towards the lower energy layer in order to sequentially apply the convolution to previously scattered proton
  [~,~,maskE]=meshgrid(ones(1,size(Fluence,2)),ones(1,size(Fluence,1)),1:N_thickn); %Energy mask. Allow selecting energy layers with a simple conditional check. Invert X and Y fluence coordinates for meshgrid quirk
              %NB: first input de meshgrid becomes second index of the outp ut matrix and vice versa

  for k = 1:N_thickn
      %Loop for each energy layer
      sigma = sigmas(k); %sigma in pxl

      if (sigma > 0) && numel(find(Fluence(:,:,k)))
        %Bother with convolution only if the scatter is significative and if there is some fluence

        FilterSize = double(sigma .* 6 + 1); %Filter size must be odd. 3*sigma takes 92% of Gaussian. *2 to get diameter
        h = fspecial('gaussian',[1,FilterSize],double(sigma)); % 1D Gaussian kernel
        FluenceConv = conv2(h,h,Fluence(:,:,k),'same'); %Following [1], it is faster to apply 2 1D convolutions for a symetrical kernel

        if ~isempty(maskLayer)

            scatterMask = zeros(smask);
            % The proton scattered into air are assumed to continue their journey with within the same energy layer
            % Only the proton scattered into a higher step need to be tracked
            %If the scattered photon spread over several narrow step, allocate to all the steps
            Fluence   = Fluence + FluenceConv .* maskLayer .* (maskE >= k+1) ; %MAttribute the scattered proton to all the lower energy layers

            for Ek = k+1  : N_thickn
              scatterMask = scatterMask + maskLayer(:,:,Ek);
            end

            Fluence(:,:,k)   = Fluence(:,:,k) + FluenceConv .* ~scatterMask; %MAttribute the un-scattered proton to the k-th energy layer

          else
            %Do not remove the contribution from the central hole of the steps
            Fluence(:,:,k) = FluenceConv;
          end
      else
          % The Moliere scatter and the beam divergence are small. The fluence remain unchanged after going through the step
          %No need to apply convolution
          %Fluence(:,:,k) = Fluence(:,:,k); %This line does not do anything usefull
      end

  end %for k
end
