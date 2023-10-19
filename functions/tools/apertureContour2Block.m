%% apertureContour2Block
% Compute the coordinates (in IEC gantry) of the pixels containing brass to define the aperture block
%
%% Syntax
% |[Abrass , BlockThickness , IsocenterToBlockTrayDistance , Rmax] = apertureContour2Block(Beam ,  Spacing  , BDL_file)|
%
%
%% Description
% |[Abrass , BlockThickness , IsocenterToBlockTrayDistance , Rmax] = apertureContour2Block(Beam ,  Spacing  , BDL_file)| Description
%
%
%% Input arguments
%
% |Beam| -_STRUCTURES_- Information about the beam
%     * |Beam.Iso2Skin| -_SCALAR_- Isocentre (mm) To Skin Distance along the proton beam axis
%     * |Beam.Layers(L).Energy| -_SCALAR_- Energy (MeV) of the L-th layer
%     * |Beam.BlockThickness| -_SCALAR_- Thickness of the aperture block in mm
%     * |Beam.IsocenterToBlockTrayDistance| -_SCALAR_- Distance (mm) from isocentre to upstream side of aperture block
%     * |Beam.BlockData|  -_SCALAR MATRIX_- |ApertureBlockData(i,:)=[x,y]|  Coordinates (mm) of the i-th point defining the contour of the aperture block projected onto the machine isocentric plane in the IEC BEAM LIMITING DEVICE coordinate system.
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size ([x,y,z] in mm) of the pixels of the CT scan in which the aperture will be inserted
%
% |BDL_file| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%
%% Output arguments
%
% |Abrass| -_SCALAR MATRIX_- |A(i,:) = [x,y,z,1]| Matrix of 4 vector [x,y,z,1] defining the position of the i-th voxel ( mm, in IC gantry CS) of the voxels belonging to the device
%
% |Rmax| -_SCALAR_- Radius (mm) of the circle circumbscribing the square aperture block
%
%
%%REFERENCE
% [1] https://nl.mathworks.com/help/images/classify-pixels-that-are-partially-enclosed-by-roi.html
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Abrass , Rmax] = apertureContour2Block(Beam ,  Spacing  , BDL_file)

  [minSize , maxSize , Rmax, BlockBorder] = findApertureBlockSize(Beam.BlockData ); % Find the minimum and maximum coordinate of the brass block
  pxlSize = min(Spacing)./3;
  sAgntr = ceil((maxSize - minSize) ./ pxlSize); %Number of pixels of the aperture mask

  fac = 3; %Make the grid larger
  Apmask  = zeros(sAgntr); %Create a mask of the air hole in the brass block

  %Get the magnification foctor from isocentre to aperture plane
  mag = magnification(Beam.IsocenterToBlockTrayDistance , BDL_file); %Get the magnification factor from isocentre plane to upstream side aperture plane
  fprintf('Isocenter to Upstream Block side Distance : %f mm \n', Beam.IsocenterToBlockTrayDistance);
  fprintf('Magnification factor X : %f \n', mag(1))
  fprintf('Magnification factor Y : %f \n', mag(2))

  for cntIdx = 1:numel(Beam.BlockData)

    %In the treatment plan, the aperture contour is defined in the isocentre plane
    [~, x0  , y0 ] = DICOM2PXLindex([] , [pxlSize , pxlSize , 1] , [minSize , 0] , true , 0 , 0 , 0); %The magnification center is at pixel index [x0,y0]
    [~ , Xp , Yp ] = DICOM2PXLindex([] , [pxlSize , pxlSize , 1] , [minSize , 0] , true , Beam.BlockData{cntIdx}(:,1) , Beam.BlockData{cntIdx}(:,2) , zeros(numel(Beam.BlockData{cntIdx}(:,2)),1));

    %In the CT san, the aperture contour is projected (de-magnified) into the aperture plane
    Xmask = double(mag(1) .* (Xp - x0) + x0 ).* fac - floor(fac./2);
    Ymask = double(mag(2) .* (Yp - y0) + y0 ).* fac - floor(fac./2);

    ApmaskTrans = poly2mask(  Ymask , Xmask , sAgntr(1).*fac,  sAgntr(2).*fac ); %poly2mask gives a transposed mask. So inverse X and Y in the input parameters
            %Because of the pixels partially enclosed in the mask, poly2mask may make the aperture mask smaller because it drops the voxels
            %on the border of the mask
            %To avoid the problem, poly2mask fills the mask in a matrix that is 3* larger than the final image
            %When the image is rescaled to smaller pixel sizee, any partially filled voxels will still be selected in the final image
    ApmaskTrans = imresize3(single(ApmaskTrans),1./fac,'box'); %'box' compute the average of the pixel intensity
    Apmask = ApmaskTrans | Apmask; %Mask of the air hole

  end

  %Dilate the border of the aperture
  % Apmask is changed from the air hole into the brass block
  Apmask = addBorder2mask(Apmask , BlockBorder , [pxlSize , pxlSize]);
  Xvec = 1:sAgntr(1);
  Yvec = 1:sAgntr(2);
  Xvec =  minSize(1) + Xvec .* pxlSize;
  Yvec =  minSize(2) + Yvec .* pxlSize;
  [Y,X] = meshgrid(Yvec,Xvec);

  % Apmask is the mask of the brass block
  Ind = find(Apmask); %Coordinate of the voxels containing brass
  Xa = X(Ind);
  Ya = Y(Ind);
  Ap = [Xa,Ya]; %X,Y coordinate of voxels containing brass
  Zres = Spacing(3) ./ 5; %Z resolution at Nyquest limit along the Z axis. This will avoid having gaps between layers duye to pixelisation


  %Create a x * 4 matrix with the coordinate of all the CT voxel to be filled by brass
  % IsocenterToBlockTrayDistance defines the upstream side of the block which is at more positive values of Z IEC gantry.
  Z = Beam.IsocenterToBlockTrayDistance-Beam.BlockThickness : Zres : Beam.IsocenterToBlockTrayDistance;
  Z2 = reshape(repmat(Z ,  size(Ap,1) ,1) , [1,size(Ap,1).* numel(Z)]);

  Abrass = ones(size(Ap,1).* numel(Z) , 4); %Indices of the voxels containing brass
  Abrass(:,1:2) = repmat(Ap , numel(Z) , 1);
  Abrass(:,3) = Z2' ;


end
