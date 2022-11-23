%% makeAperture
% Compute the chape of the aperture block for all the beams in the plan
%
%% Syntax
% |Plan = makeAperture(Plan , PTV , Body ,  Spacing , ImagePositionPatient )|
%
%
%% Description
% |Plan = makeAperture(Plan , PTV , Body ,  Spacing , ImagePositionPatient )| Description
%
%
%% Input arguments
%
% |Plan| -_STRUCTURE_- Information about the treament plan
%  * |Plan.DoseGrid| - _struct_ - Structure containing the information about the dose grid. The following data must be present:
%  * |Plan.DoseGrid.nvoxels| - _scalar_ - Total number of voxels of the dose grid, i.e., number of voxels in the CT image.
%  * |Plan.CTinfo| -_STRUCTURE_- DICOM header of the CT scan
% * |Plan.ScannerDirectory| - _STRING_ - Name of the folder containing the definition of the CT scanner properties in MCsquare in folder "plugins\openMCsquare\lib\Scanners"
% * |Plan.Beams| -_VECTOR of STRUCTURES_- Information about the different beams in the plan
%     * |Plan.Beams(b).GantryAngle| -_SCALAR_- Angle (deg) of the b-th beam in the plan
%     * |Plan.Beams(b).isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan for b-th beam
%     * |Plan.Beams(b).PatientSupportAngle| -_SCALAR_- Angle (deg) of the couch for the b-th beam in the plan
%     * |Plan.Beams(b).Iso2Skin| -_SCALAR_- Isocentre (mm) To Skin Distance along the proton beam axis
%     * |Plan.Beams(b).Layers(l).Energy| -_SCALAR_- Energy (MeV) of the l-th layer of the b-th beam
%
% |PTV| - _SCALAR MATRIX_ - Mask defining the position of the PTV |PTV(x,y,z)=1| if the voxel is inside the PTV
%
% |Body| - _SCALAR MATRIX_ - Mask defining the contour of the body. This is used to define the position of hte aperture |PTV(x,y,z)=1| if the voxel is inside the PTV
%
% |Spacing| - _SCALAR VECTOR_ - Pixel size (|mm|) of the pixels defining the aperture block
%
% |ImagePositionPatient| - _SCALAR VECTOR_ - Coordinate (in |mm|) of the first pixel of the image in the coordinate system of the image
%
%
%% Output arguments
%
% |Plan| -_STRUCTURE_- Updated treatment plan.
%   * |Plan.Beams(b).IsocenterToBlockTrayDistance| -_SCALAR_- Distance (mm) from isocentre to upstream surface of aperture block
%   * |Plan.Beams(b).Iso2Skin| -_SCALAR_- Isocentre (mm) To Skin Distance along the proton beam axis
%   * |Plan.Beams(b).BlockData|  -_SCALAR MATRIX_- |BlockData(i,:)=[x,y]|  Coordinates (mm) of the i-th point defining the contour of the aperture block projected onto the machine isocentric plane in the IEC BEAM LIMITING DEVICE coordinate system.
%   * |Plan.Beams(b).BlockThickness| -_SCALAR_- Thickness of the aperture block in mm
%   * |Plan.Beams(b).BlockMaterialID| - _STRING_ - Name of the aperture block material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%   * |Plan.Beams(b).SnoutPosition| -_SCALAR_- Distance (mm) from the isocenter to the downstream surface of tha aperture block
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Plan = makeAperture(Plan , PTV , Body ,  Spacing , ImagePositionPatient )

  for b = 1: size(Plan.Beams,2) %Loop for each beam

    if Plan.Beams(b).ApertureBlock
        %Dilate the PTV with the margin
        if (isfield(Plan.Beams(b), 'ApertureMargin'))
          if (prod(Plan.Beams(b).ApertureMargin) ~= 0)
            fprintf('Adding margin to PTV for lateral fall-off \n')
            %If a margin is defined, expand the PTV volume to be sure to have a sufficient hole in the aperture
            PTV = DilateROI(PTV , Plan.Beams(b).ApertureMargin , Spacing);
          else
            fprintf('No margin to PTV \n')
          end %f (prod(Plan.
        end %if (isfield(Plan

        %Build the aperture in IEC gantry
        %Project the aperture on a plane perpendicualr to proton axis
        %Use a sufficient density of pixels to avoid rounding errors and be sure to remove all brass pixels in hole
        fprintf('Projecting PTV on isocentre plane\n')
        [Agntr , X, Y , pxlSize] = projectROIOnIsoplane(Plan.Beams(b) , PTV , Spacing , ImagePositionPatient , min(Spacing)./3);

        sAgntr = size(Agntr);

        %Draw aperture contour to export in plan
        Apmask  = zeros(sAgntr); %Create a mask of the air hole in the brass block
        [B,L,n] = bwboundaries(Agntr); %Contours of the border between air and brass. In the isocentre plane.
                  %B{}(x,y) The first index of B is also the first index of Agntr

        for cntIdx = 1:numel(B)
          bnd = B{cntIdx};
          %In the treatment plan, the aperture contour is defined in the isocentre plane
          ApertureBlockData{cntIdx} = [X(sub2ind(size(X),bnd(:,1),bnd(:,2))) , Y(sub2ind(size(Y),bnd(:,1),bnd(:,2)))]; %List of (X,Y) points defining the contour of the aperture block
        end %for cntIdx

        %Compute the minimum distance from isocentre to aperture
        %to avoid contact with patient skin
        [~ , ~ , Rmax] = findApertureBlockSize(ApertureBlockData );
        Iso2Skin = getIsocentreToSkinDistance(Plan.Beams(b) , Body , Rmax , Spacing , ImagePositionPatient) + 2; %Isocentre (mm) To Skin Distance along the proton beam axis

        Plan.Beams(b).Iso2Skin = Iso2Skin;  %Isocentre (mm) To Skin Distance along the proton beam axis
        [IsocenterToBlockTrayDistance , BlockThickness , ~ , BlockMountingPosition ] = getIsocenterToBlockTrayDistance(Plan.Beams(b));
        Plan.Beams(b).IsocenterToBlockTrayDistance =  IsocenterToBlockTrayDistance; %Distance (mm) from isocenter to upstrea surface of aperture block
        Plan.Beams(b).BlockMountingPosition = BlockMountingPosition;
        Plan.Beams(b).BlockThickness = BlockThickness;
        Plan.Beams(b).BlockData = ApertureBlockData;
        Plan.Beams(b).BlockMaterialID = 'BRASS';
        Plan.Beams(b).SnoutPosition = IsocenterToBlockTrayDistance; %Distance isocentre to UPSTREAM surface of the aperture
        % By defining the snout position at the same reference point as the iso to block tray distance,
        % both distances are equal and independent of the block thickness


        %Draw the aperture contour on the WET map (for information)
        if Plan.showGraph
          figure(100+b)
          hold on
          for cntIdx = 1:numel(ApertureBlockData)
            plot(ApertureBlockData{cntIdx}(:,1),ApertureBlockData{cntIdx}(:,2),'-r')
            hold on
            drawnow
          end
        end

      else

        %No aperture block is requested. Compute anyway the distance from isocentre to skin surface
        fprintf('No aperture\n')
        fprintf('Projecting PTV on isocentre plane\n')
        [Agntr , X, Y , pxlSize] = projectROIOnIsoplane(Plan.Beams(b) , PTV , Spacing , ImagePositionPatient , min(Spacing)./3);
        Xvec = unique(X(:));
        Yvec = unique(Y(:));
        sAgntr = size(Agntr);
        Aidx = find(Agntr);
        [Xidx,Yidx] = ind2sub(sAgntr , Aidx);

        %Compute the minimum distance from isocentre to aperture
        %to avoid contact with patient skin
        Rmax = sqrt(max(Xvec(Xidx).^2 + Yvec(Yidx).^2));
        [~ , OrigIdx] = min(X.^2 + Y.^2,[],'all','linear'); %Find the index of the pixel at the origin of the coordinate system
        Plan.Beams(b).Iso2Skin = getIsocentreToSkinDistance(Plan.Beams(b) , Body , Rmax , Spacing , ImagePositionPatient) + 2; %Isocentre (mm) To Skin Distance along the proton beam axis
        fprintf('Isocenter to Skin Distance : %f mm \n', Plan.Beams(b).Iso2Skin)
                      %Add Iso2Skin to Beams so that getIsocenterToBlockTrayDistance can use this distance

    end %if Plan.Beams(b).ApertureBlock

  end % for b

end
