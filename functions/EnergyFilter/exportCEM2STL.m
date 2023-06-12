%% exportCEM2STL
% Export the CEM 3D binary mask to a STL file on disk
%
%% Syntax
% |exportCEM2STL(CEF3dMask , pixelSize , origin , AccessoryCode , filename)|
%
%
%% Description
% |exportCEM2STL(CEF3dMask , pixelSize , origin , AccessoryCode , filename)| Description
%
%
%% Input arguments
% |CEF3dMask| -_SCALAR MATRIX_- |CEF3dMask(x,y,z)| 3D mask of the CEF.  |CEF3dMask(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEF for beam b belongs to the CEM.
%
% |pixelSize| -_SCALAR MATRIX_- |pixelSize = [x,y,z]| Pixel size of the CEM. The CEM is specified in the plane of the CEM
%
% |origin|  -_SCALAR MATRIX_- |origin= [x,y,z]| Physical coordinate the voxel |CEF3dMask(1,1,1)| (in IEC gantry) in the plane of the CEF.
%
% |AccessoryCode| -_STRING_- (300A,00F9) An accessory identifier to be read by a device such as a bar code reader. Read from the YAML configuration file
%
% |ModulatorMountingPosition| -_STRING_- DICOM tag Modulator Mounting Position (300D,1017). Defines the face of the CEM used as reference for |IsocenterToRangeModulatorDistance|
%                         SOURCE_SIDE is using the downstream face of the block as a reference position for expressing the isocenter to block tray distance.
%                         PATIENT_SIDE is using the upstream face
%
% |filename| -_STRING_- file name (including full path and extension) of the STL file
%
%% Output arguments
%
% None
%
%TODO STL files do not specify what units their distances are in. When a program opens a STL file,
% it only knows that the model measures a certain number of units in each dimension. [3]
% [4] There is no way to convey the units used in an STL file. You do need to specify units when you generate the file,
% but you also need to know what units were used when you import it, as units are not specified in the file itself.
% In the 3D printing world, this generally doesn't cause too many problems, as everyone just specifies mm,
% since this is what most printers understand natively anyhow. STEP files for instance contain the units used in the file,
% so objects always import correctly, unless you override the units. Alas, this is not the case for STL files.
%
%% Contributors
% Authors :R. Labarbe, L. Hotoiu (open.reggui@gmail.com)


function exportCEM2STL(CEM3DMask, pixelSize , origin , AccessoryCode, ModulatorMountingPosition, filename)

    fprintf('Exporting CEM to STL... \n')
    intrpPxlSize = 0.2;

    if strcmp(ModulatorMountingPosition, 'PATIENT_SIDE')
        signZ = -1;
    elseif strcmp(PModulatorMountingPosition, 'SOURCE_SIDE')
        signZ = 1;
    else
        error('Modulator mounting position is incorrectly specified in dicom\n');
    end

    % Prepare data for isosurface/meshgrid/imagesc YX swap. The usual matlab YX axis inversion
    CEM3DMask = permute(CEM3DMask, [2,1,3]);
    CEMPixelSize = [pixelSize(2), pixelSize(1), pixelSize(3)]; 
    CEMPxlfac = round(CEMPixelSize ./ intrpPxlSize); %One input pixel is divided in multiple smaller pixels

    % Resize the original 3D mask coming for dicom plan
    CEM3Dmask_full = imresize3(CEM3DMask, size(CEM3DMask).*CEMPxlfac,'nearest');

    % (Re)Compute the new pixel size to refine STL export
    intrpPxlSize = CEMPixelSize ./ CEMPxlfac;
    fprintf('CEM pixel size for exporting in STL : (%3.1f, %3.1f, %3.1f) mm \n', intrpPxlSize(1), intrpPxlSize(2), intrpPxlSize(3));
    intrpPxlSize = intrpPxlSize(1);

    fprintf('Computing surface triangles \n')
    %FV = surf2solid(Xm, Ym, Vq, 'ELEVATION', 0, 'THICKNESS', intrpPxlSize/2); %Create vertexes on triangular grid
    FV = isosurface(CEM3Dmask_full, 0.5, 'noshare'); %Compute the vertices in index coordinates
    FV.vertices(:,1) = (FV.vertices(:,1)-1) .* CEMPixelSize(1)./CEMPxlfac(1) + origin(1) - ((CEMPxlfac(1)-1)/2).*intrpPxlSize; %Convert from index coordinates into physics coordinates
    FV.vertices(:,2) = (FV.vertices(:,2)-1) .* CEMPixelSize(2)./CEMPxlfac(2) + origin(2) - ((CEMPxlfac(2)-1)/2).*intrpPxlSize;
    FV.vertices(:,3) = (FV.vertices(:,3)-1) .* CEMPixelSize(3)./CEMPxlfac(3) + origin(3) - ((CEMPxlfac(3)-1)/2).*intrpPxlSize;
    FV.vertices(:,3) = signZ .* FV.vertices(:,3);
    fprintf('Done \n')


    fprintf('Saving STL file to %s \n',filename);
    NbDigit = ceil(-log10(intrpPxlSize));
    sizeXYZ = max(FV.vertices, [], 1) - min(FV.vertices, [], 1);
    PhysSize = round(sizeXYZ, NbDigit); %mm length of the hedgehog;
    header = [AccessoryCode, '-- Unit : mm X=', num2str(PhysSize(1)), 'mm Y=', num2str(PhysSize(2)), 'mm Z=', num2str(PhysSize(3)), 'mm'];
    if (numel(header) > 80)
        header = header(1:80); %Header must be less than 80 characters in STL standard
    end

    filepath = fileparts(filename);
    if (~exist(filepath,'dir'))
        %The folder to save the CT does not exist. Create it
        mkdir (filepath)
    end
    stlwrite(filename, FV.faces , FV.vertices , 'TITLE' , header); %Create the STL file
    fprintf('Done \n')

    %Export ridge filter structural coordinates to text file
    % coordFilename = [fullfile(pathName,Plan.Beams(b).RangeModulator.AccessoryCode) '.i'];
    % fprintf('Saving structural coordinates to .txt file at %s \n', coordFilename)
    % exportCEFtoMCNPX(coordFilename, GridLayout, SpikeCentres, BaseSize, SpikeSteps, TowerApothems, Plan.Spike.SpikeType);
    %fprintf('Done \n')
end
