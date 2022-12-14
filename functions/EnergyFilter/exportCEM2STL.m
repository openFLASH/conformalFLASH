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
% Authors : R. Labarbe (open.reggui@gmail.com)

function exportCEM2STL(CEF3dMask , pixelSize , origin , AccessoryCode , filename)

      intrpCTpxlSize = min(pixelSize);

      fprintf('Computing surface triangles \n')
      %FV = surf2solid(double(Xbl ), double(Ybl ), CEFelevation , 'ELEVATION',0); %Create vertexes on triangular grid
      FV = isosurface(CEF3dMask , 0.5,'noshare','verbose'); %Compute the vertices in index coordinates
      FV.vertices(:,1) =  FV.vertices(:,1) .* pixelSize(1) - origin(1); %Convert from index coordinates into physics coordinates
      FV.vertices(:,2) =  FV.vertices(:,2) .* pixelSize(2) - origin(2);
      FV.vertices(:,3) =  FV.vertices(:,3) .* pixelSize(3);
      fprintf('Done \n')


      fprintf('Saving STL file to %s \n',filename);
      NbDigit = ceil(-log10(intrpCTpxlSize));
      PhysSize = round(size(CEF3dMask) .* pixelSize , NbDigit); %mm length of the hedgehog;
      header = [AccessoryCode , '-- Unit : mm X=',num2str(PhysSize(1)),'mm Y=',num2str(PhysSize(2)),'mm Z=',num2str(PhysSize(3)),'mm'];
      if (numel(header) > 80)
        header = header(1:80); %Header must be less than 80 characters in STL stnadard
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
    fprintf('Done \n')

end
