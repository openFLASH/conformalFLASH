%% exportBinaryCEF
% Export the conformal Energy Modulator into a binary file
%
%% Syntax
% |binIMGinfo = exportBinaryCEF(handles, Plan, b , outname, format)|
%
%
%% Description
% |binIMGinfo = exportBinaryCEF(handles, Plan, b , outname, format)| Description
%
%
%% Input arguments
%
% |Plan| -_STRUCTURE_- Information about the treament plan
%   * |Plan.Beams(b).RangeModulator.CEM3Dmask| -_SCALAR MATRIX_- 3D mask of the CEF. |CEM3Dmask(x,y,z)=1| if the voxel at location (x,y,z)  in the plane of the CEF for beam b belongs to the CEF.
%                                               Z=0 at the base of CEF. Z increase in the smae way as Zg if the spike point toward the proton source
%   * |Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing| -_SCALAR VECTOR_- |CompensatorPixelSpacing = [x,y,z]| Pixel size (mm) in the plane of the CEF for the |CompensatorThicknessData| matrix
%   * |Plan.Beams(b).RangeModulator.ModulatorOrigin| -_SCALAR VECTOR_- Physical coordinate [x,y,z] the voxel |CompensatorThicknessData(1,1)| and |hedgehog3D(1,1,1)| for beam b  in the plane of the CEF.
%   * |Plan.Beams(b).RangeModulator.FrameOfReferenceUID| -_STRING_- DICOM UID of the frame of reference for the CEM
%   * |Plan.Beams(b).RangeModulator.ReferencedRangeModulator| -_STRING_- DICOM UID of the binary file describing the CEM
%
%
%% Output arguments
%
% |binIMGinfo| - _STRUCTURE_ -  Description
%
%
%% Contributors
% Authors : L. Hotoiu (open.reggui@gmail.com)

function binIMGinfo = exportBinaryCEF(Plan, b , outname, format)

    % Generate unique UIDs for each binary CEM file
    binIMGinfo = Plan.CTinfo;
    binIMGinfo.Spacing = Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing; % This is the pixel size in the exported CEM. Not the pixel size of the original CT scan associated to the plan
    binIMGinfo.ImagePositionPatient = Plan.Beams(b).RangeModulator.ModulatorOrigin; %This is the origin of the CEM file. Not of the original CT scan associated to the plan
    binIMGinfo.FrameOfReferenceUID = Plan.Beams(b).RangeModulator.FrameOfReferenceUID;
    binIMGinfo.SeriesInstanceUID = Plan.Beams(b).RangeModulator.ReferencedRangeModulator;
    img = Plan.Beams(b).RangeModulator.CEM3Dmask;

    % Save the binary image to disk in uint8
    binIMGinfo.OriginalHeader.Modality = '';
    fileName = [outname, '_', strrep(binIMGinfo.SeriesInstanceUID,'.','')];
    binIMGinfo = save_Image(img, binIMGinfo, fileName, format);

    % Write header data to text file
    load([fileName, '_hdr.mat']); % out is the name of the loaded header variable
    fileID = fopen([fileName, '_header.txt'], 'w');
    [~,~,machinefmt,encodingOut] = fopen(fileID);

    fprintf(fileID, '%s\t', 'Encoding');
    fprintf(fileID, ' %s', encodingOut);
    fprintf(fileID, '\n%s\t', 'Bytes order');
    fprintf(fileID, ' %s', machinefmt);
    fprintf(fileID, '\n%s\t', 'Data type');
    fprintf(fileID, ' %s', 'uint8');
    fprintf(fileID, '\n%s\t', 'CEF data reading order = ');
    fprintf(fileID, ' %s', 'output per XY slice, in +Z IECGantry direction, such that the CEM base XY slice is written last, at the top of the file');
    fprintf(fileID, '\n%s\t', 'CEM material');
    fprintf(fileID, ' %d', Plan.Spike.MaterialID);
    fprintf(fileID, '\n%s\t', 'Spacing');
    fprintf(fileID, ' %d', Plan.Beams(b).RangeModulator.Modulator3DPixelSpacing);
    fprintf(fileID, '\n%s\t', 'ImagePositionPatient');
    fprintf(fileID, ' %d', Plan.Beams(b).RangeModulator.ModulatorOrigin);
    fprintf(fileID, '\n%s\t', 'PatientOrientation');
    fprintf(fileID, ' %d', out.info.PatientOrientation);
    fprintf(fileID, '\n%s\t', 'PatientID');
    fprintf(fileID, ' %s', out.info.PatientID);
    fprintf(fileID, '\n%s\t', 'FrameOfReferenceUID');
    fprintf(fileID, ' %s', out.info.FrameOfReferenceUID);
    fprintf(fileID, '\n%s\t', 'SeriesInstanceUID');
    fprintf(fileID, ' %s', out.info.SeriesInstanceUID);
    fprintf(fileID, '\n%s\t', 'StudyInstanceUID');
    fprintf(fileID, ' %s', out.info.StudyInstanceUID);
    fprintf(fileID, '\n%s\t', 'Type');
    fprintf(fileID, ' %s', out.info.Type);
    fprintf(fileID, '\n%s\t', 'Size');
    fprintf(fileID, ' %d', out.info.Size);

    fprintf('Export done\n');
    fclose(fileID);
end
