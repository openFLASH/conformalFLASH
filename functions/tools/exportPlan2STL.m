%% exportPlan2STL
% Export the CEM inside flash plan to a STL file on disk
%
%% Syntax
% |exportPlan2STL(plan_filename)|
%
%
%% Description
% |exportPlan2STL(plan_filename)| Function exporting the CEM inside a flashplan to STL file on disk in the same folder as the dicom plan
%
%% Input arguments
%
% |plan_filename| -_STRING_- location and name of flash plan on disk
%
%% Output arguments
%
% None
%
% TODO STL files do not specify what units their distances are in. When a program opens a STL file,
% it only knows that the model measures a certain number of units in each dimension. [3]
% [4] There is no way to convey the units used in an STL file. You do need to specify units when you generate the file,
% but you also need to know what units were used when you import it, as units are not specified in the file itself.
% In the 3D printing world, this generally doesn't cause too many problems, as everyone just specifies mm,
% since this is what most printers understand natively anyhow. STEP files for instance contain the units used in the file,
% so objects always import correctly, unless you override the units. Alas, this is not the case for STL files.
%
%% Contributors
% Authors : L. Hotoiu (open.reggui@gmail.com)

function exportPlan2STL(plan_filename)

    [plan_filepath, plan_name, plan_ext] = fileparts(plan_filename);

    handles = Initialize_reggui_handles();
    handles.dataPath = plan_filepath;

    %Plan.ScannerDirectory = 'default';
    Plan.ScannerDirectory = 'D:\MATLAB\REGGUI\plugins\openMCsquare\lib\Scanners\default';
    Plan.showGraph = true;
    Plan.BDL = 'D:\MATLAB\REGGUI\plugins\openMCsquare\lib\BDL\BDL_default_UN3_Al_RangeShifter_tilted_IBA_FLASH_G0_round_UpennFlash_v02.txt';
    [handles, Plan] = parseFLASHplan(plan_filename, Plan, handles);


    %Add voxels all around the base in order to close the STL object
    CEM3Dmask1 = Plan.Beams.RangeModulator.CEM3Dmask;
    CEM3Dmask = zeros(size(CEM3Dmask1,1) + 2, size(CEM3Dmask1,2) + 2, size(CEM3Dmask1,3) + 2);
    CEM3Dmask(2:end-1, 2:end-1, 2:end-1) = CEM3Dmask1;
    origin = Plan.Beams.RangeModulator.ModulatorOrigin - Plan.Beams.RangeModulator.Modulator3DPixelSpacing;

    stl_filename = fullfile(plan_filepath, [plan_name '.stl']);
    exportCEM2STL(CEM3Dmask, Plan.Beams.RangeModulator.Modulator3DPixelSpacing, origin, Plan.Beams.RangeModulator.AccessoryCode, Plan.Beams.RangeModulator.ModulatorMountingPosition, stl_filename)
end
