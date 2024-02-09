% Contributors
% Authors : Lucian Hotoiu (open.reggui@gmail.com)
% -------------------------------------------------------------------------



% !!!!!!!! README !!!!!!!! ------------------------------------------------
%
% Prior to running the script we must rotate and crop the two images to be 
% compared. The object in the two images should have roughly the same 
% orientation, position and occupy similar areas in the image. These
% operations have to be well defined before running. 
%
% Note: If the scan is taken in the agreed-upon position then the code will 
% perform the needed operations automatically.
%
% !!!!!!!! README !!!!!!!! ------------------------------------------------



function CEFphysQA(configFile)

tic


%Inputs / Outputs

% Load the JSON file with the parameters for the computation
%-----------------------------------------------------------
config = loadjson(configFile);

scanCEF_path = config.files.scanCEF_path;
CT_air_path = config.files.ct_air_path;
RSplanFileName = config.files.RS_plan_fileName;
Plan.BDL = config.files.BDL;
Plan.ScannerDirectory = config.files.scanner_directory;
outputPath = config.files.outputPath;

PSF_FWHM = config.ct_scan_operations.PSF_FWHM; % in [mm]. Point Spread Function size as measured for the used CT scanner.
scan_mask_low_thres = config.ct_scan_operations.maskLowThres;
scan_mask_high_thres = config.ct_scan_operations.maskHighThres;

% Dshape scan orientation
% % scanCEF_path = 'C:\Users\lhotoiu\Downloads\20230628\D58_NO_bubbles\CT.1.3.12.2.1107.5.1.4.83552.30000023062822093699000003556.dcm';
permOrder1 = [1,3,2]; %Permutation 1 of the dimensions of CT scan to align it with plan
permOrder2 = [2,1,3]; %Permutation 2 of the dimensions of CT scan to align it with plan
permOrder = [permOrder1; permOrder2]; % input permutations in [perm1; perm3; etc] order
flipAxis = [2]; %Which axis index should be flipped (after permute); % input flips in [dim1; dim2; dim3] order

% Dog Ruthie Fowler scan orientation
% % scanCEF_path = 'D:\MATLAB\Data\Tests\Upenn\CanineTrial\CEMQA\geomQA\CEM_CT_RuthieFowler\CEM_CT_RuthieFowler\FOV_CEM1\2.16.840.1.114362.1.12209795.22564599136.663717879.118.3255.dcm';
% permOrder1 = [1,3,2]; %Permutation 1 of the dimensions of CT scan to align it with plan
% permOrder = [permOrder1]; % input permutations in [perm1; perm3; etc] order
% flipAxis = [3; 2]; %Which axis index should be flipped (after permute); input flips in [dim1; dim2; dim3] order

% %refCEF_path = 'E:\Lucian\Data\ScanCEFUpenn_sylvain\ctCEF505030_DCM\resample_333um\CEF_opt_resample.mhd';
% %stl_path = 'C:\Users\lhotoiu\Downloads\CEF1id_square_decimated_more.stl';
% -------------------------------------------------------------------------



% Enable additional functionalities
geomGammaIndexAnalysis = config.flags.geomGammaIndexAnalysis;
weplGammaIndexAnalysis = config.flags.weplGammaIndexAnalysis;


% Import reference CEF
source = 'Plan'; % 'CT'; 'STL';
% #########################################################################




% Initial image orientation assessment 
% #########################################################################

% Load the data
pre_handles = loadAlldatasets(scanCEF_path , CT_air_path , RSplanFileName , outputPath , Plan , permOrder, flipAxis, source);

% Define image names
scan_CEF_imageName = 'scan_CEF';
air_CT_imageName = 'air_CT';

% Display one slice of plan CEM for information
[pre_ref_cef_data, pre_ref_cef_info, ~] = Get_reggui_data(pre_handles, air_CT_imageName);

minTot = min(pre_ref_cef_data,[],'all');
maxTot = max(pre_ref_cef_data,[],'all');

%for i = 1:size(pre_ref_cef_data,3)
i = 30;
figure(10)
imagf = squeeze(pre_ref_cef_data(:,:,i));
image((imagf - minTot) ./ (maxTot  - minTot) .* 255)
title(['Reference -- Slice = ' num2str(i)])
xlabel ("Y #voxels");
ylabel ("X #voxels");
%pause
%end


% Display one scan slice for information
[pre_scan_cef_data, pre_scan_cef_info, ~] = Get_reggui_data(pre_handles, scan_CEF_imageName);

minTot = min(pre_scan_cef_data,[],'all');
maxTot = max(pre_scan_cef_data,[],'all');
%for i = 1:size(pre_scan_cef_data,3)
i=150;
figure(11)
imagf = squeeze(pre_scan_cef_data(:,:,i));
image((imagf - minTot) ./ (maxTot  - minTot) .* 255)
title(['Measurement -- Slice = ' num2str(i)])
xlabel ("Y #voxels");
ylabel ("X #voxels");
%pause
%end
% #########################################################################




% HU density analysis
% #########################################################################


% Compute grid for resampling
fprintf('Computing relevant region for probing and assesing scan object HU density \n');
pre_ref_cef_info.Spacing = pre_ref_cef_info.Spacing([1,3,2]); % invert 3 and 2 dimension according to reggui handles
pre_scan_cef_info.Spacing = pre_scan_cef_info.Spacing([1,3,2]); % invert 3 and 2 dimension according to reggui handles

grid.origin = [1; 1; 1];
grid.size = floor((size(pre_scan_cef_data)' ./ (pre_ref_cef_info.Spacing ./ pre_scan_cef_info.Spacing)));
grid.spacing = pre_ref_cef_info.Spacing;


%handles = Resample_all(handles,[1;1;1],[318;125;298],[0.6;0.5;0.6]); 


pre_handles = Resample_all(pre_handles, grid.origin, grid.size, grid.spacing);
[pre_scan_cef_data, pre_scan_cef_info, ~] = Get_reggui_data(pre_handles, scan_CEF_imageName);
%pre_scan_cef_data = medfilt3(pre_scan_cef_data, [7 7 7], 'replicate'); % very slow

% Threshold, erode and intersect
cef_mask = 'scan_cef_mask';
%pre_handles = AutoThreshold(scan_CEF_imageName, [128], cef_mask, pre_handles);
pre_handles = ManualThreshold(scan_CEF_imageName, [scan_mask_low_thres scan_mask_high_thres], cef_mask, pre_handles);
%pre_handles = Erosion(cef_mask, [2.5 2.5 2.5], 'eroded_cef_mask', pre_handles);
pre_handles = Erosion(cef_mask, [2.2 2.2 2.2], 'eroded_cef_mask', pre_handles);
[eroded_cef_mask_data, ~, ~] = Get_reggui_data(pre_handles, 'eroded_cef_mask');
%[cef_mask_data, cef_mask_info, ~] = Get_reggui_data(pre_handles, cef_mask);
intersect_data = pre_scan_cef_data .* eroded_cef_mask_data;
pre_handles = Set_reggui_data(pre_handles, 'intersect', intersect_data, pre_scan_cef_info, 'images', 1);
intersect_HU_data = nonzeros(intersect_data);



% Get mean and std of HU in scanned CEF. Do so inside eroded ROI to cut off
% image artefacts at the extremities of the object
mean_HU_scan = mean(intersect_HU_data, 'all', 'omitnan');
std_HU_scan = std(intersect_HU_data, [], 'all', 'omitnan');
min_HU_scan = min(intersect_HU_data, [], 'all', 'omitnan');
max_HU_scan = max(intersect_HU_data, [], 'all', 'omitnan');



nr_slices = size(intersect_data, 3);
offset = 15; % 50
%offset = 42; % 50
slice = offset + floor(nr_slices/10);
slice_increment = 3; %1
%slice_increment = 4; %1

figure(903);
tiledlayout(2,2);

for idx = 1:4
    nexttile
    contourf(squeeze(intersect_data(:,:,slice)),100,'LineColor','none');
    %contourf(squeeze(intersect_data(:,slice,:)),100,'LineColor','none');
    colorbar;
    %caxis([(mean_HU_scan - 3*std_HU_scan) (mean_HU_scan + 3*std_HU_scan)]);
    xticklabels(xticks * pre_handles.spacing(1))
    yticklabels(yticks * pre_handles.spacing(3))
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title(['Z plane ' num2str(slice)]);
    slice = slice+slice_increment;
end

sgtitle('Resample/eroded scan HU distribution in Z plane')



% plot volume histogram
figure(997);
%histogram(reshape(intersect_HU_data,size(intersect_HU_data,1)*size(intersect_HU_data,2)*size(intersect_HU_data,3),1,1),255);
histogram(intersect_HU_data,255);
xlabel('HU');
ylabel('# of voxels');
title('Resampled/eroded scan - HU histogram inside object only');
set(gca,'YScale','log')
% #########################################################################




% Geometrical analysis
% #########################################################################

%Load again the data as the previous data set was resampled in pre_handles
handles = loadAlldatasets(scanCEF_path , CT_air_path , RSplanFileName , outputPath , Plan , permOrder, flipAxis, source);

% Import scanned CEF
[scan_cef_data, scan_cef_info, ~] = Get_reggui_data(handles, scan_CEF_imageName);

% Save scan to disk
handles = save2Disk(handles, scan_cef_data, size(scan_cef_data), scan_cef_info, 'scan_CEF_image', outputPath);           


% Import ref CEF
ref_CEF_imageName = 'air_CT';
[ref_cef_data, ~, ~] = Get_reggui_data(handles, ref_CEF_imageName);

figure(9977)
histogram(ref_cef_data,255);
xlabel('HU');
ylabel('# of voxels');
title('HU histogram reference image, before registration');
set(gca,'YScale','log')


% Rigid registration normalised scan and reference cef
fprintf('Computing rigid registration between scan (fixed) and reference (moving) images \n');
ref_rigid_def = 'ref_rigid_def';
ref_rigid_trans = 'ref_rigid_trans';
handles = Registration_ITK_rigid_multimodal(scan_CEF_imageName, ref_CEF_imageName, ref_rigid_def, ref_rigid_trans, handles);


pxlSizeXY = scan_cef_info.Spacing(2);
pxlSizeZ = scan_cef_info.Spacing(3); % we approximate here that the voxels are isotropical. In reality there is a 0.001mm difference in X&Y compared to Z resolution of the scan.

[ref_rigid_def_data, ref_rigid_def_info, ~] = Get_reggui_data(handles, ref_rigid_def);

figure(9978)
histogram(ref_rigid_def_data,255);
xlabel('HU');
ylabel('# of voxels');
title('HU histogram reference image, after registration');
set(gca,'YScale','log')


% Apply Point Spread Function measured on CT scanner to reference CEM
sigma = PSF_FWHM/2.355; % mm
%sigma_pxl = ceil([sigma/pxlSizeXY, sigma/pxlSizeXY, sigma/pxlSizeZ]); % in pxl
sigma_pxl = sigma/pxlSizeXY;
fprintf('Convolving registered reference CEM by PSF \n');
ref_rigid_def_data = imgaussfilt3(ref_rigid_def_data, sigma_pxl, "Padding","circular","FilterDomain","auto");
handles = Set_reggui_data(handles, ref_rigid_def, ref_rigid_def_data, ref_rigid_def_info, 'images', 1);


% Save registered images to disk
handles = save2Disk(handles, ref_rigid_def_data, size(ref_rigid_def_data), ref_rigid_def_info, ref_rigid_def, fullfile(handles.dataPath));
[ref_rigid_trans_data, ref_rigid_trans_info, ~] = Get_reggui_data(handles, ref_rigid_trans);
handles = save2Disk(handles, ref_rigid_trans_data, size(ref_rigid_trans_data), ref_rigid_trans_info, ref_rigid_trans, fullfile(handles.dataPath));


% Autothresold normalised image for better spatial dimensional preservation
fprintf('Computing scan and reference binary images \n');
scan_cef_mask = 'scan_cef_mask';
%handles = AutoThreshold(scan_CEF_imageName, [128], scan_cef_mask, handles);
handles = ManualThreshold(scan_CEF_imageName, [scan_mask_low_thres scan_mask_high_thres], cef_mask, handles);
ref_cef_mask = 'ref_cef_mask';
%handles = AutoThreshold(ref_rigid_def, [128], ref_cef_mask, handles);
handles = ManualThreshold(ref_rigid_def, [-220 scan_mask_high_thres], ref_cef_mask, handles);



% Save cropped scan and ref masks to disk
[scan_cef_mask_data, scan_CEF_mask_info, ~] = Get_reggui_data(handles, scan_cef_mask);
handles = save2Disk(handles, scan_cef_mask_data, size(scan_cef_mask_data), scan_CEF_mask_info, scan_cef_mask, fullfile(handles.dataPath));
[ref_cef_mask_data, ref_CEF_mask_info, ~] = Get_reggui_data(handles, ref_cef_mask);
handles = save2Disk(handles, ref_cef_mask_data, size(ref_cef_mask_data), ref_CEF_mask_info, ref_cef_mask, fullfile(handles.dataPath));



% Compute binary difference between masks
% [scan_cef_contour,~,~] = Get_reggui_data(handles, scan_cef_mask,'images');
% [ref_cef_contour,~,~] = Get_reggui_data(handles, ref_cef_mask,'images');
% mask_difference_3D = scan_cef_contour - ref_cef_contour;
mask_difference_3D = scan_cef_mask_data - ref_cef_mask_data;


% Scan geometrical dimensions
fprintf('Computing geometrical differences (scan - ref) \n');
scan_distance_map2D = squeeze(sum(scan_cef_mask_data, 3)).* pxlSizeZ;

figure(998)
contourf(scan_distance_map2D, 100, 'LineColor', 'none');
colorbar;
%caxis([-40 40]);
xlabel('Y (mm)');
ylabel('X (mm)');
title('Scan geometrical dimensions in Z direction (mm)');
xticklabels(xticks * pxlSizeXY)
yticklabels(yticks * pxlSizeXY)



% Reference geometrical dimensions
ref_distance_map2D = squeeze(sum(ref_cef_mask_data, 3)).* pxlSizeZ;

figure(999)
contourf(ref_distance_map2D, 100, 'LineColor', 'none');
colorbar;
%caxis([-40 40]);
xlabel('Y (mm)');
ylabel('X (mm)');
title('Ref geometrical dimensions in Z direction (mm)');
xticklabels(xticks * pxlSizeXY)
yticklabels(yticks * pxlSizeXY)



% Get mean and std of geometrical differences of scan - ref
dist_diff2D_map = squeeze(sum(mask_difference_3D, 3)).* pxlSizeZ; % makes more sense to compute it on the cumulated distance than per voxel as per previous line.
mean_dist_diff_scan_ref = mean(dist_diff2D_map, 'all', 'omitnan');
std_dist_diff_scan_ref = std(dist_diff2D_map, [], 'all', 'omitnan');
min_dist_diff_scan_ref = min(dist_diff2D_map, [], 'all', 'omitnan');
max_dist_diff_scan_ref = max(dist_diff2D_map, [], 'all', 'omitnan');
fprintf('Mean distance difference scan-ref (mm): %d. \n', mean_dist_diff_scan_ref);
fprintf('STD distance difference scan-ref (mm): %d. \n', std_dist_diff_scan_ref);
fprintf('Min distance difference scan-ref (mm): %d. \n', min_dist_diff_scan_ref);
fprintf('Max distance difference scan-ref (mm): %d. \n', max_dist_diff_scan_ref);

fileID = fopen(fullfile(outputPath, 'summary.txt'),'w');
fprintf(fileID,'Mean distance difference scan-ref (mm): %d. \n', mean_dist_diff_scan_ref);
fprintf(fileID,'STD distance difference scan-ref (mm): %d. \n', std_dist_diff_scan_ref);
fprintf(fileID,'Min distance difference scan-ref (mm): %d. \n', min_dist_diff_scan_ref);
fprintf(fileID,'Max distance difference scan-ref (mm): %d. \n', max_dist_diff_scan_ref);



% Compute geometrical differences
distance_diff_map2D = squeeze(sum(mask_difference_3D, 3)).* pxlSizeZ;

figure(1000)
%axx = gca;
[~,cc] = contourf(distance_diff_map2D, 100, 'LineColor', 'none');
colorbar;
%xlabel('Y (# voxels)');
%ylabel('X (# voxels)');
title('Map of distance differences in Z direction (scan - ref)(mm)');
xlabel('Y (mm)');
ylabel('X (mm)');
xticklabels(xticks .* pxlSizeXY)
yticklabels(yticks .* pxlSizeXY)



% Save binary mask difference = scan - ref
difference3D_mask_name = 'scan-ref_cef_mask_diff';
difference3D_CT_info = handles.images.info{2}; % copy CT info from CEF if existent
handles = Set_reggui_data(handles, difference3D_mask_name, mask_difference_3D, difference3D_CT_info, 'images', 1); %Add the interpolated CT scan to the list of CT scan
handles = save2Disk(handles, mask_difference_3D, size(mask_difference_3D), difference3D_CT_info, difference3D_mask_name, fullfile(handles.dataPath));
%--------------------------------------------------------------------------




% Gamma analysis distance difference
% #########################################################################
if(geomGammaIndexAnalysis)
    options.dd = config.flags.geom_dd; % (mm)
    options.DD = config.flags.geom_DD; % (%)
    options.FI = 5; % number of points for internal interpolation
    options.global_ref = 1; % local (0) or global (1 - default) WET difference will be used
    options.threshold = 5; % (% unit) specify a WET threshold (in % of the reference WET) under which the gamma is not computed
    options.search_distance = 5; % specify a reference WET to be used for global WET computation and thresholding instead of the max WET
    fprintf('Computing distance gamma index (ref - scan) \n');
    [gamma, myMask, myIm1_resampled] = gamma_2D(ref_distance_map2D, ref_CEF_mask_info, scan_distance_map2D, scan_CEF_mask_info, [], options);
    
    
    figure(10003)
    contourf(gamma,100,'LineColor','none');
    colorbar;
    clim([0 5]); % limited to meaningful gamma index of 2. The actual values may go higher in certain areas...
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title('Gamma index (ref-scan) Z direction');
    xticklabels(xticks * pxlSizeXY)
    yticklabels(yticks * pxlSizeXY)
    
    
    % Compute passing rate
    gamma = single(gamma);
    passing_rate = length(gamma(gamma<1 & myMask>=0.5))/length(gamma(myMask>=0.5)); % global passing rate (union of masks + threshold)
    average_gamma = mean(gamma(myMask>=0.5)); % average gamma (union of masks + threshold)
    disp(['Average gamma: ',num2str(round(average_gamma,2))])
    disp(['General passing rate: ',num2str(round(passing_rate*100,2)),'%'])
end
% #########################################################################




% HU distribution analysis
% #########################################################################

% density_CEF_scan = hu_to_density(scan_cef_image, 'reggui_material_pmma_phantom.txt');
% denstiy_scan_2Dmap = sum(density_CEF_scan,2,'omitnan').* handles.spacing(2);
% mean_density_scan = mean(density_CEF_scan, 'all', 'omitnan');
% std_density_scan = std(density_CEF_scan, [], 'all', 'omitnan');
% fprintf('Mean density scanned CEF is equal to %d. \n', mean_density_scan);
% fprintf('STD density scanned CEF is equal to %d. \n', std_density_scan);


% Print density analysis values to screen
fprintf('Mean HU resampled/eroded scanned CEF is equal to %d. \n', mean_HU_scan);
fprintf('STD HU resampled/eroded scanned CEF is equal to %d. \n', std_HU_scan);
fprintf('Min HU resampled/eroded scanned CEF is equal to %d. \n', min_HU_scan);
fprintf('Max HU resampled/eroded scanned CEF is equal to %d. \n', max_HU_scan);

fprintf(fileID,'\nMean HU resampled/eroded scanned CEF is equal to %d. \n', mean_HU_scan);
fprintf(fileID,'STD HU resampled/eroded scanned CEF is equal to %d. \n', std_HU_scan);
fprintf(fileID,'Min HU resampled/eroded scanned CEF is equal to %d. \n', min_HU_scan);
fprintf(fileID,'Max HU resampled/eroded scanned CEF is equal to %d. \n', max_HU_scan);
fclose(fileID);


% Threshold scan image to remove noise, artefacts, outliers
fprintf('Thresholding out noise/artifacts/outliers \n');
%scan_cef_data(scan_cef_data > min_HU_scan & scan_cef_data < (mean_HU_scan-3*std_HU_scan)) = mean_HU_scan;
%scan_cef_data(scan_cef_data > (mean_HU_scan+3*std_HU_scan)) = mean_HU_scan;
%handles = Set_reggui_data(handles, scan_CEF_imageName, scan_cef_data, scan_cef_info,'images',1);


nr_slices = size(scan_cef_data, 3);
offset = 60; %200
slice = offset + floor(nr_slices/10);
slice_increment = ceil((nr_slices - slice)./9);


figure(1003);
tiledlayout(2,2);
for idx = 1:4
    nexttile
    contourf(squeeze(scan_cef_data(:,:,slice)),100,'LineColor','none');
    %contourf(squeeze(scan_cef_data(:,slice,:)),100,'LineColor','none');
    colorbar;
    %caxis([(mean_HU_scan - 3*std_HU_scan) (mean_HU_scan + 3*std_HU_scan)]);
    xticklabels(xticks * pxlSizeXY)
    yticklabels(yticks * pxlSizeXY)
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title(['Z plane ' num2str(slice)]);
    slice = slice+slice_increment;

end

sgtitle('Original scan HU distribution in Z plane')



% plot volume histogram
figure(1007);
histogram(reshape(scan_cef_data,size(scan_cef_data,1)*size(scan_cef_data,2)*size(scan_cef_data,3),1,1),255);
xlabel('HU');
ylabel('# of voxels');
title('Original scan HU histogram');
set(gca,'YScale','log')
%--------------------------------------------------------------------------



% Threshold scan image to remove noise, artefacts, outliers
fprintf('Thresholding in noise/artefacts/outliers \n');
outlier_scan_CEF_imageName = 'outliers_scan_CEF';
outlier_scan_cef_data = scan_cef_data;
outlier_scan_cef_data(outlier_scan_cef_data < -900) = NaN; % exclude air
outlier_scan_cef_data(outlier_scan_cef_data > (mean_HU_scan-3*std_HU_scan) & scan_cef_data < (mean_HU_scan+3*std_HU_scan)) = NaN; % exclude CEM densities
handles = Set_reggui_data(handles, outlier_scan_CEF_imageName, outlier_scan_cef_data, scan_cef_info,'images',1);


nr_slices = size(scan_cef_data, 3);
offset = 130; %200
slice = offset + floor(nr_slices/10);
slice_increment = ceil((nr_slices - slice)./9);


figure(1004);
tiledlayout(2,2);
for idx = 1:4
    nexttile
    contourf(squeeze(outlier_scan_cef_data(:,:,slice)),100,'LineColor','none');
    %contourf(squeeze(outlier_scan_cef_data(:,slice,:)),100,'LineColor','none');
    colorbar;
    %caxis([(mean_HU_scan - 3*std_HU_scan) (mean_HU_scan + 3*std_HU_scan)]);
    xticklabels(xticks * pxlSizeXY)
    yticklabels(yticks * pxlSizeXY)
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title(['Z plane ' num2str(slice)]);
    slice = slice+slice_increment;
end

sgtitle('Outliers/noise original scan HU distribution in Z plane')



% plot volume histogram
figure(1008);
histogram(reshape(outlier_scan_cef_data,size(outlier_scan_cef_data,1)*size(outlier_scan_cef_data,2)*size(outlier_scan_cef_data,3),1,1),255);
xlabel('HU');
ylabel('# of voxels');
title('Outliers/noise scan HU histogram');
set(gca,'YScale','log')
% #########################################################################




% WEPL Gamma index analysis
% #########################################################################
if (weplGammaIndexAnalysis)
    [def_ref_cef_data, def_ref_cef_info, ~] = Get_reggui_data(handles, ref_rigid_def);
    %def_ref_cef_data(def_ref_cef_data < mean_HU_scan) = mean_HU_scan;
    %def_ref_cef_data(def_ref_cef_data > max_HU_scan) = max_HU_scan;
    %handles = Set_reggui_data(handles, ref_CEF_imageName, def_ref_cef_data, def_ref_cef_info,'images',1);
    
    
    % Compute WEPLs
    fprintf('Computing scan/reference/difference WET maps \n');
    WEPL_Z_CEF_scan = hu_to_we(scan_cef_data, 'reggui_material_pmma_phantom.txt');
    WEPL_scan_2Dmap = sum(WEPL_Z_CEF_scan, 3, 'omitnan').* pxlSizeZ;
    
    handles = save2Disk(handles, WEPL_Z_CEF_scan, size(WEPL_Z_CEF_scan), difference3D_CT_info, 'scan_CEF_WEPL_3D', fullfile(handles.dataPath));
    
    
    WEPL_Z_CEF_ref = hu_to_we(def_ref_cef_data, 'reggui_material_pmma_phantom.txt');
    WEPL_ref_2Dmap = sum(WEPL_Z_CEF_ref, 3, 'omitnan').* pxlSizeZ;
    
    handles = save2Disk(handles, WEPL_Z_CEF_ref, size(WEPL_Z_CEF_ref), difference3D_CT_info, 'ref_CEF_WEPL_3D', fullfile(handles.dataPath));
    
    WEPL_diff_3Dmap = (WEPL_Z_CEF_scan - WEPL_Z_CEF_ref);
    WEPL_diff_2Dmap = WEPL_scan_2Dmap - WEPL_ref_2Dmap;
    
    handles = save2Disk(handles, WEPL_diff_3Dmap, size(WEPL_diff_3Dmap), difference3D_CT_info, 'scan_ref_CEF_WEPL_diff_3D', fullfile(handles.dataPath));
    
    
    %WEPL_diff_2Dmap = sum(WEPL_diff, 3, 'omitnan').* handles.spacing(2);
    figure(1005);
    contourf(squeeze(WEPL_diff_2Dmap),100,'LineColor','none');
    colorbar;
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title('WET difference Z direction (scan - ref)(mm)');
    xticklabels(xticks * pxlSizeXY)
    yticklabels(yticks * pxlSizeXY)
    
    %WEPL_scan_2Dmap = sum(WEPL_Z_CEF_scan, 3, 'omitnan').* handles.spacing(2);
    figure(1009);
    contourf(squeeze(WEPL_scan_2Dmap),100,'LineColor','none');
    colorbar;
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title('WET scan Z direction (mm)');
    xticklabels(xticks * pxlSizeXY)
    yticklabels(yticks * pxlSizeXY)

    
    %WEPL_opt_2Dmap = sum(WEPL_Z_CEF_ref, 3, 'omitnan').* handles.spacing(2);
    figure(1011);
    contourf(squeeze(WEPL_ref_2Dmap),100,'LineColor','none');
    colorbar;
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title('WET ref Z direction (mm)');
    xticklabels(xticks * pxlSizeXY)
    yticklabels(yticks * pxlSizeXY)


    % Compute gamma index in between the scan and the ref WEPL maps
    options.dd = config.flags.wepl_dd; % (mm)
    options.DD = config.flags.wepl_DD; % (%)
    options.FI = 5; % number of points for internal interpolation
    options.global_ref = 1; % local (0) or global (1 - default) WET difference will be used
    options.threshold = 5; % (% unit) specify a WET threshold (in % of the reference WET) under which the gamma is not computed
    options.search_distance = 5; % specify a reference WET to be used for global WET computation and thresholding instead of the max WET
    fprintf('Computing WET gamma index (ref - scan) WEPL \n');
    [gamma, myMask, myIm1_resampled] = gamma_2D(WEPL_ref_2Dmap,difference3D_CT_info,WEPL_scan_2Dmap,difference3D_CT_info,[],options);
    
    
    figure(1113)
    contourf(gamma,100,'LineColor','none');
    colorbar;
    clim([0 2]); % limited to meaningful gamma index of 2. The actual values may go higher in certain areas...
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title('Gamma index (ref-scan WEPL) Z direction');
    xticklabels(xticks * pxlSizeXY)
    yticklabels(yticks * pxlSizeXY)
    
    
    % Compute passing rate
    gamma = single(gamma);
    passing_rate = length(gamma(gamma<1 & myMask>=0.5))/length(gamma(myMask>=0.5)); % global passing rate (union of masks + threshold)
    average_gamma = mean(gamma(myMask>=0.5)); % average gamma (union of masks + threshold)
    disp(['Average gamma: ',num2str(round(average_gamma,2))])
    disp(['General passing rate: ',num2str(round(passing_rate*100,2)),'%'])
end

% Start reggui with script handles
%start_reggui_GUI(handles);

toc
end
% =========================================================================




% Loads the datasets from disk
% #########################################################################
function pre_handles = loadAlldatasets(scanCEF_path, CT_air_path, RSplanFileName, outputPath, Plan, permOrder, flipAxis, source)

    % Initialise reggui pre-handles
    pre_handles = Initialize_reggui_handles();
    pre_handles.dataPath = outputPath;
            
    switch source
        case 'CT'
            ref_CEF_imageName = 'ref_CEF';
            [ref_CEF_directory, ref_CEF_fileName, EXT] = fileparts(refCEF_path);
            pre_handles = Import_data(ref_CEF_directory, [ref_CEF_fileName EXT], 1, ref_CEF_imageName, pre_handles);
            [pre_ref_cef_data, pre_ref_cef_info, ~] = Get_reggui_data(pre_handles, ref_CEF_imageName);

        case 'STL'
            air_CT_imageName = 'air_CT';
            [air_CT_directory, air_CT_fileName, EXT] = fileparts(CT_air_path);
            pre_handles = Import_image(air_CT_directory, [air_CT_fileName EXT], 1, air_CT_imageName, pre_handles);
            [pre_air_CT_data, pre_air_CT_info, ~] = Get_reggui_data(pre_handles, air_CT_imageName);
            
            fprintf('Read STL file and generate 3D mask...\n');
            pre_handles.spacing = [0.5, 0.5, 0.5];
            stlObject = stlread(stl_path);
            CEMmask3D = stl2mask(stlObject, pre_handles.spacing);
            
            Plan.Spike.MaterialID = 'Accura_ClearVue';
            Plan.ScannerDirectory = 'default';
            
            % Export CEF in CT
            HU_air =  getMaterialPropCT('Schneider_Air', Plan.ScannerDirectory); %Hounsfield unit associated to air in the material file
            HU_CEF = getMaterialPropCT(Plan.Spike.MaterialID, Plan.ScannerDirectory); %HU of the CEF
            %CEF = double(CEMmask3D);
            CEF = single(CEMmask3D);
            CEF(CEF == 1) = HU_CEF;
            CEF(CEF == 0) = HU_air;
            
            % Instead of saving to disk you can add the new image to handles
            % with set_reggui_data and carry on as if imported from file
            %pre_handles = Set_reggui_data(pre_handles, 'CEF_in_CT', CEF, pre_air_CT_info, 'images', 1);
            handles2 = save2Disk(pre_handles, CEF, size(CEF), pre_air_CT_info, 'CEF_in_CT', outputPath);
            CEF = []; %free memory

        case 'Plan'
            air_CT_imageName = 'air_CT';
            [air_CT_directory, air_CT_fileName, EXT] = fileparts(CT_air_path);
            pre_handles = Import_image(air_CT_directory, [air_CT_fileName EXT], 1, air_CT_imageName, pre_handles);
            [pre_air_CT_data, pre_air_CT_info, ~] = Get_reggui_data(pre_handles, air_CT_imageName);
            
            Plan.Spike.MaterialID = 'Accura_ClearVue';
            Plan.showGraph = 'false';
            Plan.CTname = air_CT_imageName;
            [handles, Plan] = parseFLASHplan(RSplanFileName , Plan , pre_handles);
            
            pre_handles.origin(2) = -Plan.Beams.RangeModulator.IsocenterToRangeModulatorDistance; %Move CT scan to the position of the CEM, to save memory space
            PixelSize = Plan.Beams.RangeModulator.Modulator3DPixelSpacing;
            
            Plan.Beams.GantryAngle = 0;
            Plan.Beams.Beam.PatientSupportAngle =0;
            minField = Plan.Beams.RangeModulator.ModulatorOrigin;
            maxField = Plan.Beams.RangeModulator.ModulatorOrigin + Plan.Beams.RangeModulator.Modulator3DPixelSpacing .* size(Plan.Beams.RangeModulator.CEM3Dmask);
            Zdistal = Plan.Beams.RangeModulator.IsocenterToRangeModulatorDistance - Plan.Beams.RangeModulator.Modulator3DPixelSpacing(3) .* size(Plan.Beams.RangeModulator.CEM3Dmask,3) - 5;
            
            %Create a high res CT scan
            pre_handles  = createHighResCT(pre_handles , air_CT_imageName , air_CT_imageName , Plan.Beams , Plan.Beams.RangeModulator.Modulator3DPixelSpacing , -1050 , minField , maxField , Zdistal , pre_air_CT_info);
            [Plan, pre_handles] = setCEMinCT(pre_handles, Plan, air_CT_imageName);
            
            %Rotate CT so that Z axis is paralell to spikes, i.e. paralell to Zg
            pre_handles.images.data{2} = permute(pre_handles.images.data{2},[1,3,2]);
            %pre_handles.images.data{2} = flip(pre_handles.images.data{2}, 1);
            pre_handles.spacing = pre_handles.spacing([1,3,2]);
            pre_handles.origin = pre_handles.origin([1,3,2]);
            pre_handles.size = size(pre_handles.images.data{2});
            pre_handles.images.info{2}.Spacing = pre_handles.spacing;
            pre_handles.images.info{2}.ImagePositionPatient = pre_handles.origin;
                  
            [pre_ref_cef_data, pre_ref_cef_data_info, ~] = Get_reggui_data(pre_handles, air_CT_imageName);
            pre_handles.origin(3) = -round(pre_handles.size(3) ./2); %Move origin back to middle of CEM
            pre_ref_cef_data_info.PatientOrientation = pre_air_CT_info.PatientOrientation;
            pre_handles = Set_reggui_data(pre_handles, air_CT_imageName, pre_ref_cef_data, pre_ref_cef_data_info, 'mydata');
            
            % Save reference CEM image to disc
            pre_handles = save2Disk(pre_handles, pre_ref_cef_data, size(pre_ref_cef_data), pre_ref_cef_data_info, 'ref_CEF_image', outputPath);
            
            % Clear the images from handles as the reference needs to be 
            % imported as data (not image) for rigid registration
            pre_handles = Remove_all_images(pre_handles,1);
    end
    
    % Import scanned CEF 
    scan_CEF_imageName = 'scan_CEF';
    [scan_CEF_directory, scan_CEF_fileName, EXT] = fileparts(scanCEF_path);
    pre_handles = Import_image(scan_CEF_directory, [scan_CEF_fileName EXT],1,scan_CEF_imageName,pre_handles); %Load CEM image into handles.data
    [pre_scan_cef_data, pre_scan_cef_info, ~] = Get_reggui_data(pre_handles, scan_CEF_imageName);

    % Rotate the measured CT scan to place it in same orientation as the plan CEM
    % Permute the dimension of CEM image and flip some dimension to get the same orientation for the reference CT and the mesurement CT
    if ~isempty(flipAxis)
        for i = 1:size(flipAxis, 1)
            pre_handles.images.data{2} = flip(pre_handles.images.data{2}, flipAxis(i,:));
        end
    end
    if ~isempty(permOrder)
        for i = 1:size(permOrder, 1)
            pre_handles.images.data{2} = permute(pre_handles.images.data{2} , permOrder(i,:));
            pre_handles.images.info{2}.Spacing = pre_handles.images.info{2}.Spacing(permOrder(i,:));
        end
    end
    spc = pre_handles.images.info{2}.Spacing;
    pre_handles.images.info{2}.ImagePositionPatient = (- round(size(pre_handles.images.data{2}) ./2) .* spc')'; %Move origin back to middle of CEM 
end
% =========================================================================



function [crop_image, crop_image_info, handles] = cropDData(image_name, scan_mask_low_thres, scan_mask_high_thres, handles)

    image_mask = [image_name '_mask'];
    %handles = AutoThreshold(image_name, [128], image_mask, handles);
    handles = ManualThreshold(image_name, [scan_mask_low_thres scan_mask_high_thres], image_mask, handles);
    handles = Resample_all(handles, image_mask, [0], handles.spacing,'from_mask');
    [crop_image, crop_image_info, ~] = Get_reggui_data(handles, image_name);
end
% =========================================================================