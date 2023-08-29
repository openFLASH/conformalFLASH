% Contributors
% Authors : :Lucian Hotoiu (open.reggui@gmail.com)
% -------------------------------------------------------------------------

tic
clearvars;
close all;

% !!!!!!!! README !!!!!!!! ------------------------------------------------
% Prior to running the script rotate and crop the two images to be compared.
% The object in the two images should have roughly the same orientation,
% position and occupy similar areas in the image.
% -------------------------------------------------------------------------

%Inputs / Outputs
scanCEF_path = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\CT_CEM\reggui_20230201_zz_CEF\reggui_20230201_zz_CEF_0001.dcm';
permOrder = [1,3,2]; %Permutation of the dimensiosn of CT scan to align it wiht plan
flipAxis = 3; %Which axis index should be flipped (after permute)

%refCEF_path = 'E:\Lucian\Data\ScanCEFUpenn_sylvain\ctCEF505030_DCM\resample_333um\CEF_opt_resample.mhd';
%stl_path = 'C:\Users\lhotoiu\Downloads\CEF1id_square_decimated_more.stl';

CT_air_path = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\reggui_CT_air\reggui_CT_air_0001.dcm';

RSplanFileName = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\D-58\RP_D58.dcm'
Plan.BDL = 'D:\programs\openREGGUI\flash\openMCsquare\lib\BDL\BDL_default_UN1_G0_Al_RangeShifter_tilted.txt';
Plan.ScannerDirectory = 'D:\programs\openREGGUI\REGGUI\plugins\openMCsquare\lib\Scanners\default';

%outputPath = 'E:\Lucian\Data\ScanCEFUpenn_sylvain\OUTPUT2_RG35';
outputPath = 'D:\programs\openREGGUI\REGGUI_userdata\raystation\CEM_compare';

% -------------------------------------------------------------------------


%Load the data
pre_handles = loadAlldatasets(scanCEF_path , CT_air_path , RSplanFileName , outputPath , Plan , permOrder , flipAxis);
air_CT_imageName = 'air_CT';
scan_CEF_imageName = 'scan_CEF';
[pre_ref_cef_data, pre_ref_cef_info, ~] = Get_reggui_data(pre_handles, air_CT_imageName);
[pre_scan_cef_data, pre_scan_cef_info, ~] = Get_reggui_data(pre_handles, scan_CEF_imageName);

% Compute grid for resampling
fprintf('Computing relevant region for probing and assesing scan object HU density \n');
grid.origin = [1; 1; 1];

grid.size = floor((size(pre_scan_cef_data)' ./ (pre_ref_cef_info.Spacing ./ pre_scan_cef_info.Spacing)));
grid.spacing = pre_ref_cef_info.Spacing;

pre_handles = Resample_all(pre_handles, grid.origin, grid.size, grid.spacing);
[pre_scan_cef_data, pre_scan_cef_info, ~] = Get_reggui_data(pre_handles, scan_CEF_imageName);
%pre_scan_cef_data = medfilt3(pre_scan_cef_data, [7 7 7], 'replicate'); % very slow

% Threshold, erode and intersect
cef_mask = 'scan_cef_mask';
pre_handles = AutoThreshold(scan_CEF_imageName, [128], cef_mask, pre_handles);
pre_handles = Erosion(cef_mask, [4 4 4], 'eroded_cef_mask', pre_handles);
[eroded_cef_mask_data, eroded_cef_mask_info, ~] = Get_reggui_data(pre_handles, 'eroded_cef_mask');
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
offset = 20; % 50
slice = offset + floor(nr_slices/10);
slice_increment = 2; %1

figure(903);
tiledlayout(3,3);

for idx = 1:9
    nexttile
    contourf(squeeze(intersect_data(:,:,slice)),100,'LineColor','none');
    %contourf(squeeze(intersect_data(:,slice,:)),100,'LineColor','none');
    colorbar;
    caxis([(mean_HU_scan - 3*std_HU_scan) (mean_HU_scan + 3*std_HU_scan)]);
    xticklabels(xticks * pre_handles.spacing(3))
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



% -------------------------------------------------------------------------


%Load again the data
%as the previous data set was resampled in pre_handles
handles = loadAlldatasets(scanCEF_path , CT_air_path , RSplanFileName , outputPath , Plan , permOrder , flipAxis);

% Import scanned CEF
[scan_cef_data, scan_cef_info, ~] = Get_reggui_data(handles, scan_CEF_imageName);


% Normalise scan CEF hi-res image
% min_intensity_scan_cef = min(min(min(scan_cef_data, [], 'omitnan')));
% max_intensity_scan_cef = max(max(max(scan_cef_data, [], 'omitnan')));
% norm_scan_scan_cef_data = (scan_cef_data - min_intensity_scan_cef)./ (max_intensity_scan_cef - min_intensity_scan_cef);
%
% handles = Set_reggui_data(handles, scan_CEF_imageName, norm_scan_scan_cef_data, scan_CEF_info, 'images', 1);



% % Import opt CEF
ref_CEF_imageName = 'air_CT';
[ref_cef_data, ref_cef_info, ~] = Get_reggui_data(handles, ref_CEF_imageName);


% Normalise reference CEF low-res image
% min_intensity_ref_cef = min(min(min(ref_cef_data)));
% max_intensity_ref_cef = max(max(max(ref_cef_data)));
% norm_ref_cef_data = (ref_cef_data - min_intensity_ref_cef)./ (max_intensity_ref_cef - min_intensity_ref_cef);
%
% handles = Set_reggui_data(handles, ref_CEF_imageName, norm_ref_cef_data, ref_CEF_info, 'images', 1);



% Rigid registration normalised scan and reference cef
fprintf('Computing rigid registration between scan (fixed) and reference (moving) images \n');
ref_rigid_def = 'ref_rigid_def';
ref_rigid_trans = 'ref_rigid_trans';
handles = Registration_ITK_rigid_multimodal(scan_CEF_imageName, ref_CEF_imageName, ref_rigid_def, ref_rigid_trans,handles);


% Save registered images to disk
[ref_rigid_def_data, ref_rigid_def_info, ~] = Get_reggui_data(handles, ref_rigid_def);
handles = save2Disk(handles, ref_rigid_def_data, size(ref_rigid_def_data), ref_rigid_def_info, ref_rigid_def, fullfile(handles.dataPath));
[ref_rigid_trans_data, ref_rigid_trans_info, ~] = Get_reggui_data(handles, ref_rigid_trans);
ref_rigid_trans_info.PatientOrientation = 'HFS';
handles = save2Disk(handles, ref_rigid_trans_data, size(ref_rigid_trans_data), ref_rigid_trans_info, ref_rigid_trans, fullfile(handles.dataPath));


% Autothresold normalised image for better spatial dimensional preservation
fprintf('Computing scan and reference binary images \n');
scan_cef_mask = 'scan_cef_mask';
handles = AutoThreshold(scan_CEF_imageName, [128], scan_cef_mask, handles);
ref_cef_mask = 'ref_cef_mask';
handles = AutoThreshold(ref_rigid_def, [128], ref_cef_mask, handles);



%Autothreshold does not work in scripting in this particular case for the
%registered ref cef
% Load it from disk instead
% defCEFmask_path = 'C:\Users\lhotoiu\Downloads\reggui_def_cef_mask\reggui_def_cef_mask_0001.dcm';
% [def_CEFmask_directory, def_CEFmask_fileName, EXT] = fileparts(defCEFmask_path);
% handles = Import_data(def_CEFmask_directory, [def_CEFmask_fileName EXT], 1, ref_cef_mask, handles);




% Save cropped scan and ref masks to disk
[scan_cef_mask_data, scan_CEF_mask_info, ~] = Get_reggui_data(handles, scan_cef_mask);
handles = save2Disk(handles, scan_cef_mask_data, size(scan_cef_mask_data), scan_CEF_mask_info, scan_cef_mask, fullfile(handles.dataPath));
[ref_cef_mask_data, ref_CEF_mask_info, ~] = Get_reggui_data(handles, ref_cef_mask);
handles = save2Disk(handles, ref_cef_mask_data, size(ref_cef_mask_data), ref_CEF_mask_info, ref_cef_mask, fullfile(handles.dataPath));


% Compute distance distribution between the optimised CEF and the scanned CEF
% [indicators, distances1D, distances3D] = Contour_distance(scan_cef_mask, ref_cef_mask, handles); %0.333 (limit for which distance calculation still works) is the isotropic size of the mesh to compute the distances


% Crop sides of scan image to remove artefacts
%scan_cef_data_crop = imcrop3(scan_cef_data, [1 1 1 500 570 700]); % the cuboid inverts x and y coordinates of the initial image to crop
%scan_cef_data_crop = imcrop3(scan_cef_data_crop, [20 5 20 110 180 190]); % the cuboid inverts x and y coordinates of the initial image to crop
%handles = Set_reggui_data(handles, scan_CEF_imageName, scan_cef_data_crop, scanCEF_info, 'images', 1);


% Resample opt image to same spacing as scan - nearest interpolation to avoid smoothing
% cur_sz = size(opt_cef_data);
% target_sz = cur_sz .* (ref_CEF_info.Spacing./scan_CEF_info.Spacing)';
% [X, Y, Z] = meshgrid(linspace(1,cur_sz(2),target_sz(2)),linspace(1,cur_sz(1),target_sz(1)),linspace(1,cur_sz(3),target_sz(3)));
% X = single(X);
% Y = single(Y);
% Z = single(Z);
% opt_cef_data_resample = interp3(opt_cef_data,X,Y,Z,'nearest');
% optCEF_info = scanCEF_info;
% %handles = Set_reggui_data(handles, 'opt_resampled', opt_cef_data_resample, optCEF_info, 'images', 1); %Add the interpolated CT scan to the list of CT scan
% handles = Set_reggui_data(handles, opt_CEF_imageName, opt_cef_data_resample, optCEF_info, 'images', 1); %Add the interpolated CT scan to the list of CT scan

% opt_cef_data_crop = imcrop3(opt_cef_data_resample, [1 1 1 500 570 700]); % the cuboid inverts x and y coordinates of the initial image to crop
% handles = Set_reggui_data(handles, opt_CEF_imageName, opt_cef_data_crop, optCEF_info, 'images', 1);



% Compute binary difference between masks
% [scan_cef_contour,~,~] = Get_reggui_data(handles, scan_cef_mask,'images');
% [ref_cef_contour,~,~] = Get_reggui_data(handles, ref_cef_mask,'images');
% mask_difference_3D = scan_cef_contour - ref_cef_contour;
mask_difference_3D = scan_cef_mask_data - ref_cef_mask_data;

pxlSizeZ = scan_cef_info.Spacing(3); % we approximate here that the mm are isotropical. In reality there is a 0.001mm difference in X&Y compared to Z resolution of the scan.
fprintf('Computing geometrical differences (scan - ref) \n');

% Scan geometrical dimensions
scan_distance_map2D = squeeze(sum(scan_cef_mask_data, 3)).* pxlSizeZ;

figure(998)
contourf(scan_distance_map2D, 100, 'LineColor', 'none');
colorbar;
%caxis([-40 40]);
xlabel('Y (mm)');
ylabel('X (mm)');
title('Scan geometrical dimensions in Z direction (mm)');
xticklabels(xticks * pxlSizeZ)
yticklabels(yticks * pxlSizeZ)


% Reference geometrical dimensions
ref_distance_map2D = squeeze(sum(ref_cef_mask_data, 3)).* pxlSizeZ;

figure(999)
contourf(ref_distance_map2D, 100, 'LineColor', 'none');
colorbar;
%caxis([-40 40]);
xlabel('Y (mm)');
ylabel('X (mm)');
title('Ref geometrical dimensions in Z direction (mm)');
xticklabels(xticks * pxlSizeZ)
yticklabels(yticks * pxlSizeZ)


% Get mean and std of geometrical differences of scan - ref
dist_diff3D = mask_difference_3D .* pxlSizeZ;
dist_diff2D = squeeze(sum(mask_difference_3D, 3)).* pxlSizeZ; % makes more sense to compute it on the cumulated distance than per voxel as per previous line.
mean_dist_diff_scan_ref = mean(dist_diff2D, 'all', 'omitnan');
std_dist_diff_scan_ref = std(dist_diff2D, [], 'all', 'omitnan');
min_dist_diff_scan_ref = min(dist_diff2D, [], 'all', 'omitnan');
max_dist_diff_scan_ref = max(dist_diff2D, [], 'all', 'omitnan');
fprintf('Mean distance difference scan-ref (mm): %d. \n', mean_dist_diff_scan_ref);
fprintf('STD distance difference scan-ref (mm): %d. \n', std_dist_diff_scan_ref);
fprintf('Min distance difference scan-ref (mm): %d. \n', min_dist_diff_scan_ref);
fprintf('Max distance difference scan-ref (mm): %d. \n', max_dist_diff_scan_ref);



% Compute geometrical differences
distance_diff_map2D = squeeze(sum(mask_difference_3D, 3)).* pxlSizeZ;

figure(1000)
contourf(distance_diff_map2D, 100, 'LineColor', 'none');
colorbar;
%caxis([-40 40]);
xlabel('Y (mm)');
ylabel('X (mm)');
title('Map of distance differences in Z direction (scan - ref)(mm)');
xticklabels(xticks * pxlSizeZ)
yticklabels(yticks * pxlSizeZ)

% Save binary mask difference = scan - ref
difference3D_mask_name = 'scan-ref_cef_mask_diff';
difference3D_CT_info = handles.images.info{2}; % copy CT info from CEF if existent
handles = Set_reggui_data(handles, difference3D_mask_name, mask_difference_3D, difference3D_CT_info, 'images', 1); %Add the interpolated CT scan to the list of CT scan
difference3D_CT_info.PatientOrientation = 'HFS';
handles = save2Disk(handles, mask_difference_3D, size(mask_difference_3D), difference3D_CT_info, difference3D_mask_name, fullfile(handles.dataPath));


% Save 3D distance between contours = scan - opt
difference3D_dist_heatmap_name = 'scan-ref_cef_dist_diff_map';
handles = Set_reggui_data(handles, difference3D_dist_heatmap_name, dist_diff3D, difference3D_CT_info, 'images', 1); %Add the interpolated CT scan to the list of CT scan
handles = save2Disk(handles, dist_diff3D, size(dist_diff3D), difference3D_CT_info, difference3D_dist_heatmap_name, fullfile(handles.dataPath));

%--------------------------------------------------------------------------



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
tiledlayout(3,3);
for idx = 1:9
    nexttile
    contourf(squeeze(scan_cef_data(:,:,slice)),100,'LineColor','none');
    %contourf(squeeze(scan_cef_data(:,slice,:)),100,'LineColor','none');
    colorbar;
    caxis([(mean_HU_scan - 3*std_HU_scan) (mean_HU_scan + 3*std_HU_scan)]);
    xticklabels(xticks * pxlSizeZ)
    yticklabels(yticks * pxlSizeZ)
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
offset = 60; %200
slice = offset + floor(nr_slices/10);
slice_increment = ceil((nr_slices - slice)./9);


figure(1004);
tiledlayout(3,3);
for idx = 1:9
    nexttile
    contourf(squeeze(outlier_scan_cef_data(:,:,slice)),100,'LineColor','none');
    %contourf(squeeze(outlier_scan_cef_data(:,slice,:)),100,'LineColor','none');
    colorbar;
    %caxis([(mean_HU_scan - 3*std_HU_scan) (mean_HU_scan + 3*std_HU_scan)]);
    xticklabels(xticks * pxlSizeZ)
    yticklabels(yticks * pxlSizeZ)
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



%--------------------------------------------------------------------------


[def_ref_cef_data, def_ref_cef_info, ~] = Get_reggui_data(handles, ref_rigid_def);
%def_ref_cef_data(def_ref_cef_data < mean_HU_scan) = mean_HU_scan;
%def_ref_cef_data(def_ref_cef_data > max_HU_scan) = max_HU_scan;
%handles = Set_reggui_data(handles, ref_CEF_imageName, def_ref_cef_data, def_ref_cef_info,'images',1);


% Compute WEPLs
fprintf('Computing scan/reference/difference WET maps \n');
WEPL_Z_CEF_scan = hu_to_we(scan_cef_data, 'reggui_material_pmma_phantom.txt');
WEPL_scan_2Dmap = sum(WEPL_Z_CEF_scan, 3, 'omitnan').* pxlSizeZ;

handles = save2Disk(handles, WEPL_scan_2Dmap, size(WEPL_scan_2Dmap), difference3D_CT_info, 'scan_CEF_WEPL_2D', fullfile(handles.dataPath));


WEPL_Z_CEF_ref = hu_to_we(def_ref_cef_data, 'reggui_material_pmma_phantom.txt');
WEPL_ref_2Dmap = sum(WEPL_Z_CEF_ref, 3, 'omitnan').* pxlSizeZ;

handles = save2Disk(handles, WEPL_ref_2Dmap, size(WEPL_ref_2Dmap), difference3D_CT_info, 'ref_CEF_WEPL_2D', fullfile(handles.dataPath));


WEPL_diff = (WEPL_Z_CEF_scan - WEPL_Z_CEF_ref).* pxlSizeZ;
WEPL_diff_2Dmap = WEPL_scan_2Dmap - WEPL_ref_2Dmap;

handles = save2Disk(handles, WEPL_diff, size(WEPL_diff), difference3D_CT_info, 'scan_ref_CEF_WEPL_diff_3D', fullfile(handles.dataPath));


%WEPL_diff_2Dmap = sum(WEPL_diff, 3, 'omitnan').* handles.spacing(2);
figure(1005);
contourf(squeeze(WEPL_diff_2Dmap),100,'LineColor','none');
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
title('WET difference Z direction (scan - ref)(mm)');
xticklabels(xticks * pxlSizeZ)
yticklabels(yticks * pxlSizeZ)

%WEPL_scan_2Dmap = sum(WEPL_Z_CEF_scan, 3, 'omitnan').* handles.spacing(2);
figure(1009);
contourf(squeeze(WEPL_scan_2Dmap),100,'LineColor','none');
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
title('WET scan Z direction (mm)');
xticklabels(xticks * pxlSizeZ)
yticklabels(yticks * pxlSizeZ)

%WEPL_opt_2Dmap = sum(WEPL_Z_CEF_ref, 3, 'omitnan').* handles.spacing(2);
figure(1011);
contourf(squeeze(WEPL_ref_2Dmap),100,'LineColor','none');
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
title('WET ref Z direction (mm)');
xticklabels(xticks * pxlSizeZ)
yticklabels(yticks * pxlSizeZ)

%--------------------------------------------------------------------------


% Compute gamma index in between the scan and the ref WEPL maps
%options.dd = 0.4; % (mm)
%options.DD = 1; % (%)
%options.FI = 5; % number of points for internal interpolation
%options.global_ref = 1; % local (0) or global (1 - default) WET difference will be used
%options.threshold = 5; % (% unit) specify a WET threshold (in % of the reference WET) under which the gamma is not computed
%options.search_distance = 5; % specify a reference WET to be used for global WET computation and thresholding instead of the max WET
%fprintf('Computing WET gamma index (ref - scan) \n');
%[gamma, myMask, myIm1_resampled] = gamma_2D(WEPL_ref_2Dmap,difference3D_CT_info,WEPL_scan_2Dmap,difference3D_CT_info,[],options);


% figure(1113)
% contourf(gamma,100,'LineColor','none');
% colorbar;
% caxis([0 2]); % limited to meaningul gamma index of 2. The actual values may go higher in certain areas...
% xlabel('X (mm)');
% ylabel('Y (mm)');
% title('Gamma index (ref-scan) Z direction');
% xticklabels(xticks * pxlSizeZ)
% yticklabels(yticks * pxlSizeZ)


% Compute passing rate
% gamma = single(gamma);
% passing_rate = length(gamma(gamma<1 & myMask>=0.5))/length(gamma(myMask>=0.5)); % global passing rate (union of masks + threshold)
% average_gamma = mean(gamma(myMask>=0.5)); % average gamma (union of masks + threshold)
% disp(['Average gamma: ',num2str(round(average_gamma,2))])
% disp(['General passing rate: ',num2str(round(passing_rate*100,2)),'%'])


% Start reggui with script handles
%start_reggui_GUI(handles);

toc


%---------------------------------------------------------
% Loads the datasets from disk
% Rotate the measured CT scan to place it in same orientation as the plan CEM
%---------------------------------------------------------
function pre_handles = loadAlldatasets(scanCEF_path , CT_air_path , RSplanFileName , outputPath , Plan , permOrder , flipAxis)

  % Initialise reggui pre-handles
  pre_handles = Initialize_reggui_handles();
  pre_handles.dataPath = outputPath;


  % Import reference CEF
  %source = 'CT';
  %source = 'STL';
  source = 'Plan';

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

          Plan.Spike.MaterialID = 'ABS_Resin';
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

          Plan.Spike.MaterialID = 'ABS_Resin';
          Plan.showGraph = 'false';
          Plan.CTname = air_CT_imageName;
          [handles, Plan] = parseFLASHplan(RSplanFileName , Plan , pre_handles);

          pre_handles.origin(2) = -Plan.Beams.RangeModulator.IsocenterToRangeModulatorDistance; %Move CT scan to the position of the CEM, to save memory space
          PixelSize = Plan.Beams.RangeModulator.Modulator3DPixelSpacing ;

          Plan.Beams.GantryAngle = 0;
          Plan.Beams.Beam.PatientSupportAngle =0;
          minField = Plan.Beams.RangeModulator.ModulatorOrigin;
          maxField = Plan.Beams.RangeModulator.ModulatorOrigin + Plan.Beams.RangeModulator.Modulator3DPixelSpacing .* size(Plan.Beams.RangeModulator.CEM3Dmask);
          Zdistal = Plan.Beams.RangeModulator.IsocenterToRangeModulatorDistance - Plan.Beams.RangeModulator.Modulator3DPixelSpacing(3) .* size(Plan.Beams.RangeModulator.CEM3Dmask,3) - 5;

          %Create a high res CT scan
          pre_handles  = createHighResCT(pre_handles , air_CT_imageName , air_CT_imageName , Plan.Beams , Plan.Beams.RangeModulator.Modulator3DPixelSpacing , -1050 , minField , maxField , Zdistal , pre_air_CT_info);
          [Plan, pre_handles] = setCEMinCT(pre_handles , Plan , air_CT_imageName );

          %Rotate CT so that Z axis is paralell to spikes, i.e. paralell to Zg
          pre_handles.images.data{2} = permute(pre_handles.images.data{2},[1,3,2]);
          pre_handles.spacing = pre_handles.spacing([1,3,2]);
          pre_handles.origin  = pre_handles.origin([1,3,2]);
          pre_handles.size = size(pre_handles.images.data{2});
          pre_handles.images.info{2}.Spacing = pre_handles.spacing;
          pre_handles.images.info{2}.ImagePositionPatient = pre_handles.origin;


          [pre_ref_cef_data, pre_ref_cef_info, ~] = Get_reggui_data(pre_handles, air_CT_imageName);
          pre_handles.origin(3) = -round(pre_handles.size(3) ./2); %Move origin back to middle of CEM

  %Display one slice for information
  minTot = min(pre_ref_cef_data,[],'all');
  maxTot = max(pre_ref_cef_data,[],'all');

   %for i = 1:size(pre_ref_cef_data,3)
    i = 50;
    figure(10)
    imagf = squeeze(pre_ref_cef_data(:,:,i));
    image((imagf - minTot) ./ (maxTot  - minTot) .* 255)
    title(['Reference -- Slice = ' num2str(i)])
     %pause
   %end


  end

  % Import scanned CEF to apply median filter and erode
  scan_CEF_imageName = 'scan_CEF';
  [scan_CEF_directory, scan_CEF_fileName, EXT] = fileparts(scanCEF_path);
  pre_handles = Import_data(scan_CEF_directory, [scan_CEF_fileName EXT],1,scan_CEF_imageName,pre_handles); %Load CEM image into handles.data

  %permute the dimension of CEM image and flip some dimension
  %in order to get the smae orientation for the reference CT and the meausrmenet CT
  if ~isempty(permOrder)
    pre_handles.mydata.data{2} = permute(pre_handles.mydata.data{2} , permOrder);
    pre_handles.mydata.info{2}.Spacing = pre_handles.mydata.info{2}.Spacing(permOrder);
  end
  if flipAxis
    pre_handles.mydata.data{2} = flipdim(pre_handles.mydata.data{2}, flipAxis);
  end
  spc = pre_handles.mydata.info{2}.Spacing;
  pre_handles.mydata.info{2}.ImagePositionPatient = - round(size(pre_handles.mydata.data{2}) ./2) .* spc'; %Move origin back to middle of CEM

  pre_handles = Data2image(scan_CEF_imageName,scan_CEF_imageName,pre_handles,0); %Load image in handles.data
  [pre_scan_cef_data, pre_scan_cef_info, ~] = Get_reggui_data(pre_handles, scan_CEF_imageName); %resample the image and move it from handles.data to hansdles.image

  %Display one slice for information
  minTot = min(pre_scan_cef_data,[],'all');
  maxTot = max(pre_scan_cef_data,[],'all');
  %for i = 1:size(pre_scan_cef_data,3)
    i=50;
    figure(11)
    imagf = squeeze(pre_scan_cef_data(:,:,i));
    image((imagf - minTot) ./ (maxTot  - minTot) .* 255)
    title(['Measurement -- Slice = ' num2str(i)])
    %pause
  %end

end
