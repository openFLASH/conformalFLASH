%% Pij2Dose
% Compute the 3D dose distribution from the influence matrix Pij and thwe weight |w| of the spots.
% The voxel of the |dose| are ordered along the DICOM patient CS.
%
% NOTE: In order to reduce memory usage, MIROpt has a list of "nominal voxels" and "robust voxels" where the dose influnece matrix is going to be stored. The corresponding variables are:
% Plan.OptROIVoxels_nominal and Plan.OptROIVoxels_robust
% These voxels include the structures selected for optimization in nominal and robust scenarios respectively,
% so the dose outside these regions is not stored because it will not be used for the objective function.
%That is why, at the end, we should perform a final dose calculation with MCsquare to have the final dose on the whole body. The filtering of the Influence matrix with these voxels is done in the ReadPijMatrix.

%
%% Syntax
% |dose = Pij2Dose(Pij , w , DoseGridSize)|
%
%
%% Description
% |dose = Pij2Dose(Pij , w , DoseGridSize)| Description
%
%
%% Input arguments
% |Pij| -_SCALAR VECTOR_- Dose influence matrix: |P(vox,spot)| The dose contribution to voxel |vox| of the spot number |spot|
%
% |w| - _SCALAR VECTOR_ - |w(i)| optimised weight of the i-th spot
%
% |DoseGridSize| -_SCALAR VECTOR_- [Nx, Ny, Nz] Number of pixels of the Dose map along each dimension
%
%% Output arguments
%
% |dose| - _SCALAR MATRIX_ -  |dose(x,y,z)| Dose at the voxel (x,y,z) in the DICOM patient CS
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function dose = Pij2Dose(Pij , w , DoseGridSize)

dose = zeros(DoseGridSize);

D = Pij * w' ; %dose contribution of all Bragg peaks contributing to the beamlet into the whole ROI
dose(1:prod(DoseGridSize))= D;

dose = flipdim(dose,3); %The Zaxis (patient) is inverted when created by the Pij matrices. This is a left handed CS where the Z-axis is inverted with respect to the DICOM patient CS
%Flip the Z axis to get the DICOM patient CS

end
