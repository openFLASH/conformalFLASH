%% getBEVsize
% Compute the X,Y dimension of the beam's eye view in the IEC gantry CS
% based on the coordiantes of the voxels contained in the target (expressed in DICOM CS)
%
%% Syntax
% |res = help_header(im1,im2)|
%
%
%% Description
% |res = help_header(im1,im2)| Description
%
%
%% Input arguments
% |x| -_SCALAR VECTOR_- |x(i)| X coordinate (mm) in DICOM CS of the i-th point in the target
%
% |y| -_SCALAR VECTOR_- |y(i)| Y coordinate (mm) in DICOM CS of the i-th point in the target
%
% |z| -_SCALAR VECTOR_- |z(i)| Z coordinate (mm) in DICOM CS of the i-th point in the target
%
% |Beam| -_STRUCTURES_- Information about the beam
%     * |Beam.gantry_angle| -_SCALAR_- Gantry Angle (deg)
%     * |Beam.table_angle| -_SCALAR_- Couch Angle (deg)
%     * |Beam.isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
%
%
%% Output arguments
%
% |bev_size| - _SCALAR VECTOR_ - [x,y] dimension (mm) of the beam's eye view
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function bev_size = getBEVsize(X,Y,Z , Beam)

  X = X - Beam.isocenter(1);
  Y = Y - Beam.isocenter(2);
  Z = Z - Beam.isocenter(3);
  P = [X(:) , Y(:) , Z(:) , ones(numel(X),1)]'; %The 8 apexes of the paralellipipedic CT scan. In DICOM coordinates

  M = matDICOM2IECgantry(Beam.gantry_angle , Beam.table_angle);
  Pbev = M * P; %Coordinate of the 8 apexes in IC gantry CS

  P2 = max(Pbev(1:2,:),[],2); %Two points at the opposite of the diagonal of the BEV
  P1 = min(Pbev(1:2,:),[],2);
  bev_size = P2 - P1;

end
