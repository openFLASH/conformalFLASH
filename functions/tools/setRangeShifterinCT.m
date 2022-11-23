%% setRangeShifterinCT
% Add the range shifter in the CT scan for all beams with tag "RSinfo".
% The cross section width of the range shifter is large enough to cover the whole area of the aperture block
% If |minField| and |maxField| are provided, only the cross section defined by these parameter is inserted in the CT.
% The CT scan is the image with name |CTname| stored in |handles|
% Expand the CT if it is too small.
% Save the updated CT in the |handles| with name |CTname|
% Update the |Plan| structure with the update CT size
%
%% Syntax
% |[Plan , handles] = setRangeShifterinCT(handles , Plan , CTname , minSize , maxSize)|
%
%
%% Description
% |[Plan , handles] = setRangeShifterinCT(handles , Plan , CTname , minSize , maxSize)| Description
%
%
%% Input arguments
%
% |handles| -_STRUCTURE_- REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |Plan| -_STRUCTURE_- Information about the treament plan
%  * |Plan.ScannerDirectory| - _STRING_ - Name of the folder containing the definition of the CT scanner properties in MCsquare in folder "plugins\openMCsquare\lib\Scanners"
%  * |Plan.DoseGrid| - _struct_ - Structure containing the information about the dose grid. The following data must be present:
%  * |Plan.DoseGrid.nvoxels| - _scalar_ - Total number of voxels of the dose grid, i.e., number of voxels in the CT image.
%  * |Plan.CTinfo| -_STRUCTURE_- DICOM header of the CT scan
%  * |Plan.Beams| -_VECTOR of STRUCTURES_- Information about the different beams in the plan
%     * |Plan.Beams(b).GantryAngle| -_SCALAR_- Gantry Angle (deg)
%     * |Plan.Beams(b).PatientSupportAngle| -_SCALAR_- Couch Angle (deg)
%     * |Plan.Beams(b).isocenter| -_SCALAR VECTOR_- [x,y,z] Coordiantes (mm) of the isocentre in the planning CT scan
%     * |Plan.Beams(b).BlockData| Coordinates (mm) of the i-th point defining the contour of the aperture block projected onto the machine isocentric plane in the IEC BEAM LIMITING DEVICE coordinate system for the b-th beam.
%     * |Plan.Beams(b).RSinfo| -_STRUCTURE_- Information about the range shifter
%        * |Plan.Beams(b).RSinfo.RangeShifterMaterial| - _STRING_ - Name of the range shifter material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%        * |Plan.Beams(b).RSinfo.RSslabThickness| -_SCALAR_- mm Thickness of the individual slabes of the range shifter
%        * |Plan.Beams(b).RSinfo.NbSlabs| -_SCALAR_- Number of slabs of range shifter
%
% |CTname| -_STRING_- Name of the CT image in handles.images into which the device is to be inserted
%
% |minField| -_SCALAR VECTOR_- [OPTIONAL, only needed if part of the CEM is to be used] [X, Y] Coordinate (mm, in IEC gantry) of [-x,-y] the corner of the field
%
% |maxField| -_SCALAR VECTOR_- [OPTIONAL, only needed if part of the CEM is to be used] [X, Y] Coordinate (mm, in IEC gantry) of [+x,+y] the corner of the field
%
%
%% Output arguments
%
% |handles| -_STRUCTURE_- Updated REggui data handle. The CT scan is stored in |handles| in the image with name |Plan.CTname|.
%
% |Plan| -_STRUCTURE_- Updated structure defining a multi energy layer PBS plan
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function [Plan , handles] = setRangeShifterinCT(handles , Plan , CTname , minSize , maxSize)

  if nargin < 4
    minSize = [];
    maxSize = [];
  end

  CT = Get_reggui_data(handles,CTname,'images'); %Update the CT scan with the aperture block in handles
  ImagePositionPatient = handles.origin;
  PixelSize = min(handles.spacing);


  HUair =  getMaterialSPR('Schneider_Air' , Plan.ScannerDirectory) + 1; %Hounsfield unit associated to air in the material file

  for b = 1: size(Plan.Beams,2) %Loop for each beam

    if isfield(Plan.Beams , 'RSinfo')
            %Add the range shifter only if there is one defined
            if isempty(minSize) | isempty(maxSize)
                [minSize , maxSize] = findApertureBlockSize(Plan.Beams(b).BlockData );
            end

            Xrs =  minSize(1):PixelSize./2:maxSize(1);
            Yrs =  minSize(2):PixelSize./2:maxSize(2);
            [Yrs,Xrs] = meshgrid(Yrs,Xrs);

            Ap = [Xrs(:),Yrs(:)]; %X,Y coordinate of voxels containing brass

            HUrangeshifter =  getMaterialSPR(Plan.Beams(b).RSinfo.RangeShifterMaterial , Plan.ScannerDirectory); %HU and relative stopping power of the range shifter
            HUrangeshifter = HUrangeshifter + 1;


            Ars = []; %Indices of the voxels containing brass
            %Divide the step by /2 in order to avoid sampling issue leading to missing planes in the range shifter
            for slab = Plan.Beams(b).RSinfo.NbSlabs : -1 : 1
                fprintf('Slab # %d \n',slab)
                for step = PixelSize:PixelSize./2:Plan.Beams(b).RSinfo.RSslabThickness(slab)
                  Zg = Plan.Beams.RSinfo.IsocenterToRangeShifterDistance + Plan.Beams(b).RSinfo.RSslabThickness(Plan.Beams(b).RSinfo.NbSlabs) + ... %distance isocentre to upstream side of RS slab
                        Plan.Beams(b).RSinfo.SlabOffset(slab)  + ... %Jump to the upstream side of the slab number |slab|
                        - step ; %Add a layer to the slab from upstream to downstream surface
                                      %Z in IEC gantry at which this layer of CEF is located
                                      % Plan.Beams.RSinfo.IsocenterToRangeShifterDistance defines the DOWNSTREAM surface. The other slices of the slab are at small Zg distance. Hence + step
                  Ars = [Ap , ones(size(Ap,1),1).* Zg , ones(size(Ap,1),1)]; %Coordinate of the voxels of the slab of the range shifter
                  [CT , ImagePositionPatient] = insertDeviceInCT(CT , Ars , HUrangeshifter , Plan.Beams(b) , handles.spacing , ImagePositionPatient , HUair);
                end % for step
            end %for slab
      end % if isfield(Plan.Beams , 'RSinfo')
  end %  for beam

  %Update the handles and plan
  [handles , Plan] = updateAllImages(handles , Plan , CT , ImagePositionPatient , HUair  , CTname);

end
