%% makeRangeShifter4snout
% Compute the number of slabs of range shifter and the residual range shifting in CEM
% for the currently selected snout
%
%% Syntax
% |Plan = makeRangeShifter4snout(Plan)|
%
%
%% Description
% |Plan = makeRangeShifter4snout(Plan)| Description
%
%
%% Input arguments
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are
%   * |Plan.BDL| -_STRING_- Beam data library. Name of the folder in REGGUI\plugins\openMCsquare\lib\BDL
%   * |Plan.Beams(b).RSinfo.R_max| -_SCALAR_- Range (cm) in water to reach the distal surface of the PTV
%   * |Plan.Spike.MaterialID| - _STRING_ - Name of the CEM material, as defined in the file "plugins\openMCsquare\lib\Materials\list.dat"
%   * |Plan.Spike.min_thickness| -_SCALAR_- Thickness (mm) in of the base on which the spikes are built. The base has the same |R_WET| as the spikes
%
%
%% Output arguments
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are
%   * |Plan.Beams(b).CEFbaseWET| -_SCALAR_- WET (mm) of the base of the CET acting as a range shifter
%   * |Plan.Beams(b).CEFbaseThickness| -_SCALAR_- Thickness (mm) of the base of the CEF
%   * |Plan.Beams(b).RSinfo| -_STRUCT_- Information about the range shifter
%       * |Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance| -_SCALAR_- Isocenter to DOWNSTEAM edge of range shifter (mm).
%       * |Plan.Beams(b).RSinfo.NbSlabs| -_SCALAR_- Number of slabs of range shifter
%       * |Plan.Beams(b).RSinfo.RSslabThickness|  -_SCALAR VECTOR_- |RSslabThickness(s)| Thickness (mm) of the s-th slab of the range shifter
%       * |Plan.Beams(b).RSinfo.RangeShifterMaterial| -_STRING_- Material of the range shifter
%       * |Plan.Beams(b).RSinfo.RangeShifterWET| -_SCALAR_- water equivalent thickness (mm) of the range shifter
%       * |Plan.Beams(b).RSinfo.SlabOffset| -_SCALAR VECTOR_-  Offset (mm) between the last slab and i-th slab
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Plan = makeRangeShifter4snout(Plan)

if size(Plan.Beams,2) > 1
  error('Plan cannot have more than 1 beam')
end

  for b = 1: size(Plan.Beams,2) %Loop for each beam
    %Define the range shifter thickness and the CEF base thickness
    [RangeShifterWET , R_beam , CEF_WET , NbSlabs , CEFbaseThickness , RSslabThickness , RangeShifterMaterial] = getRangeShifterWet(Plan, Plan.BDL, Plan.Beams(b).RSinfo.R_max , Plan.Spike.MaterialID , Plan.Spike.MinThickness);

    if (CEFbaseThickness < Plan.Spike.MinThickness)
      %The residual range is smaller than the minimum thickness of the CEF
      CEFbaseThickness
      Plan.Spike.MinThickness
      error('The base of the hedgehog is too thin')
      %TODO The solution is to remove one slab of range shifter and add a big thickness of CEF
    end
    fprintf('Max range machine : %3.2f cm \n', R_beam)
    fprintf('WET to PTV        : %3.2f cm \n', Plan.Beams(b).RSinfo.R_max)
    fprintf('Range shifter WET : %3.2f cm \n', RangeShifterWET)
    fprintf('Number of range shifter slabs : %d \n',NbSlabs)
    fprintf('CEF           WET : %3.2f cm \n', CEF_WET)
    fprintf('Thickness of CEF base : %3.2f mm \n',CEFbaseThickness)
    if NbSlabs > 0
        Plan.Beams(b).RSinfo.NbSlabs = NbSlabs;
        Plan.Beams(b).RSinfo.RSslabThickness = RSslabThickness;
        Plan.Beams(b).RSinfo.RangeShifterMaterial = RangeShifterMaterial; %REad range shifter material from snout information
        Plan.Beams(b).RSinfo.RangeShifterWET = RangeShifterWET .* 10; %water equivalent thickness (mm) of the range shifter

        %Define the position des range shifter slabs from the snout configuration
        param = getMachineParam(Plan.BDL); %GEt information about geometry of FLASH snout
        Plan.Beams(b).RSinfo.SlabOffset = param.snout.RangeShifterOffset(1:Plan.Beams(b).RSinfo.NbSlabs) - param.snout.RangeShifterOffset(1) +  Plan.Beams(b).RSinfo.RSslabThickness(1); %Offset from |IsocenterToRangeShifterDistance| and the upstream side of the i-th slab

        %Compute isocentre to range shifter distance
        Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance = getIsocenterToRangeShifterDistance(Plan.Beams(b) );
        fprintf('Isocenter To downstream side of Range Shifter : %3.2f mm\n',Plan.Beams(b).RSinfo.IsocenterToRangeShifterDistance)
      else
        fprintf('There is no range shifter slab \n')
        if isfield(Plan.Beams(b) , 'RSinfo')
          Beam = rmfield(Plan.Beams(b), 'RSinfo');
          Plan.Beams = Beam; %REmove the RSinfo struct from Beams
                %TODO This works only if there is one single beam
          Plan.Beams.NumberOfRangeShifters = 0;
        end
      end

      Plan.Beams(b).CEFbaseWET = CEF_WET .* 10; %WET (mm) of the base of the CET acting as a range shifter
      Plan.Beams(b).CEFbaseThickness = CEFbaseThickness; %Thickness (mm) of the base of the CEF

  end


end
