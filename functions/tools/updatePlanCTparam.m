%% updatePlanCTparam
% Update |Plan.DoseGrid| with the information from |handles|
%
%% Syntax
% |Plan  = updatePlanCTparam(handles , Plan  )|
%
%
%% Description
% |Plan  = updatePlanCTparam(handles , Plan  )| Description
%
%
%% Input arguments
% |handles| -_STRUCTURE_- REggui data handle
%
% |Plan| -_STRUCTURE_- Information about the treament plan
%
%
%% Output arguments
%
% |Plan| -_STRUCTURE_- Information about the treament plan
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function Plan  = updatePlanCTparam(handles, Plan)

  CT = Get_reggui_data(handles,Plan.CTname,'images');
  Plan.DoseGrid.size = size(CT);
  Plan.DoseGrid.resolution = handles.spacing ; % mm
  Plan.DoseGrid.nvoxels = prod(Plan.DoseGrid.size);

  Plan.CTinfo.ImagePositionPatient = handles.origin; %Update origin of the CT scan
  Plan.CTinfo.Spacing = handles.spacing;
end
