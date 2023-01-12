%% weight2charge
% Converts weight (MU computed by MCsquare) to proton charge
%
%% Syntax
% |plan = weight2charge(plan, bdl, newEnergy)|
%
%
%% Description
% |plan = weight2charge(plan, bdl, newEnergy)| Description
%
%
%% Input arguments
% |plan| -_STRUCTURE_- plan
%
% |bdl| -_STRING_- path to the BDL file
%
% |newEnergy| -_SCALAR_- [optional] energy (MeV) to normalize the charge (for weight w at energy E, charge = charge(w, E)/charge(1, newEnergy)
%
%
%% Output arguments
%
% plan - plan with charges instead of weights
%
%
%% Contributors
% Authors : S. Deffet (open.reggui@gmail.com)


function plan = weight2charge(plan, bdl, newEnergy)

	[pluginPath , MCsqExecPath , BDLpath , MaterialsPath , ScannersPath] = get_MCsquare_folders();
	BDLfullPath = fullfile(BDLpath,bdl);
	T = textread(BDLfullPath,'%s','delimiter','\n');
  T{length(T)+1} = -1;
  mask_header = find(~cellfun(@isempty,strfind(T(1:length(T)-1),'NominalEnergy')));
  BDL_parameters = importdata(BDLfullPath, ' ', mask_header);
  BDL_e = BDL_parameters.data(:, 1);
  BDL_mu = BDL_parameters.data(:, 4);

  if nargin>2
      prot2 = interp1(BDL_e, BDL_mu, newEnergy, 'linear', 'extrap'); %proton/MU at energy |newEnergy|
  else
      prot2 = 1;
  end

  if isfield(plan, 'Beams')
      for b=1 : length(plan.Beams)
          for l=1 : length(plan.Beams(b).Layers)
              en = plan.Beams(b).Layers(l).Energy;
              w = plan.Beams(b).Layers(l).SpotWeights; %Weight per fraction of the spot

              prot = interp1(BDL_e, BDL_mu, en, 'linear', 'extrap'); %proton/MU at energy of the layer
              plan.Beams(b).Layers(l).SpotWeights = w .* prot./ prot2;
          end
      end
  elseif isfield(plan{1}, 'spots')
      for b=1 : length(plan)
          for s=1 : length(plan{b}.spots)
              en = plan{b}.spots(s).energy;
              w = plan{b}.spots(s).weight;

              prot = interp1(BDL_e, BDL_mu, en, 'linear', 'extrap');
              plan{b}.spots(s).weight = w .* prot ./ prot2;
          end
      end
  else
      error('Unknown format for plan');
	end
end
