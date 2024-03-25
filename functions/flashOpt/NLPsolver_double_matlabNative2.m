%% NLPsolver_hybrid_double_matlabNative
% Perform a robust optimisation of the weight so as to meet the objective function for all scenario.
% The optimisation is:
% minimize [t,w]  t
% subject to  t ≥ f(D(w, s)) ∀ s ∈ S [CONSTRAINTS]
%             w ≥ 0
%
% where S is the set of error scenarios
% and w are the spot weights
% t is an auxillary variable
%
% Note that the constraints are in fact what is usually defined as the objective function in case of a non robust optimisation
%
% The optimisation follows the Minimax optimisation of Fredriksson et al [5]
% The minimax optimization minimizes the objective function in the worst case scenario.
% By considering only physically realizable scenarios, correlations between voxels are preserved
% and the extra conservativeness that re- sults from independent handling of voxels is avoided.
%
%
%% Syntax
% |[x, info] = NLPsolver_hybrid_double_matlabNative(ROI,Plan,Outputs,OptConfig,x0,s_nominal)|
%

%% Description
% [x, info] = NLPsolver_hybrid_double_matlabNative(ROI,Plan,Outputs,OptConfig,x0,s_nominal)| Description
%
%
%% Input arguments
% |ROI| - _struct_ - MIROpt structure containing information about all
% volumes in the RTSTRUCT file. The following data must be present in the
% structure:
%
% * |ROI(i).name| - _string_ - Name of the i-th ROI in the RTSTRUCT list. This must be encoded for all structures, i. e., i = 1, ..., N; where N is the total number of ROIs.
% * |ROI(i).mask1D| - _array_ - Logical column vector storing a binary mask for ROI i (voxels inside the volume of interest are equal to 1, and those outside are equal to zero).
% * |ROI(i).mask3D.value|- _array_ - |ROI(i).mask3D.value(x,y,z)=1| if the voxel at (x,y,z) is located inside the RT struct
% * |ROI(i).nvoxels| - _scalar_ - Total number of voxels in the ROI i.
% * |ROI(i).voxelsID| - _array_ - Column vector containing the linear indices of the non-zero voxels in |ROI(i).mask1D|. At this stage this filed is empty.
%
% |Plan| - _struct_ - MIROpt structure where all the plan parameters are
% stored. The following data must be present in the structure:
%
% * |Plan.DoseGrid| - _struct_ - Structure containing the information about the dose grid. The following data must be present:
% * |Plan.DoseGrid.nvoxels| - _scalar_ - Total number of voxels of the dose grid, i.e., number of voxels in the CT image.
%
% * |Plan.optFunction| - _struct_ - Structure containing the information about the objective functions set by the user on the dose to the target volume and organs at risk. The following data must be present for each objective function with index i:
%     * |Plan.optFunction(i).ROIname| - _string_ - Name of the ROI to which the objective function must be applied.
%     * |Plan.optFunction(i).name| - _string_ - Name of the objective function type. The possible types supported so far are: 'min', 'max', 'min_mean' and 'max_mean'.
%     * |Plan.optFunction(i).Dref| - _scalar_ -  Reference dose in Gy
%     * |Plan.optFunction(i).Vref| - _scalar_ -  Reference volume in 1% only needed for minDVH and maxDVH functions
%     * |Plan.optFunction(i).impw| - _scalar_ -  Importance weight for the current objective function (impw > 0)
%     * |Plan.optFunction(i).robust| - _scalar_ - Boolean value indicating if the objective function must be considered in all scenarios (if set as robust, equal to 1), or only in the nominal case (if set as non-robust, equal to zero).
%
% * |Plan.NbrRandomScenarios| - _scalar_ - Number of scenarios used to simulate the random setup errors.
% * |Plan.NbrRangeScenarios| - _scalar_ - Number of scenarios used to simulate the range errors.
% * |Plan.Nbr4DScenarios| - _scalar_ - Number of scenarios used to simulate the errors related to breathing motion.
% * |Plan.NbrSystSetUpScenarios| - _scalar_ - Number of scenarios used to simulate the systematic setup errors.
% * |Plan.PlanExistentFile| - _string_ - Path pointing to the location of the file containing the spot doses. If empty, the spot doses (dose influence matrix) will be computed from scratch.
%
% |OptConfig| - _struct_ - Structure containing information needed to run
% the optimization. The following data must be present in the structure:
%
% * |OptConfig.max_iter| - _scalar_ - Maximum number of iterations. If this number is reached, the optimization will stop.
% * |OptConfig.plotTargetDVH| - _scalar_ - Boolean value indicating if the DVH for the target must be plotted during the optimization (if equal to 1) or not (if equal to zero).
% * |OptConfig.BeamletMatrixPrecision| - _string_ - Parameter indicating the precision that must be used for the dose influence matrix. It has two values: 'd' if double precision and 'f' if float precision. The second one is recommended to reduce memory usage.
%
% * |OptConfig.mixedOptimization| - _struct_ - Structure containing information related to the mixed optimization. The following data must be present in the structure:
% * |mixedOptimization.ON| - _scalar_ - Boolean value indicating if the mixed optimization is on (if equal to 1) or off (if equal to zero).
% * |mixedOptimization.MC_corrections| - _scalar_ - Maximum number of corrections to be applied in the mixed optimization.
% * |mixedOptimization.iter| - _scalar_ - Maximum number of iterations to perform in between corrections, i.e., in the inner loop of the mixed optimization.
% * |mixedOptimization.RandSetUpError| - _array_ - Random setup errors are sampled for a Gaussian probability distribution. Therefore, this is a row (1x3) vector containing the standard deviation in x, y, and z DICOM axis to be applied in order to correct for random errors in the mixed optimization.
% * |mixedOptimization.Opt4D| - _struct_ - Structure containing information regarding the simulation of the breathing motion in order to correct for it using the mixed optimization. This feature is not supported so far.
%
% * |s_nominal| - _scalar vector_ - Index for the global nominal case. There will be a s_nominal(i4D) for each 4D scenario, which corresponds to shift=0, range=0 and sigma(random)=0 in case there are several sigmas or either the given sigma(random) if there is only one.

%% Output arguments
% |x| - _SCALAR VECTOR_ - |x(i)| optimised wieght of the i-th spot.
%                           D(i) = P*x(2:end) gives the TOTAL dose at the i-th voxel.
%                           The optimizer tries and make x(1)=worstCase  as small as possible (because funcs.objective returns x(1)). X(1) defines the worst acceptable sum of deviation of each objective for each one scenario
%
% |info| - _STRUCTURE_ - Information returned by the |ipopt| optimizer
%
% |Plan| -_STRUCTURE_- Updated PLan structure with pre-computed values
%   *|Plan.SpotTrajectoryInfo| -_STRUCTURE_- [OPTIONAL: if provided, the spot trajectory is not optimised] Pre-computed spot trajectory. Used to accelerate computations.
%                         |SpotTrajectoryInfo{b}| Parameter for beam |b|
%     * |Plan.SpotTrajectoryInfo.beam{b}.sobpSequence| -_SCALAR VECTOR_-  Order of the indices of |spot| to sort the spots. |OrderedSOBP = spot(sobpSequence,:)|
%     * |Plan.SpotTrajectoryInfo.beam{b}.Nmaps| -_STRUCTURE_-  Topological information on the initial spot sequence
%       * |Nmaps.NeighbourghMap| -_SCALAR VECTOR_- |NeighbourghMap(d,i)| d=# of delivered spot; i=# of impacted spot; |NeighbourghMap(d,i)|  = 0 if there is an impact
%       * |Nmaps.NeighbourghWeightMap| -_SCALAR VECTOR_- |NeighbourghWeightMap(d,i)| fraction of the dose of spot d (# of delivered spot) that is also delivered at spot i (i=# of impacted spot)
%       * |Nmaps.NeighbourghTimeMap| -_SCALAR MATRIX_- |NeighbourghTimeMap(d,i)| Time (ms) required to go from spot d to spot i
%     * |Plan.SpotTrajectoryInfo.sobpPosition{b}| - CELL VECTOR_ - beamletPosition{b}(i,:) = [x,y] The i-th spot of the b-th beam is located at position [x,y] in the BEV coordinates
%     * |Plan.SpotTrajectoryInfo.weight2spot| - SCALAR MATRIX_ - |weight2spot(weightIdx,:) = [b,BeamLetNb,l]| The spot at index |weightIdx| is related to the b-th beam and the SOBP # |BeamLetNb| in layer # l
%
%
%% REFERENCES
% [1] https://github.com/coin-or/Ipopt
% [2] https://coin-or.github.io/Ipopt/
% [3] https://projects.coin-or.org/Ipopt/wiki/MatlabInterface
% [4] https://dial.uclouvain.be/pr/boreal/object/boreal:189076
% [5] https://aapm.onlinelibrary.wiley.com/doi/epdf/10.1118/1.3556559
% [6] https://github.com/e0404/matRad/issues/408
% [7] https://openreggui.org/git/open/REGGUI/issues/20
%
%% Contributors
% Author(s): Rudi Labarbe (for the FLASH part of the function)

function [x, info , Plan] = NLPsolver_double_matlabNative2(ROI,Plan,OptConfig, warm_start_in,s_nominal,handles)

  loaded = false;
  if isfield(Plan,'PlanExistentFile')
    if (exist(fullfile(Plan.PlanExistentFile,'x.mat'),'file') && exist(fullfile(Plan.PlanExistentFile,'info.mat'),'file'))
      fprintf('Loading previously optimised weights ...')
      data = load(fullfile(Plan.PlanExistentFile,'x'));
      x = data.x;
      data = load(fullfile(Plan.PlanExistentFile,'info'));
      info = data.info;
      fprintf('done \n')
      loaded = true;
    end
  end

  if exist(fullfile(Plan.PlanExistentFile,'spotSequence.mat'),'file')
    fprintf('Loading previously optimised trajectory ...')
    data = load(fullfile(Plan.PlanExistentFile,'spotSequence'));
    Plan.SpotTrajectoryInfo = data.GBL_SpotTrajectoryInfo;
    fprintf('done \n')
  end

  if (~loaded)
    %No previous computation available. Recompute optimisation

    % Total number of scenarios
    Plan.nscenarios = Plan.NbrRangeScenarios * Plan.NbrSystSetUpScenarios * Plan.NbrRandomScenarios * Plan.Nbr4DScenarios;

    %..................
    % global variables
    %..................
    GLB_store_fobj = 0;
    GLB_backUp.scenarioDose = struct('D',sparse(zeros(Plan.DoseGrid.nvoxels,1)));
    GLB_backUp.w = ones(1,Plan.nominalSpots)*(-1);% to force to calculate the dose at the first iteration
    GLB_wBackUp_counter = 1;
    GLB_wBackUp_Niter = 25; % save w values and plot DVH for target volume every Niter iterations

    eV = 1.6021766208e-19 ; % J/eV
    param = getMachineParam(Plan.BDL);
    if isfield(Plan , 'Inozzle')
      GLB_alpha = 1000 .* MU_to_NumProtons(1, param.MAXenergy) .* eV ./ (Plan.Inozzle .* 1e-9); % ms / MU :   aplha * w = T
    end

    %Get the (x,y) position of each spot in the weight vector
    GLB_DRa = []; % Average percentile dose rate
    GLB_DRm = []; % Median  percentile dose rate
    GLB_DADRm = []; % Median dose averaged dose rate
    GLB_DRADm = [];
    GLB_SpotTiming=[];  %|SpotTiming{b}(i)| Time (ms) taken at the i-th pixel to deliver the |percentile| of the dose for the b-th beam

    %Loading pre-computed projection if they exists
    if isfield(Plan,'PlanExistentFile') && (exist(fullfile(Plan.PlanExistentFile,'planRidge.mat'),'file'))
          fprintf('Loading projection data \n')
          data = load(fullfile(Plan.PlanExistentFile,'planRidge.mat'));
          Plan.bevPTV = data.Plan.bevPTV;
          Plan.bev_x = data.Plan.bev_x;
          Plan.bevOAR = data.Plan.bevOAR;
    end

    %If scanAlgo is defined, then get an estimate of beam current from scanAlgo
    if isfield(Plan, 'scanAlgoGW')
      fprintf('----- Getting Inozzle from scanAlgo gateway: %s \n',Plan.scanAlgoGW.scanalgoGateway_IP)
      Plan.Inozzle  = getCurrentFromScanAlgo(Plan.scanAlgoGW , Plan.BDL);
      fprintf('Nozzle current from scanAlgo : %f nA \n',Plan.Inozzle);
    end

    %Check whether the objective function requires the dose rate
    %If so, then we need to compute trajectory information
    flagUseDoseRate = 0; %Flag is true if anyone function requires dose rate computation
    for optFidx = 1:length(Plan.optFunction)
     flagUseDoseRate = flagUseDoseRate + UseDoseRate(Plan.optFunction(optFidx).ID);
   end
   if flagUseDoseRate
     %The trajectory is optimized once.
     %Only the spot weights will be optimised afterwards
     GBL_SpotTrajectoryInfo = optimizeTrajectory(Plan  , ROI);
     Plan.SpotTrajectoryInfo = GBL_SpotTrajectoryInfo;
   end

    if (OptConfig.plotTargetDVH == 1)
        h = figure;
        title('DVH for target volume');
        xlabel('Dose (Gy)');
        ylabel('% Volume');
        axis([0,100,0,100]);
    end

    % test
    for s = 1:Plan.NbrSystSetUpScenarios
        Plan.SystSetUpScenario(s).ws = zeros(Plan.virtualSpots,1);
        Plan.SystSetUpScenario(s).ws(Plan.SystSetUpScenario(s).wID) = 1;
    end

    if warm_start_in(1).status == 1

        disp(' Using warm start info to continue previous optimization...');

        % Initialize the dual point.
        options.zl     = warm_start_in.zl';
        options.zu     = warm_start_in.zu';
        options.lambda = warm_start_in.lambda;
        % Initialize the primal variables
        x0 = warm_start_in.varstruct.x;


        % Set IPOPT options

        % WARM START
        options.ipopt.warm_start_init_point = 'yes';
        options.ipopt.warm_start_bound_push = 1e-6; % default = 0.001
        options.ipopt.warm_start_mult_bound_push = 1e-6; % default = 0.001

        % Barrier parameter
        options.ipopt.mu_strategy = 'monotone';% default = 'monotone'
        options.ipopt.mu_init = warm_start_in.varstruct.mu; % default 0.1, recommended for warm start: 1e-6

    else

        disp(' Starting optimization from scratch...');

        % The starting point. (ROW vector)
        w0 = ones(1,Plan.nominalSpots)*0.1; %Initial value of the weights

        % INITIALIZE t = objective function for the worst case scenario
        checkBackUp([0 w0]); %Compute the dose for the weight w0 and store results in backUp.scenarioDose (a global variable)
        fobj = zeros(Plan.nscenarios,1);
        s = 1;
        for i_4D = 1:Plan.Nbr4DScenarios
            for rr = 1:Plan.NbrRandomScenarios
                for ss = 1:Plan.NbrSystSetUpScenarios % syst. setup scenarios = ss
                    for rs = 1:Plan.NbrRangeScenarios % range scenario = rs

                        for i = 1:length(Plan.optFunction) %Loop for every objective
                            if (Plan.optFunction(i).robust == 1)
                                fobj(s) = fobj(s) + evalFunction(Plan.optFunction(i),GLB_backUp.scenarioDose(s).D, i_4D , rr , rs , i);
                            else
                                if (s == s_nominal(i_4D))
                                    fobj(s) = fobj(s) + evalFunction(Plan.optFunction(i),GLB_backUp.scenarioDose(s).D, i_4D , rr , rs , i);
                                else
                                    fobj(s) = fobj(s) + 0;
                                end %if (s
                            end %if (Plan.optFunction(i)
                        end %for i = 1
                        s = s + 1;
                    end

                end
            end
        end


        % take worst case value
        t0 = max(fobj); %The objective is defined as the worst unmet constraint

        %initialize variables
        %t0 is worst case deviation to the objective that is allowed (for all objective)
        %w0 is a vector defining the weight of each spot
        x0 = [t0 w0];

        % Options

        % Initialization
        options.ipopt.bound_frac = 0.01; %default 0.01
        options.ipopt.bound_push = 0.001; %default 0.01

    end

    %% Set the IPOPT options.

    % INFO (from G.Janssens)
    % En ce qui concerne le système IBA, il y avait une limite de 12 MU/spot,
    % mais depuis la dernière version du scanning controller software cette
    % limite a disparu. Si la dose par spot est plus grande que 12 MU, le spot
    % est répété. Parfois, les cliniciens peuvent eux-mêmes limiter la dose par
    % spot dans le TPS. La limite inférieure est de 0.01 MU/spot pour la
    % dedicated nozzle, et 0.02 MU/spot pour l’universal nozzle.

    % minMU = 0.02; % minimum weight value for Universal Nozzle 0.02 MU
    minMU = 0;
    options.lb  = [0 minMU*ones(1,Plan.nominalSpots)];         % Lower bounds on x

    % maxMU = 12; %
    maxMU = Inf;
    options.ub  = [Inf maxMU*ones(1, Plan.nominalSpots)];      % Upper bounds on x

    % The bound on constraint functions.
    options.cl  = zeros(1, Plan.nscenarios);           % Lower bounds on constraints
    options.cu  = ones(1, Plan.nscenarios)*(Inf);      % Upper bounds on constraints
    %The objective is that the value is POSITIVE

    % Output
    options.ipopt.print_level = 5; %maximum print level=12
    options.ipopt.output_file = fullfile(Plan.output_path,'ipopt_output_file.txt');
    options.ipopt.file_print_level = 5;
    options.ipopt.print_info_string = 'yes';
    options.ipopt.print_timing_statistics = 'yes';

    % Termination
    options.ipopt.tol              = 1e-8; % default 1e-8
    options.ipopt.max_iter         = OptConfig.max_iter; % default 3000
    options.ipopt.dual_inf_tol     = 1; % default 1
    options.ipopt.constr_viol_tol  = 0.0001; % default 0.0001
    options.ipopt.compl_inf_tol    = 0.0001; % default 0.0001


    % Quasi-Newton
    options.ipopt.hessian_approximation = 'limited-memory'; % default 'exact'
    options.ipopt.limited_memory_update_type = 'bfgs';
    % By using a big number we are forcing to use complete-memory bfgs since we
    % take into accout all the previous iterations to compute the hessian
    % update (BLGorissen advice)
    options.ipopt.limited_memory_max_history =1000;
    %options.ipopt.limited_memory_max_history =20;

    % DERIVATIVE CHECK
    %options.ipopt.derivative_test = 'first-order'; % Derivatives are checked so far and there are ok, values for relative delta. range from 1e-2 to 1e-4



    %% Print IPOPT options

    save(fullfile(Plan.PlanExistentFile,'NLPsolver_options.mat'),'options', '-v7.3');


    %% Callback functions.

    funcs.objective         = @objective; %Calculates the objective function at the current point.
    funcs.constraints       = @constraints; %t evaluates the constraint functions at the current point for all scenarios. It takes one input, x= [worst case , weights]. The return value is a vector of length equal to the number of constraints (= one per scenario)
    funcs.gradient          = @gradient; %Computes the gradient of the objective at the current point.
    funcs.jacobian          = @jacobian; % Evaluates the Jacobian of the constraints at the current point.
    funcs.jacobianstructure = @jacobianstructure; % It takes no inputs. The return value is a sparse matrix whereby an entry is nonzero if and only if the Jacobian of the constraints is nonzero at ANY point.
    %funcs.iterfunc          = @getOptVars; %callback routine that is called once per algorithm iteration.
                            %TODO [6][7] The new version of IPOPT has a new API which does not support varstruct.x anymore. iterfunc is optiopnal. Let's desactivate it.
                            % This will also remove the warm start.
    %
    % The optimizer  tries and make x(1)=worstCase  as small as possible (because funcs.objective returns x(1)). X(1) defines the worst acceptable sum of deviation of each objective for each one scenario
    % Then funcs.constraints makes sure that the sum of deviations of all objectives are below x(1).  The optimizer makes sure that the following constraints are met : 0 <= x(1) - SUM(RoptFunction(i)) <= +Inf
    % The constraint is met for all scenarios
    % The optimizer is changing the values of X = [worstCase , weights]. It tries, at the same time, to lower the value of the worst case (x(1)) AND to change the value of the weights  (X(2:end))so as to meet the constraints
    %
    % At each iteration, the computation of the new dose map is done by a call to checkBackUp(x) made by funcs.constraints. The new dose map is stored in a globa lvariable that is the used by funcs.constraints

    %% Define warm_start
    warm_start_out = struct('status',[], 'lambda',[], 'zl',[],'zu',[], 'varstruct',{});

    %% Run IPOPT.
    tic;
    [x, info] = ipopt(x0,funcs,options);
    t=toc;

    disp(['Optimization time = ',num2str(t)]);


    % save warm_start_out
    warm_start_out(1).zl = info.zl;
    warm_start_out(1).zu = info.zu;
    warm_start_out(1).lambda = info.lambda;

    save(fullfile(Plan.PlanExistentFile,'warm_start_out'),'warm_start_out');
    save(fullfile(Plan.PlanExistentFile,'x'),'x')
    save(fullfile(Plan.PlanExistentFile,'info'),'info')
    if flagUseDoseRate
      save(fullfile(Plan.PlanExistentFile,'spotSequence'),'GBL_SpotTrajectoryInfo')
    end

end %if (~loaded)

%==============================
% Below are several functions defined inside the body of the main function
% This allows the use of global variables
%==============================

% -----------------------------------------------------------------------
% Compute the dose in whole volume for the wieght(x(2:end))
% The resulting dose is stored in global variable backUp.scenarioDose(s).D
% The weight for each scenario are sotred in global variable Plan.SystSetUpScenario(ss).ws(Plan.SystSetUpScenario(ss).wID)
% The value of global variable GLB_backUp.w  is updated with current value of weights
% -----------------------------------------------------------------------
function checkBackUp(x)

w = x(2:end); %During the optimization, the spot weight is the weight for ALL fractions.
      %It is at the end of the optimization that the spot weight are recomputed per fractions
      % In SpotWeightsOptimization at line 85, the weight are divided by the number of fractions.

% If weights vector has changed, compute the new dose
if any(w ~= GLB_backUp.w)
    s = 1;
    for i_4D = 1:Plan.Nbr4DScenarios
        for rr = 1:Plan.NbrRandomScenarios
            for ss = 1:Plan.NbrSystSetUpScenarios % syst. setup scenarios = ss
                Plan.SystSetUpScenario(ss).ws(Plan.SystSetUpScenario(ss).wID) = w;
                for rs = 1:Plan.NbrRangeScenarios % range scenario = rs
                    % 1D vector for the total dose in CT
                    GLB_backUp.scenarioDose(s).D = Plan.Scenario4D(i_4D).RandomScenario(rr).RangeScenario(rs).P * Plan.SystSetUpScenario(ss).ws; %Compute the dose map for this scenario
                                  %This is the total dose for ALL fractions
                    for optFidx = 1:length(Plan.optFunction)
                      if UseDoseRate(Plan.optFunction(optFidx).ID)
                        if ~isfield(Plan.optFunction(optFidx),'Vref')
                          %Use the default percentile
                          percentile = 0.01; %Default percentile for dose rate
                        else
                          percentile = (1 - Plan.optFunction(optFidx).Vref) ./2;
                        end

                          %The trajectory is already optimised. No need to do it again
                          [ GLB_DRa{i_4D,rr,rs,optFidx} , ~ , GLB_DRm{i_4D,rr,rs,optFidx} , ~ , GLB_DADRm{i_4D,rr,rs,optFidx} , ~,  GLB_DRADm{i_4D,rr,rs,optFidx}, GLB_SpotTiming{i_4D,rr,rs,optFidx}]  = getDRa(GBL_SpotTrajectoryInfo , Plan.SystSetUpScenario(ss).ws(Plan.SystSetUpScenario(ss).wID) , Plan.Scenario4D(i_4D).RandomScenario(rr).RangeScenario(rs).P , Plan , Plan.optFunction(optFidx).Dref  , ROI(Plan.optFunction(optFidx).ROIindex).mask1D , sparse(GLB_backUp.scenarioDose(s).D), Plan.optFunction(optFidx).DMF, Plan.optFunction(optFidx).DR50, percentile , Plan.optFunction(optFidx).ROIname, []);
                      end
                    end
                    s = s + 1;
                end
            end
        end
    end
    GLB_backUp.w = w;
end

end

% ------------------------------------------------------------------
% OBJECTIVE FUNCTION
% ------------------------------------------------------------------
function f = objective (x)
f = x(1);
%The objective is the first element of the x vector
% The optimizer will try to make x(1) as small as possible. This is defining the worst deviation to the constraints
% Indeed, the constraints are computed as difference between x(1) (=wost case) and the value of the constraint function X(1) - C(w)
% The optimiser makes sure that 0 < X(1) - C(w) < +Inf i.e. the difference is POSITIVE, i.e. all the constaints are below the maximum deviation allowed by x(1)
end


% ------------------------------------------------------------------
% CONSTRAINTS function
% ------------------------------------------------------------------
% For each scenario, this functon computes the SUM of the deviation of each contraint to its target X(1)
% Then, for each scenario, it computes the difference between the maximum allowed (sum of ) deviation and the value of this sum for the scenario
% The function return a vector c(s): the difference (Cmax - sum(constraint(s)) to the max deviation to the constraint allowed for the scenario s
% The constraints are defined so that all c(s) are positive. Therefore, the sum of all deivations to the objective will be smaller that the imposed worst case % for all sceario
%
% INPUT  : x = [worst , weights]. X(1) is the budget for errors
% OUTPUT : c(scenarios)
% ------------------------------------------------------------------
function c = constraints (x)

checkBackUp(x); %When the constraint function is called, the dose in whole volume is recomputed by |checkBackUp| and the results are saved in the global variable |backUp.scenarioDose(s).D|

c = ones(Plan.nscenarios,1)*x(1); %Worst case objective
s = 1; %Running index on the different scenarios
for i_4D = 1:Plan.Nbr4DScenarios
    for rr = 1:Plan.NbrRandomScenarios
        for ss = 1:Plan.NbrSystSetUpScenarios % syst. setup scenarios = ss
            for rs = 1:Plan.NbrRangeScenarios % range scenario = rs

                for i = 1:length(Plan.optFunction) %Loop for every ojective function
                    if (Plan.optFunction(i).robust == 1)
                        c(s) = c(s) - evalFunction(Plan.optFunction(i),GLB_backUp.scenarioDose(s).D, i_4D , rr , rs , i);
                    else
                        if (s == s_nominal(i_4D))
                          %Compute the value of the i-th objective for dose map |GLB_backUp.scenarioDose(s).D|
                          % Subtract the value of the deviation of the i-th objective from the budget of errors
                          %The optimizer makes sure that c(s) remains positive for all scearios s
                            c(s) = c(s) - evalFunction(Plan.optFunction(i),GLB_backUp.scenarioDose(s).D , i_4D , rr , rs , i);
                        %else
                        %    c(s) = c(s) - 0;
                      end %if (s ==
                    end %if (Plan.optFunction(i)
                end %for i
                s = s + 1;
            end

        end
    end
end
end


% ------------------------------------------------------------------
%Gradient of the objective function at point x
% ------------------------------------------------------------------
function g = gradient (x)
  %The objective function depends only on the first component of the vector
g = [1 zeros(1,length(x) - 1)];
end

% ------------------------------------------------------------------
% The return value is a sparse matrix whereby an entry is nonzero
% if and only if the Jacobian of the constraints
% is nonzero at ANY point.
% ------------------------------------------------------------------
function J = jacobianstructure ()
J = sparse(ones(Plan.nscenarios,Plan.nominalSpots + 1));
end

% ------------------------------------------------------------------
% Jacobian of the constraint at point weight = x
%It takes one input, x. The output must always be an M x N
%      sparse matrix, where M is the number of scenario and N is the
%      number of variables. Type HELP SPARSE for more information on
%      constructing sparse matrices in MATLAB.
%
% J(Constraints_{s} , weight_i) = dConstraints_{s}(x) / dweight_i
% There is one constraint per scenario s.
% ------------------------------------------------------------------
function J = jacobian (x)
checkBackUp(x);

% J(Constraints_{s} , weight_i) = d Constraints_{s}(x) / dweight_i
J = zeros(Plan.nscenarios,Plan.nominalSpots + 1);

% One row for each scenario, one column for each variable.
s = 1;
for i_4D = 1:Plan.Nbr4DScenarios
    for rr = 1:Plan.NbrRandomScenarios
        for ss = 1:Plan.NbrSystSetUpScenarios % syst. setup scenarios = ss
            for rs = 1:Plan.NbrRangeScenarios % range scenario = rs

                J(s,1) = 1; % derivative with respect to the auxiliar variable t
                % initialize fder
                fder = 0;
                for i = 1:length(Plan.optFunction)
                    if (Plan.optFunction(i).ID ~= 8 && Plan.optFunction(i).ID ~= 9 && Plan.optFunction(i).ID ~= 10)

                        tmp = zeros(Plan.DoseGrid.nvoxels,1);

                        % Jacobian for a dose constraint.
                        if (Plan.optFunction(i).robust == 1)
                            tmp(ROI(Plan.optFunction(i).ROIindex).mask1D) = evalDerivative(Plan.optFunction(i),GLB_backUp.scenarioDose(s).D,i , Plan.Scenario4D(i_4D).RandomScenario(rr).RangeScenario(rs).P(:,Plan.SystSetUpScenario(ss).wID));
                            %tmp(pixel) = evalDerivative = d Constraints_{s}(x) / dpixel
                        else
                            if (s == s_nominal(i_4D))
                                tmp(ROI(Plan.optFunction(i).ROIindex).mask1D) = evalDerivative(Plan.optFunction(i),GLB_backUp.scenarioDose(s_nominal(i_4D)).D,i , Plan.Scenario4D(i_4D).RandomScenario(rr).RangeScenario(rs).P(:,Plan.SystSetUpScenario(ss).wID));
                            end
                        end

                        % evalDerivative computes d Constraints_{s}(x) / dpixel)
                        % Multiply  (d Constraints_{s}(x) / dpixel) * P(pixel , weight))
                        tmp = tmp' * Plan.Scenario4D(i_4D).RandomScenario(rr).RangeScenario(rs).P(:,Plan.SystSetUpScenario(ss).wID);
                        %tmp(w) = evalDerivative = d Constraints_{s}(x) / d_w

                    else
                        % Jacobian for a dose rate constraint
                        if (Plan.optFunction(i).robust == 1)
                            tmp = evalDerivative(Plan.optFunction(i),GLB_backUp.scenarioDose(s).D,i , Plan.Scenario4D(i_4D).RandomScenario(rr).RangeScenario(rs).P(:,Plan.SystSetUpScenario(ss).wID) , GLB_SpotTiming{i_4D,rr,rs,i});
                            %tmp(w) = evalDerivative = d Constraints_{s}(x) / d_w
                        else
                            if (s == s_nominal(i_4D))
                                tmp = evalDerivative(Plan.optFunction(i),GLB_backUp.scenarioDose(s_nominal(i_4D)).D,i , Plan.Scenario4D(i_4D).RandomScenario(rr).RangeScenario(rs).P(:,Plan.SystSetUpScenario(ss).wID) , GLB_SpotTiming{i_4D,rr,rs,i});
                            end
                        end
                        % evalDerivative computes d Constraints_{s}(x) / d_weight_i
                        % tmp is OK. No need for further processing
                    end
                    fder = fder - tmp; %derivative of the sum is the sum of the derivatives
                    %fder = sum_opt_func (d Constraints_{s}(x) / dpixel)
                end
                %d Constraints_{s}(x) / d_weight_i = sum_pixel(sum_opt_func (d Constraints_{s}(x) / dpixel) * P(pixel , weight))
                %  * sum_pixel is implemented by "tmp' * Pij" which uses the MAtlab * matrix multiplication with 2 vectors  (fder(pxl) is a 1D vector)
                %  * sum_opt_func is implemented by the "for i = 1:length(Plan.optFunction)" loop and the accumulation in the fder variable by the line "fder = fder - tmp"
                % Indeed: pixel = P(pixel , weight) .* weight => d_pixel / d_weight = P(pixel , weight)
                %J(s,2:end) = fder' * Plan.Scenario4D(i_4D).RandomScenario(rr).RangeScenario(rs).P(:,Plan.SystSetUpScenario(ss).wID);
                J(s,2:end) = fder;
                s = s + 1;
            end
        end
    end
end

J =  sparse(J);
end


% -----------------------------------------------------------------
% Call back function called at periodic time during the iterations
% -----------------------------------------------------------------
function next = getOptVars (t, f, varstruct)
% INPUT:
% t = current iteration number
% f = current objective function value
% varstruct = struct containing the fields: x, inf_pr,inf_du,mu,d_norm,regularization_size,alpha_du,alpha_pr and ls_trials.

% if t == 0
%     varstruct = warm_start_in.varstruct;
% end % This changes the variables only here inside the function, but the information is lost when we leave it..which is not useful
% ASK IN IPOPT FORUM HOW TO CHANGE THE VALUE OF THESE VARIABLES FOR
% WARM_START.
warm_start_out(1).varstruct = varstruct;

% Terminate the optimzation if requirements* are fullfilled
% if (abs(old_f - f) <= 5e-10)
%         continue = false;
% else
%     old_f = f;
    next = true;
% end

% * the difference between the objective function of two consecutive
% iterations is smaller than a given value.

if (t == GLB_wBackUp_counter*GLB_wBackUp_Niter)

    GLB_wBackUp_counter = GLB_wBackUp_counter + 1;

    %save w to file
    w = varstruct.x(2:end);
    save(fullfile(Plan.PlanExistentFile,'w'),'w');

    % Plot target DVH in nominal case
    if (OptConfig.plotTargetDVH == 1)
        d = GLB_backUp.scenarioDose(s_nominal).D(logical(ROI(Plan.TargetROI_ID).mask1D));
        int_bins = min(d):0.005:max(d); %bin centers
        int_hist = hist(d,int_bins); % in the future we have to take into account the relative volume of each voxel (SEE plotQvhGeneral in DPBNmodule)
        DVH = 100/sum(int_hist)*flip(cumsum(flip((int_hist'),1)),1);
        h = plot(int_bins',DVH,  'g-');
        drawnow;
    end

end

if (OptConfig.plotObjective == 1)

    o = 0;
    for i = 1:length(Plan.optFunction)
        o = o + evalFunction(Plan.optFunction(i),GLB_backUp.scenarioDose(s_nominal).D); %TODO more input parameters
    end

    GLB_store_fobj(t+1) = o;

    % plot objective function output
    figHandles = get(0,'Children');
    if ~isempty(figHandles)
        IdxHandle = strcmp(get(figHandles,'Name'),'Progress of Optimization');
    else
        IdxHandle = [];
    end

    if any(IdxHandle)
        figOpt = figHandles(IdxHandle);
        AxesInfigOpt = findall(figOpt,'type','axes');
        set(AxesInfigOpt,'NextPlot', 'replacechildren')
        children = get(AxesInfigOpt,'children');
        delete(children);
    else
        figOpt = figure('Name','Progress of Optimization','NumberTitle','off','Color',[1 1 1]);
        hold on, grid on, grid minor,
        AxesInfigOpt = findall(figOpt,'type','axes');
    end

    defaultFontSize = 14;
    set(AxesInfigOpt,'YScale','log');
    title(AxesInfigOpt,'Progress of Optimization','LineWidth',defaultFontSize),
    xlabel(AxesInfigOpt,'# iterations','Fontsize',defaultFontSize),ylabel(AxesInfigOpt,'objective function value','Fontsize',defaultFontSize)

    % draw updated axes
    plot(AxesInfigOpt,0:1:t,GLB_store_fobj,'xb','LineWidth',1.5);
    drawnow
end

end

% ----------------------------------------------------------------------
% Compute the value of the CONSTRAINTS
%
% INPUT
% |optFunction| - _struct_ - Structure containing the information about the objective functions set by the user on the dose to the target volume and organs at risk. The following data must be present for each objective function with index i:
% * |optFunction.ROIname| - _string_ - Name of the ROI to which the objective function must be applied.
% * |optFunction.name| - _string_ - Name of the objective function type. The possible types supported so far are: 'min', 'max', 'min_mean' and 'max_mean'.
% * |optFunction.Dref| - _scalar_ -  Reference dose in Gy
% * |optFunction.Vref| - _scalar_ -  Reference volume in 1% only needed for minDVH and maxDVH functions
% * |optFunction.impw| - _scalar_ -  Importance weight for the current objective function (impw > 0)
% * |optFunction.robust| - _scalar_ - Boolean value indicating if the objective function must be considered in all scenarios (if set as robust, equal to 1), or only in the nominal case (if set as non-robust, equal to zero).
%
% |D| -_SCALAR VECTOR_- Dose at each pixel in the whole voume
%
% |i_4D| -_SCALAR_- Index of the scenario for 4D motion
% |rr| -_SCALAR_-  Index of the scenario for random range error
% |rs| -_SCALAR_- Index of the scenario for random shift error
%
% OUTPUT
%
% |fvalue| -_SCALAR_- Value of the largest deviation to the constraint defined in |optFunction|
%
% GLOBAL variable:
% ROI

function fvalue = evalFunction(optFunction,D , i_4D , rr , rs , optFidx)

% This function evaluates the penalty function associated to one ROI for one scenario (specified by the dose D)

ROIindex = optFunction.ROIindex;
d = D(ROI(ROIindex).mask1D); %A 1D vector with as many element as pixels in ROI
fprintf('Obj. %d : Dose in %s : %f <= <D>=%f Gy <= %f \n',optFidx,ROI(ROIindex).name,min(d),mean(d),max(d));

if  (optFunction.ID == 1) % 'min'

    delta = max(0, optFunction.Dref - d); %COmpute the delta between the refernece dose and the dose at the voxel
          %max(0,) is clipping negative value to zero. There will be only positive value in the vector
          %delta is a 1D vector with as many elements as pixels in ROI
    fvalue = (optFunction.impw/ROI(ROIindex).nvoxels) * sum(delta.*delta); % faster
    %The cost is delta.^2 where delta is the  difference between the target dose and the current dose at each voxels in structure
    % The cost is weighted by user defined |optFunction.impw| parameter.
    %The cost is divided by the inverse of the size of the structure: this is the AVERAGE ofthe square difference

elseif (optFunction.ID == 2) % 'max'

    delta = max(0, d - optFunction.Dref);
    fvalue = (optFunction.impw/ROI(ROIindex).nvoxels) * sum(delta.*delta); % faster

elseif (optFunction.ID == 3) % 'min_mean'

    fvalue = optFunction.impw * max(0,(optFunction.Dref - mean(d)));

elseif (optFunction.ID == 4) % 'max_mean'

    fvalue = optFunction.impw * max(0,(mean(d) - optFunction.Dref));

elseif (optFunction.ID == 5) % 'minDVH'

    diff = optFunction.Dref - d;
    dvh_values = sort(d,'descend');
    % compute current dose at Vref
    Dc = dvh_values(max([1 round(optFunction.Vref * ROI(ROIindex).nvoxels)])); %
    diff = diff(d < optFunction.Dref | d > Dc);
    fvalue = (optFunction.impw/ROI(ROIindex).nvoxels) * (diff'*diff);

elseif (optFunction.ID == 6) % 'maxDVH'

     diff = d - optFunction.Dref;
     dvh_values = sort(d,'descend');
     % compute current dose at Vref
     Dc = dvh_values(max([1 round(optFunction.Vref * ROI(ROIindex).nvoxels)])); %
     diff = diff(d > optFunction.Dref | d < Dc);
     fvalue = (optFunction.impw/ROI(ROIindex).nvoxels) * (diff'*diff); % PROBLEM

elseif (optFunction.ID == 7) % 'fall-off'
    % equivalent to a 'max' (ID == 2) objective but with a voxel-wise Dref
    % accounting for a dose gradient in the associated ROI.
    dref = optFunction.Dref(ROI(ROIindex).mask1D);
    delta = max(0, d - dref);
    fvalue = (optFunction.impw/ROI(ROIindex).nvoxels) * sum(delta.*delta); % faster

elseif (optFunction.ID == 8) % 'minimum average of the percentile dose rate constraint'
    fvalue = obj_DR(GLB_DRa{i_4D,rr,rs,optFidx} , optFunction.DRref  , optFunction.impw , ROI(ROIindex).nvoxels);

 elseif (optFunction.ID == 9) % 'minimum median of the percentile dose rate constraint'
      fvalue = obj_DR(GLB_DRm{i_4D,rr,rs,optFidx} , optFunction.DRref  , optFunction.impw , ROI(ROIindex).nvoxels);

 elseif (optFunction.ID == 10) % 'minimum median of the dose averaged dose rate constraint'
      fvalue = obj_DR(GLB_DADRm{i_4D,rr,rs,optFidx} , optFunction.DRref  , optFunction.impw , ROI(ROIindex).nvoxels);

end %if

end %evalFunction


%-------------------------------------------------------
%Evaluate the objective function at given dose rate
%-------------------------------------------------------
function fvalue = obj_DR(dr , DRref , impw , nvoxels)
  delta = max(0, DRref - dr); %Compare to reference dose rate. If dr is higher than constraint, it is OK: delta=0 to remove this contribution from the gof function
  wDR = find(dr > 0); % Find the beams for which a dose rate was computed. Only those must be included in the constraint
  fvalue = (impw/nvoxels)*sum(delta(wDR).*delta(wDR)); %The scalar product sums the constraint for all the beams (for which a dose rate is computed)
end


% ---------------------------------------------------------------------
% Compute the derivative of the CONSTRAINTS
% evalDerivative = d Constraints_{s}(x) / d_pixel intensity (for dose constraints)
% or
% evalDerivative = d Constraints_{s}(x) / d_w (for dose constraints)
%
% Derivative of the constraint function |optFunction| for the s-th scenario (represented by the dose map D)
% with respect to the pixel intensity
%
% INPUT
% |optFunction| - _VECTOR of STRUCT_ - Structure containing the information about the objective functions set by the user on the dose to the target volume and organs at risk. The following data must be present for each objective function with index i:
%
% |D| -_SCALAR VECTOR_- |D(pxl)| Intensity (= dose in Gy) of the pxl-th pixel for the s-th scenario
%
% |optFidx| -_INTEGER_- Index of |optFunction(optFidx)| for which the derivative is computed
%
% |Pij| -_SCALAR MATRIX_- Dose influence matrix. |Pij(i,w)| dose delivered at i-th pixel by w-th spot
%
% |T| -_SCALAR VECTOR_- |T(i)| Time (ms) taken at the i-th pixel to deliver the |percentile| of the dose
%
% OUTPUT
% |fder| -_SCALAR VECTOR_- For dose constraints : |fder(idx)| Derivative of the constraint function with repsect to :
%                 * the intensity (=dose) in the pxl-th pixel: d Constraints_{s}(x) / d_Ipxl. idx = pixel index
%                 * the spot weight : d Constraints_{s}(x) / d_w. idx = weight index
% ---------------------------------------------------------------------

function fder = evalDerivative(optFunction , D , optFidx , Pij , T )

    ROIindex = optFunction.ROIindex;
    d =D(ROI(ROIindex).mask1D); %d is a 1D vector with the number of pixels in the ROI mask


    if (optFunction.ID == 1) % 'min'

        diff = optFunction.Dref - d;
        diff (diff < 0) = 0;
        fder = (optFunction.impw/ROI(ROIindex).nvoxels)*(-2 * diff);

        % See evalFunction : constraint = (impw/ nvoxels) * sum_pxl(Delta .* Delta) = sum_pxl [  (impw/ nvoxels) * Delta .* Delta ]
        %	Delta is a function that depends on the pixel intensity. Therefore : d Constraints_{s}(x) (pxl)/ d_pixel  = (impw/ nvoxels) * 2 .* Delta * d_Delta/ d_pixel
        %	Therefore : d Constraints_{s}(x) / d_pixel  = (impw/ nvoxels) * 2 .* sum(Delta) * d(Dref – d)/ d_pixel  = (impw/ nvoxels) * 2 .*Delta * (-1)


    elseif (optFunction.ID == 2) % 'max'

        diff = d - optFunction.Dref;
        diff(diff < 0) = 0;
        fder = (optFunction.impw/ROI(ROIindex).nvoxels)*(2 * diff);

    elseif (optFunction.ID == 3) % 'min_mean'

        fder = (optFunction.impw/ROI(ROIindex).nvoxels)*ones(ROI(ROIindex).nvoxels,1);

    elseif (optFunction.ID == 4) % 'max_mean'

        fder = (optFunction.impw/ROI(ROIindex).nvoxels)*ones(ROI(ROIindex).nvoxels,1);

    elseif (optFunction.ID == 5) % 'minDVH'

        diff = optFunction.Dref - d;
        dvh_values = sort(d,'descend');
        % compute current dose at Vref
        Dc = dvh_values(max([1 round(optFunction.Vref*ROI(ROIindex).nvoxels)]));
        diff(d > optFunction.Dref | d < Dc) = 0; % need to have the same size as mask1D
        fder = (optFunction.impw/ROI(ROIindex).nvoxels)*(-2* diff);

    elseif (optFunction.ID == 6) % 'maxDVH'

        diff = d - optFunction.Dref;
        dvh_values = sort(d,'descend');
        % compute current dose at Vref
        Dc = dvh_values(max([1 round(optFunction.Vref*ROI(ROIindex).nvoxels)]));
        diff(d > optFunction.Dref | d < Dc) = 0; % need to have the same size as mask1D
        fder = (optFunction.impw/ROI(ROIindex).nvoxels)*(2* diff); %PROBLEM

    elseif (optFunction.ID == 7) % 'fall-off'
        % equivalent to a 'max' (ID == 2) objective but with a voxel-wise Dref
        % accounting for a dose gradient in the associated ROI.
        diff = d - optFunction.Dref(ROI(ROIindex).mask1D);
        diff(diff < 0) = 0;
        fder = (optFunction.impw/ROI(ROIindex).nvoxels)*(2 * diff);

      elseif (optFunction.ID == 8) % 'minimum average of the percentile dose rate constraint'
        Tested =  D > optFunction.Dref; %find the voxels receiving a total dose  > than the threshold and that are inside the ROI. The DR computation will occur only in those pixels
        Tested =  Tested .* ROI(ROIindex).mask1D;
        if ~isempty(find(Tested))
            fder = dC_dw(GLB_DRa{i_4D,rr,rs,optFidx} , optFunction.DRref , optFunction.impw , ROI(ROIindex).nvoxels , Pij(find(Tested),:) , T , D(find(Tested)));
        else
            %There are no pixel meeting the minimum dose requirement
            fder = zeros(1,size(Pij,2));
        end

      elseif (optFunction.ID == 9) % 'minimum median of the percentile dose rate constraint'
        Tested =  D > optFunction.Dref; %find the voxels receiving a total dose  > than the threshold and that are inside the ROI. The DR computation will occur only in those pixels
        Tested =  Tested .* ROI(ROIindex).mask1D;
        if ~isempty(find(Tested))
            fder = dC_dw(GLB_DRm{i_4D,rr,rs,optFidx} , optFunction.DRref , optFunction.impw , ROI(ROIindex).nvoxels , Pij(find(Tested),:) , T , D(find(Tested)));
        else
            %There are no pixel meeting the minimum dose requirement
            fder = zeros(1,size(Pij,2));
        end

      elseif (optFunction.ID == 10) % 'minimum median of the dose averaged dose rate constraint'
        Tested =  D > optFunction.Dref; %find the voxels receiving a total dose  > than the threshold and that are inside the ROI. The DR computation will occur only in those pixels
        Tested =  Tested .* ROI(ROIindex).mask1D;
        if ~isempty(find(Tested))
            fder = dCdadr_dw(GLB_DADRm{i_4D,rr,rs,optFidx} , optFunction.DRref , optFunction.impw , ROI(ROIindex).nvoxels , Pij(find(Tested),:) , T , D(find(Tested)));
        else
            %There are no pixel meeting the minimum dose requirement
            fder = zeros(1,size(Pij,2));
        end

    end

end

%---------------------------------------------------------------------------
%Compute derivative of the constraint function with respect to the spot weight
% d_Constr / d_w
% in the case of the average of the percentile dose rate
%
% INPUT
% |D| -_SCALAR VECTOR_- |D(pxl)| Intensity (= dose in Gy) of the pxl-th pixel for the s-th scenario
%
% |Pij| -_SCALAR MATRIX_- Dose influence matrix. |Pij(i,w)| dose delivered at i-th pixel by w-th spot
%
% |T| -_SCALAR VECTOR_- |T{b}(i)| Time (ms) taken at the i-th pixel to deliver the |percentile| of the dose for the b-th beam
%
%---------------------------------------------------------------------------
function fder = dC_dw(dr , DRref , impw , nvoxels , Pij , T , D)

  Nb_W = size(Pij , 2); %Number of weights
  delta = DRref - dr;
  delta (delta < 0) = 0; % <0 means that the constraint is respected. This term is not included in the objective function

  fder = zeros(1,Nb_W); %fder(w) = sum_beam (d_Const / d_w)

  for b = 1:numel(T)
    %Loop for each beam. Different spots are involved in different beams and
    %the same i-th pixel will have a different T for each beam
    wNZ = find(T{b}); % Find the pixels with non zero T. This will avoid problems with 1/T
    fder = fder - 2 .* (delta(b) .* impw ./ nvoxels)  .*   ( (1./T{b}(wNZ)) * Pij(wNZ,:)     - sum (GLB_alpha .* D(wNZ) ./ T{b}(wNZ)'.^2 , 1) ) ;
    %dC / d_w = sum_beam ()
  end
end

%---------------------------------------------------------------------------
%Compute derivative of the constraint function with respect to the spot weight
% d_Constr / d_w
% in the case of the average of the DADR
%---------------------------------------------------------------------------
function fder = dCdadr_dw(dr , DRref , impw , nvoxels , Pij , T , D)

  Nb_W = size(Pij , 2); %Number of weights
  delta = DRref - dr;
  delta (delta < 0) = 0; % <0 means that the constraint is respected. This term is not included in the objective function

  fder = zeros(1,Nb_W); %fder(w) = sum_beam (d_Const / d_w)

  for b = 1:numel(T)
    %Loop for each beam. Different spots are involved in different beams and
    %the same i-th pixel will have a different T for each beam
    wNZ = find(T{b}); % Find the pixels with non zero T. This will avoid problems with 1/T
    maxD = max(D(wNZ));
    fder = fder - 2 .* (delta(b) .* impw ./ (nvoxels .* maxD))  .*   ( 2.* (D(wNZ)'./T{b}(wNZ)) * Pij(wNZ,:)     - sum (GLB_alpha .* D(wNZ).^2 ./ T{b}(wNZ)'.^2 , 1) ) ;
    %dC / d_w = sum_beam ()
  end
end

%---------------------------------------------------------------------------
%Check whether the speicfied objective function requires the computation of dose rate
%---------------------------------------------------------------------------
function value = UseDoseRate(ID)
  value = ID == 8 || ID == 9 || ID == 10;
end

% ---------------------------------------------------------------------

% ----------------------------------------------------------------------
end %function [x, info , spotSequence] = NLPsolver_double_matlabNative2
%All the other functions are declared inside the body of NLPsolver_double_matlabNative2 so that it is possible to use the global variable
