tic
%addpath ./src/
%addpath ./src/accelerate

%refinement parameters

pj_filename              = 'data/projections.mat';
angle_filename           = 'data/angles.mat';
support_filename         = 'data/support.mat';
OR                       = 3; % oversampling ratio
evaluator_type           = 'backprojection'; % backprojection or ccl
num_refinements          = 5; % number of loops to run
angle_search_range       = -3:.2:3; % phi, theta, psi will be searched around current one with this search range
translation_search_range = -3:1:3; % only used for ccl -- backprojection searches all at once
use_parallel             = true; % use parfor where possible
compare_func             = @alignByNormXCorr; % evaluation function for each projection
maximize                 = true; % is the error metric returned by compare_func to be maximized or minimized?
save_intermediates       = true;
filename_intermediates   = 'results/refinement_results';
reference_angle          = [0 0 0]; % reconstruction will be reoriented so projection 1 has these euler angles. This is mostly useful in simulations
%% GENFIRE parameters

numIterations = 100; 
oversamplingRatio = 3;
griddingMethod = 1; 
allowMultipleGridMatches = 1;
constraintEnforcementMode = 1; 
interpolationCutoffDistance =.7; 
constraintPositivity = 1;
constraintSupport = 1;
bin_factor = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run refinement

% initilize REFINEMENT object with RECONSTRUCTOR, in this case with GENFIRE
REFINEMENT = Refiner('GENFIRE');

% set REFINEMENT parameters
REFINEMENT = REFINEMENT.set_parameters('ang_search_range',angle_search_range,...
    'translation_search_range',translation_search_range, ...
    'compare_func', compare_func, ... 
    'oversampling_ratio',OR, ... 
    'maximize',maximize,'use_parallel',use_parallel, ...
    'FullEvolutionRecord',1, ...
    'num_refinements',num_refinements, ...
    'evaluator_type',evaluator_type, ...
    'save_intermediates',save_intermediates, ...
    'filename_intermediate_results',filename_intermediates, ...
    'bin_factor',bin_factor, ...
    'RefineReferenceAngletoSet',reference_angle);

% set RECONSTRUCTOR parameters (in this case GENFIRE)
REFINEMENT.RECONSTRUCTOR.filename_Projections = pj_filename;
REFINEMENT.RECONSTRUCTOR.filename_Angles = angle_filename;
REFINEMENT.RECONSTRUCTOR.filename_Support = support_filename; 
REFINEMENT.RECONSTRUCTOR.filename_Results = './results/REFINEMENT_finalREconstruction.mat';
REFINEMENT.RECONSTRUCTOR.numIterations = numIterations; 
REFINEMENT.RECONSTRUCTOR.oversamplingRatio = oversamplingRatio;
REFINEMENT.RECONSTRUCTOR.griddingMethod = griddingMethod; 
REFINEMENT.RECONSTRUCTOR.allowMultipleGridMatches = allowMultipleGridMatches;
REFINEMENT.RECONSTRUCTOR.constraintEnforcementMode = constraintEnforcementMode; 
REFINEMENT.RECONSTRUCTOR.interpolationCutoffDistance = interpolationCutoffDistance; 
REFINEMENT.RECONSTRUCTOR.constraintPositivity = constraintPositivity;
REFINEMENT.RECONSTRUCTOR.constraintSupport = constraintSupport;


% initilize RECONSTRUCTOR, read data, and check data validity
REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.readFiles();
REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.CheckPrepareData();

% import projections and angles from RECONSTRUCTOR to REFINEMENT
REFINEMENT = REFINEMENT.get_projections_from_RECONSTRUCTOR();  
REFINEMENT = REFINEMENT.get_angles_from_RECONSTRUCTOR();

% run refinement
REFINEMENT = REFINEMENT.refineControl();


% get refined result
refined_angles = REFINEMENT.refineAngles;


save('REFINERESULT_new.mat','refined_angles','REFINEMENT');


hold off
sprintf('Completed in %d seconds',toc)
