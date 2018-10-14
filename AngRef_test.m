%addpath ./src/
tic

pj_filename              = 'data/data_Turbu4_1129/projections.mat';
angle_filename           = 'data/data_Turbu4_1129/angles.mat';
filenameFinalRefine      = 'results/AngRef_res';
OR                       = 2; % oversampling ratio
evaluator_type           = 'backprojection'; % backprojection or ccl
num_refinements          = 10; % number of loops to run
angle_search_range       = -3:.2:3; % phi, theta, psi will be searched around current one with this search range
translation_search_range = -3:1:3; % only used for ccl -- backprojection searches all at once
use_parallel             = true; % use parfor where possible
compare_func             = @alignByNormXCorr; % evaluation function for each projection
maximize                 = true; % is the error metric returned by compare_func to be maximized or minimized?
save_intermediates       = true;
filename_intermediates   = 'results/refinement_results';
reference_angle          = [0 0 0]; % reconstruction will be reoriented so projection 1 has these euler angles. This is mostly useful in simulations
% k_fold_value             = 4;

%% run refinement

FFTorDFT = 2;
doThetaOnly=0;
vector3 = [1 0 0]; 

DOGPU = 1;
%%
% initilize REFINEMENT object with RECONSTRUCTOR, in this case with GENFIRE

% REFINEMENT = Refiner_Yao_Kcv('GENFIRE');
REFINEMENT = Refiner_YY('GENFIRE');

% set REFINEMENT parameters
REFINEMENT = REFINEMENT.set_parameters('ang_search_range',-2:0.4:2,...
    'compare_func', @alignByNormXCorr, ... 
    'forwardProjection_func',@calculate3Dprojection_interp_wrapperZYX, ...
    'oversampling_ratio',OR, ... 
    'maximize',true,'use_parallel',true, ...
    'FullEvolutionRecord',1, ...
    'num_refinements',num_refinements, ...
    'evaluator_type','backprojection','RefineZeroCenterFlag',1);
% REFINEMENT = REFINEMENT.set_parameters('FullEvolutionRecord',1);
% REFINEMENT = REFINEMENT.set_parameters('noise_sigma',noise_sigma);
% REFINEMENT = REFINEMENT.set_parameters('num_refinements',num_refinements);

% set RECONSTRUCTOR parameters (in this case GENFIRE)
REFINEMENT.RECONSTRUCTOR.filename_Projections = pj_filename;
REFINEMENT.RECONSTRUCTOR.filename_Angles = angle_filename;

REFINEMENT.RECONSTRUCTOR.numIterations = 50; 
REFINEMENT.RECONSTRUCTOR.pixelSize = .5; 
REFINEMENT.RECONSTRUCTOR.oversamplingRatio = 2;
REFINEMENT.RECONSTRUCTOR.griddingMethod = 2; 
REFINEMENT.RECONSTRUCTOR.allowMultipleGridMatches = 1;
REFINEMENT.RECONSTRUCTOR.constraintEnforcementMode = 3; 
REFINEMENT.RECONSTRUCTOR.interpolationCutoffDistance =.2; 
REFINEMENT.RECONSTRUCTOR.constraintPositivity = 1;
REFINEMENT.RECONSTRUCTOR.constraintSupport = 1;
REFINEMENT.RECONSTRUCTOR.ComputeFourierShellCorrelation = 0; 
REFINEMENT.RECONSTRUCTOR.numBins = 50;
REFINEMENT.RECONSTRUCTOR.percentValuesForRfree = 0.00;
REFINEMENT.RECONSTRUCTOR.numBinsRfree = 35;
REFINEMENT.RECONSTRUCTOR.doCTFcorrection = 0;
REFINEMENT.RECONSTRUCTOR.CTFThrowOutThreshhold = 0;
REFINEMENT.RECONSTRUCTOR.DFT_doGPU = DOGPU;


% initilize RECONSTRUCTOR, read data, and check data validity
REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.readFiles();
REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.CheckPrepareData();

% import projections and angles from RECONSTRUCTOR to REFINEMENT
REFINEMENT = REFINEMENT.get_projections_from_RECONSTRUCTOR();  
REFINEMENT = REFINEMENT.get_angles_from_RECONSTRUCTOR();


% run refinement
REFINEMENT = REFINEMENT.refineControl();

refined_angles = REFINEMENT.AngleEvolution(:,:,end);
refined_projs = REFINEMENT.refineFullProjections;

AngleEvolution = REFINEMENT.AngleEvolution;
ShiftEvolution = REFINEMENT.ShiftEvolution;

final_Rec = REFINEMENT.RECONSTRUCTOR.reconstruction;
final_errK = REFINEMENT.RECONSTRUCTOR.errK;

save(filenameFinalRefine,'refined_angles','refined_projs', 'final_Rec','final_errK',...
    'AngleEvolution','ShiftEvolution');

% [Rfactor,Rarray,~]=Tian_calc_Rfactor_realspace...
%     (REFINEMENT.RECONSTRUCTOR.reconstruction,REFINEMENT.RECONSTRUCTOR.InputProjections,REFINEMENT.RECONSTRUCTOR.InputAngles);
toc
