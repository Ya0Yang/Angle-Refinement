t1 = clock;

addpath('src/')
addpath('splinterp/')
addpath('nufft2d-1.3.2/')

pj_filename              =  'data_turb3/pj_new_Turbu3_1015.mat';
angle_filename           =  'data_turb3/angles_usedin_refine_1015.mat';

recon_filename           =  'results/refinement_results.mat';
%% GENFIRE parameters

OR = 2;  %oversampling_ratio for backprojection
num_refinements = 5;
FFTorDFT = 2;
doThetaOnly=0;
vector3 = [1 0 0]; 

DOGPU = 0; 
%%
REFINEMENT = Refiner('GENFIRE');

% set REFINEMENT parameters
REFINEMENT = REFINEMENT.set_parameters('ang_search_range',-1:0.5:1,...
    'compare_func', @alignByNormXCorr, ... 
    'forwardProjection_func',@calculate3Dprojection_interp_wrapperZYX, ...
    'oversampling_ratio',OR,... 
    'maximize',true, 'use_parallel',0, ...
    'FullEvolutionRecord',1, ...
    'num_refinements',num_refinements, ...
    'evaluator_type','backprojection','RefineZeroCenterFlag',1);

REFINEMENT.RECONSTRUCTOR.vector3 = vector3; 

REFINEMENT.RECONSTRUCTOR.filename_Projections = pj_filename;
REFINEMENT.RECONSTRUCTOR.filename_Angles = angle_filename;

REFINEMENT.RECONSTRUCTOR.numIterations = 50; 
REFINEMENT.RECONSTRUCTOR.monitorR_loopLength = REFINEMENT.RECONSTRUCTOR.numIterations;
REFINEMENT.RECONSTRUCTOR.monitor_R = 1;

REFINEMENT.RECONSTRUCTOR.oversamplingRatio = 3;
REFINEMENT.RECONSTRUCTOR.griddingMethod = FFTorDFT;  %1 FFT; 2DFT
REFINEMENT.RECONSTRUCTOR.constraintEnforcementMode = 3; 
REFINEMENT.RECONSTRUCTOR.interpolationCutoffDistance =.2; 

REFINEMENT.RECONSTRUCTOR.DFT_doGPU = DOGPU;
%%
% initilize RECONSTRUCTOR, read data, and check data validity
REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.readFiles();
REFINEMENT.RECONSTRUCTOR = REFINEMENT.RECONSTRUCTOR.CheckPrepareData();

% import projections and angles from RECONSTRUCTOR to REFINEMENT
REFINEMENT = REFINEMENT.get_projections_from_RECONSTRUCTOR();  
REFINEMENT = REFINEMENT.get_angles_from_RECONSTRUCTOR();
REFINEMENT.measuredProjections = REFINEMENT.RECONSTRUCTOR.InputProjections;
% run refinement

REFINEMENT = REFINEMENT.refineControl();


t2 = clock;
fprintf('Completed in %.02f seconds\n',etime(t2,t1))

% get refined result
save(recon_filename, 'REFINEMENT');