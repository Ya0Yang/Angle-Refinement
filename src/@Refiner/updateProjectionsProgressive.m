%%  updateProjections %%

%%updates input projections and angles based upon results of refinement
%%inputs:
%%  original_projections - original N x N x num_projections array that produced the model
%%  refinement_pars      - output of refineOrientation

%%outputs:
%%  projections   - new projections
%%  euler_angles  - new euler angles

%% Author: Alan (AJ) Pryor, Jr.
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015-2016. All Rights Reserved.

function obj = updateProjectionsProgressive(obj,currProjsNum,batchRefinement)

projectionSeq=obj.projectionSeq;

num_proj = batchRefinement;
num_fulProjs = size(obj.refineFullProjections,3);

bin_factor = obj.bin_factor;
euler_angles = zeros(num_fulProjs,3);

original_projections = obj.refineFullProjections;

[dimx_binned,dimy_binned,~] = size(obj.refineProjections);
ncx_binned = round((dimx_binned+1)/2);
ncy_binned = round((dimy_binned+1)/2);

ShiftAr = zeros(num_fulProjs,2);

for proj_num    = currProjsNum+1:currProjsNum+num_proj  %be very careful: now pj_num is not the real projection number but the number in the sequence
    
    centers_x   = obj.x_centers(proj_num,:);
    centers_y   = obj.y_centers(proj_num,:);
    metrics     = obj.metrics(proj_num,:);
    phis        = obj.phis(proj_num,:);
    thetas      = obj.thetas(proj_num,:);
    psis        = obj.psis(proj_num,:);
    bayes_probs = obj.bayes_probs(proj_num,:);
    if obj.maximize
        best_ind = find(metrics==max(metrics(:)),1);
    else
        best_ind = find(metrics==min(metrics(:)),1);
    end

    best_center_x = centers_x(best_ind);
    best_center_y = centers_y(best_ind);
    euler_angles(proj_num,1) = phis(best_ind);
    euler_angles(proj_num,2) = thetas(best_ind);
    euler_angles(proj_num,3) = psis(best_ind);
    
    
    shiftX = (ncx_binned - best_center_x)*bin_factor;
    shiftY = (ncy_binned - best_center_y)*bin_factor;
    
    if obj.Rscanres < 1
        original_projections(:,:,proj_num) = real(My_FourierShift(original_projections(:,:,proj_num),shiftX, shiftY));
    else
        original_projections(:,:,proj_num) = circshift(original_projections(:,:,proj_num),[shiftX, shiftY]);
    end
    ShiftAr(proj_num,:) = [shiftX, shiftY];
    obj.refineFullProjections(:,:,proj_num) = original_projections(:,:,proj_num); %only update the selective projections and angles
    obj.refineAngles(proj_num,:) = euler_angles(proj_num,:);
    
end
% obj.refineFullProjections = original_projections;
% obj.refineAngles = euler_angles;

if obj.FullEvolutionRecord==1
  obj.AngleEvolution(projectionSeq(1:num_fulProjs),:,end+1) = obj.refineAngles;
  obj.ShiftEvolution(projectionSeq(1:num_fulProjs),:,end+1) = ShiftAr;
end