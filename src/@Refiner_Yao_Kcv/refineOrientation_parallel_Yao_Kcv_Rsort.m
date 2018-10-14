
function obj = refineOrientation_parallel_Yao_Kcv_Rsort(obj)

k_fold_value = obj.k_fold_value;

% if length(unique(size(obj.refineModel))) ~= 1
%     error('Model should be a cubic array')
% else
%     [dim1,dim2,~] = size(obj.refineModel);
% end

% pad model for projection calculation
%padding = round(dim*(obj.oversampling_ratio-1)/2);
%modelK = my_fft(padarray(obj.refineModel,[padding, padding, padding]));

% in case of using real-space forward projector
% if obj.RealProjection == 1
%   modelK = obj.refineModel;
% else
%   modelK = my_fft(My_paddzero(obj.refineModel,round(size(obj.refineModel)*obj.oversampling_ratio)));
% end

num_proj = size(obj.refineProjections,3);
num_calculations = length(obj.phi_search_range) * length(obj.theta_search_range) * length(obj.psi_search_range);

metrics       = zeros(num_proj,num_calculations);
phis          = zeros(num_proj,num_calculations);
thetas        = zeros(num_proj,num_calculations);
psis          = zeros(num_proj,num_calculations);
x_centers     = zeros(num_proj,num_calculations);
y_centers     = zeros(num_proj,num_calculations);
bayes_probs   = zeros(num_proj,num_calculations);

if strcmp(obj.evaluator_type,'backprojection')
    EVALUATOR     = Backprojection_Evaluator(obj);
elseif strcmp(obj.evaluator_type,'ccl')
    EVALUATOR     = CCL_Evaluator(obj);
else
    error('Unrecognized evaluator type %s',evaluator_type)
end

% separate arrays to avoid unnecessary broadcast variables
projections            = obj.refineProjections;
phiangles              = obj.refineAngles(:,1);
thetaangles            = obj.refineAngles(:,2);
psiangles              = obj.refineAngles(:,3);
phi_search_range       = obj.phi_search_range;
theta_search_range     = obj.theta_search_range;
psi_search_range       = obj.psi_search_range;
compare_func           = @obj.compare_func;
forwardProjection_func = @obj.forwardProjection_func;
FPmask                 = obj.FPmask;

% cropInd1 =  (1:dim1)-ncx + nc_padded1;
% cropInd2 =  (1:dim2)-ncy + nc_padded2;

window_half_size = obj.window_half_size;
Rscanres = obj.Rscanres;
Rmethod = obj.Rmethod;

EVALUATOR = EVALUATOR.prepData(); % todo: implement for backprojection type classes

raw_Rec = obj.refineModel;
pj = obj.refineFullProjections;
euler_angles = obj.refineAngles;

[~,Rarray,~]=Tian_calc_Rfactor_realspace(raw_Rec,pj,euler_angles);
[~,Large_R_arr_num] = sort(Rarray,'descend');
obj.R_arr = Large_R_arr_num;

if obj.FullEvolutionRecord==1
  obj.R_evolution(:,:,end+1) = Rarray;
end

x_centers_parfor    = zeros(k_fold_value,num_calculations);
y_centers_parfor    = zeros(k_fold_value,num_calculations);
metrics_parfor      = zeros(k_fold_value,num_calculations);
phis_parfor         = zeros(k_fold_value,num_calculations);
thetas_parfor       = zeros(k_fold_value,num_calculations);
psis_parfor         = zeros(k_fold_value,num_calculations);
bayes_probs_parfor  = zeros(k_fold_value,num_calculations);


parfor pj_ind = 1:k_fold_value
    pj_num = Large_R_arr_num(pj_ind);
    fprintf('Refining projection #%d/%d\n',pj_num,num_proj)
%     pj            = obj.refineProjections(:,:,pj_num);
    phi           = phiangles(pj_num);
    theta         = thetaangles(pj_num);
    psi           = psiangles(pj_num);
%     noise_sigma   = obj.noise_sigma(:,:,pj_num);
    calc_count    = 1;
    
    FSind = find(FPmask(:,:,pj_num));
    
    if Rmethod == 1
        VinCell = {window_half_size,Rscanres,FSind};
    else
        VinCell = {};
    end
    
    
    
    % because of MATLAB's for rules have to create separate variables
    % for storing results and collect them later
    tmp_phis        = zeros(1,num_calculations);
    tmp_thetas      = zeros(1,num_calculations);
    tmp_psis        = zeros(1,num_calculations);
    tmp_metrics     = zeros(1,num_calculations);
    tmp_centers_x   = zeros(1,num_calculations);
    tmp_centers_y   = zeros(1,num_calculations);
    tmp_bayes_probs = zeros(1,num_calculations);

    
    for cur_phi = (phi+phi_search_range)
            for cur_theta = (theta+theta_search_range)
                    for cur_psi = (psi+psi_search_range)
                        tmp_phis(1,calc_count)    = cur_phi;
                        tmp_thetas(1,calc_count)  = cur_theta;
                        tmp_psis(1,calc_count)    = cur_psi;
                        [metric, suggested_center_x, suggested_center_y] = EVALUATOR.evaluateOrientation(pj_num,cur_phi,cur_theta,cur_psi);
                        tmp_metrics(1,calc_count)   = metric;
                        tmp_centers_x(1,calc_count) = suggested_center_x;
                        tmp_centers_y(1,calc_count) = suggested_center_y;
                        calc_count = calc_count + 1;
                    end
            end
    end
    
    x_centers_parfor(pj_ind,:)        = tmp_centers_x;
    y_centers_parfor(pj_ind,:)        = tmp_centers_y;
    metrics_parfor(pj_ind,:)          = tmp_metrics;
    phis_parfor(pj_ind,:)             = tmp_phis;
    thetas_parfor(pj_ind,:)           = tmp_thetas;
    psis_parfor(pj_ind,:)             = tmp_psis;  
    bayes_probs_parfor(pj_ind,:)      = tmp_bayes_probs; 
end

pj_k_set = Large_R_arr_num(1:k_fold_value);

x_centers(pj_k_set,:)        = x_centers_parfor;
y_centers(pj_k_set,:)        = y_centers_parfor;
metrics(pj_k_set,:)          = metrics_parfor;
phis(pj_k_set,:)             = phis_parfor;
thetas(pj_k_set,:)           = thetas_parfor;
psis(pj_k_set,:)             = psis_parfor;
bayes_probs(pj_k_set,:)      = bayes_probs_parfor;

% collect results
obj.x_centers   = x_centers;
obj.y_centers   = y_centers;
obj.metrics     = metrics;
obj.phis        = phis;
obj.thetas      = thetas;
obj.psis        = psis;
obj.bayes_probs = bayes_probs;
end
