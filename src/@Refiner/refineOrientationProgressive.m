
function obj = refineOrientationProgressive(obj,currProjsNum,batchRefinement)

num_proj = batchRefinement;
num_calculations = length(obj.phi_search_range) * length(obj.theta_search_range) * length(obj.psi_search_range);

num_fulProjs = size(obj.refineFullProjections,3);

metrics       = zeros(num_fulProjs,num_calculations);
phis          = zeros(num_fulProjs,num_calculations);
thetas        = zeros(num_fulProjs,num_calculations);
psis          = zeros(num_fulProjs,num_calculations);
x_centers     = zeros(num_fulProjs,num_calculations);
y_centers     = zeros(num_fulProjs,num_calculations);
bayes_probs   = zeros(num_fulProjs,num_calculations);

if strcmp(obj.evaluator_type,'backprojection')
    EVALUATOR     = Backprojection_Evaluator(obj);
elseif strcmp(obj.evaluator_type,'ccl')
    EVALUATOR     = CCL_Evaluator(obj);
else
    error('Unrecognized evaluator type %s',evaluator_type)
end

% separate arrays to avoid unnecessary broadcast variables
phiangles              = obj.refineAngles(:,1);
thetaangles            = obj.refineAngles(:,2);
psiangles              = obj.refineAngles(:,3);
phi_search_range       = obj.phi_search_range;
theta_search_range     = obj.theta_search_range;
psi_search_range       = obj.psi_search_range;


EVALUATOR = EVALUATOR.prepData(); % todo: implement for backprojection type classes

parfor pj_num = currProjsNum+1:currProjsNum+num_proj  %be very careful: now pj_num is not the real projection number but the number in the sequence
    
    fprintf('Refining projection #%d\n',pj_num)
    
    phi           = phiangles(pj_num);
    theta         = thetaangles(pj_num);
    psi           = psiangles(pj_num);
    
    calc_count    = 1;

    
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
    
    x_centers(pj_num,:)        = tmp_centers_x;
    y_centers(pj_num,:)        = tmp_centers_y;
    metrics(pj_num,:)          = tmp_metrics;
    phis(pj_num,:)             = tmp_phis;
    thetas(pj_num,:)           = tmp_thetas;
    psis(pj_num,:)             = tmp_psis;  
    bayes_probs(pj_num,:)      = tmp_bayes_probs; 
end
% collect results
obj.x_centers   = x_centers;
obj.y_centers   = y_centers;
obj.metrics     = metrics;
obj.phis        = phis;
obj.thetas      = thetas;
obj.psis        = psis;
obj.bayes_probs = bayes_probs;
end
