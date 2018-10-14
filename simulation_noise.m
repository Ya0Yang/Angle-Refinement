close all
clc
clear
% This script calculates projections and reconstructs with perturbed
% angles which are then recovered by refinement
tic
addpath([(pwd) '/src/'])
addpath([(pwd) '/src/accelerate'])
% model_filename   = 'models/vesicle.mat';
model_filename   = 'models/FePt_Vol_64_asymmetric.mat';
pj_filename      = 'data/projections.mat';
angle_filename   = 'data/angles.mat';
support_filename = 'data/support.mat';

if ~isdir('data')
    mkdir data
end
if ~isdir('results')
    mkdir 
end

random_seed     = 46; % for reproducibility
rng(random_seed);
MW_size         = 40;
tomo_angles     = round(MW_size/2):2:(180-round(MW_size/2)); % true tilt angles
tomo_angles(2:end)     = tomo_angles(2:end) + rand(1,numel(tomo_angles(2:end)));
true_phis       = (rand(1,numel(tomo_angles))*360);true_phis(1)=0;
true_psis       = (rand(1,numel(tomo_angles))*360);true_psis(1)=0;
OR              = 2; % oversampling ratio used for padding model in 3D for projection calculation
jitter          = 4; % euler angles will be perturbed by up to += jitter/2
jitter_translate= 2; % x and y translations will be perturbed by +- jitter_translate/2
noise_level     = 1;
poiss_scaling   = 6.5; % ~5% poisson noise for FePt model

%% calculate projections
model       = importdata(model_filename); 
model(end+1,end+1,end+1) = 0;
[dim,~,~]   = size(model);
nc          = round((dim+1)/2); n2 = nc-1;
padding     = round((OR-1)*dim/2);
modelK      = my_fft(padarray(model,[padding, padding, padding]));
nc_padded   = round((size(modelK,1)+1)/2);
num_proj    = length(tomo_angles);
pj          = zeros(dim,dim,num_proj);
noise_sigma = zeros(size(modelK,1),size(modelK,1),num_proj);
%%

central_vec = (1:dim) - nc + nc_padded;
true_shifts = round(rand(num_proj,2) * jitter_translate) -jitter_translate/2;
for pj_num = 1:num_proj
    if mod(pj_num,1)==0
        fprintf('Calculating projection #%d/%d\n',pj_num,num_proj)
    end
   [tmp, kSlice]    = calculate3Dprojection_interp(modelK,true_phis(pj_num),tomo_angles(pj_num),true_psis(pj_num)); 
   pj(:,:,pj_num)   = tmp(central_vec,central_vec);
   pj(:,:,pj_num)   = circshift(pj(:,:,pj_num),[true_shifts(pj_num,1), true_shifts(pj_num,2)]);
   tmp = poissrnd(poiss_scaling*pj(:,:,pj_num))/poiss_scaling;
   tmp(isnan(tmp)) = 0;
   pj(:,:,pj_num)   = tmp + randn(size(pj(:,:,pj_num)))*noise_level;
   noise_sigma(:,:,pj_num) = smooth3D(abs(kSlice/10),.25); % this is just for testing purposes now, need to figure out real SNR later
end
noise_sigma(noise_sigma==0) = 1e-2;

%% perturb angles
true_angles = [true_phis',tomo_angles',true_psis'];
% euler_angles = zeros(num_proj,3);
euler_angles = true_angles;

euler_angles(:,1) = true_phis + (jitter*rand(length(tomo_angles),1)' - jitter/2);
euler_angles(:,2) = tomo_angles + (jitter*rand(length(tomo_angles),1)' - jitter/2);
euler_angles(:,3) = true_psis + (jitter*rand(length(tomo_angles),1)' - jitter/2);
starting_angles = euler_angles;
% euler_angles(3:1:end,:) = true_angles(3:1:end,:);
support = ones(dim,dim,dim);

%%
save(pj_filename,'pj')
save(angle_filename,'euler_angles')
save(support_filename,'support')
REFINEMENT_Main
