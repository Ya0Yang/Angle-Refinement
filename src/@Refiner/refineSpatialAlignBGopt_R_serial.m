%%  refineSpatialAlign_R %%

%% realign the spatial alignment of projections

%% Author: Yongsoo Yang
%% Jianwei (John) Miao Coherent Imaging Group
%% University of California, Los Angeles
%% Copyright (c) 2015-2016. All Rights Reserved.

function obj = refineSpatialAlignBGopt_R_serial(obj)
[dimx,dimy,num_proj] = size(obj.refineFullProjections);

refineAngles  = obj.refineAngles;
refineFullProjs = obj.refineFullProjections;
refineProjs = obj.refineProjections;

realignedProjs = zeros(size(refineFullProjs));

refineModelK = My_FFTN(My_paddzero(obj.refineModel,size(obj.refineModel)*obj.oversampling_ratio));

RscanRad = obj.RscanRad;
Rres = obj.Rscanres_separateR;
bin_factor = obj.bin_factor;
ShiftAr = zeros(num_proj,2);
BGforAddAr = zeros(num_proj,1);

for pj_num = 1:num_proj
    fprintf('Spatial aligning projection #%d/%d\n',pj_num,num_proj)

    phi           = refineAngles(pj_num,1);
    theta         = refineAngles(pj_num,2);
    psi           = refineAngles(pj_num,3);
    
    img1 = refineProjs(:,:,pj_num);
    [img2, ~] =calculate3Dprojection_interp(refineModelK,phi,theta,psi);
    img2 = My_stripzero( img2,size(img1));
    [d1, d2, BGforAdd, ~] = My_Rscan_subpixel_align_img1mask_bgFit(img1,img2,RscanRad,Rres)

    %shiftX = (ncx_binned - best_center_x)*bin_factor;
    %shiftY = (ncy_binned - best_center_y)*bin_factor;
    
    if Rscanres < 1
        realignedProjs(:,:,pj_num) = real(My_FourierShift(refineFullProjs(:,:,pj_num)+BGforSubtract*(refineFullProjs(:,:,pj_num)~=0),-1*d1*bin_factor, -1*d2*bin_factor));
    else
        realignedProjs(:,:,pj_num) = circshift(refineFullProjs(:,:,pj_num)+BGforSubtract*(refineFullProjs(:,:,pj_num)~=0),[-1*d1*bin_factor, -1*d2*bin_factor]);
    end
    ShiftAr(pj_num,:) = [-1*d1*bin_factor, -1*d2*bin_factor];
    BGforAddAr(pj_num) = BGforAdd;
end

obj.refineFullProjections = realignedProjs;


if obj.FullEvolutionRecord==1
  
  obj.ShiftEvolution_Rscan(:,:,end+1) = ShiftAr;
  obj.ShiftEvolution_BGopt(:,end+1) =  BGforAddAr;
end


end