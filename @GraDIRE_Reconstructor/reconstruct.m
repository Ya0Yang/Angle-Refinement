%% GraDIRE reconstruction

function obj = reconstruct(obj)
projections = obj.InputProjections;
Num_pj      = obj.NumProjs;
step_size   = obj.step_size;
dimx        = obj.Dim1;
dimy        = obj.Dim2;
dtype       = obj.dtype;
Rot_x       = obj.Rot_x;
Rot_y       = obj.Rot_y;
iterations  = obj.numIterations;

sigma       = obj.sigma;
if sigma
    kernel  = obj.kernel;
end

xj = obj.xj;
yj = obj.yj;
zj = obj.zj;

sum_rot_pjs = obj.sum_rot_pjs; 

rec = zeros(dimy,dimx,dimy,dtype);
rec_big = zeros(obj.n2_oversampled,obj.n1_oversampled,obj.n2_oversampled,dtype);
ind_V	= My_volumn_index(size(rec_big),size(rec));
if obj.initial_model == 1
    rec = obj.Support;
    rec_big(ind_V(1,1):ind_V(1,2), ind_V(2,1):ind_V(2,2), ind_V(3,1):ind_V(3,2) ) = rec;
end
shape_limit = obj.shape_limit;
if shape_limit == 1
    shape_constrain = obj.shape_constrain;
    rec(~shape_constrain) = 0;
    rec_big(ind_V(1,1):ind_V(1,2), ind_V(2,1):ind_V(2,2), ind_V(3,1):ind_V(3,2) ) = rec;
end

dt      = (step_size/Num_pj/dimx);

fprintf('GraDIRE: Reconstructing... \n\n');

if obj.monitor_R == 1
    monitorR_loopLength = obj.monitorR_loopLength;
    errR_arr = zeros(1,floor(iterations./monitorR_loopLength));
    Rarr_record = zeros(Num_pj,floor(iterations./monitorR_loopLength));
    Rarr2_record = zeros(Num_pj,floor(iterations./monitorR_loopLength));
end

if obj.use_parallel
  parforArg = Inf;
else
  parforArg = 0;
end

% parfor (pj_num = 1:num_proj, parforArg)
   
for iter=1:iterations
    recK  = double(my_fft(rec_big));    
    
    % smoothing
    if sigma
        recK = recK.*kernel;
        rec  = real(my_ifft(recK));
        rec = croppedOut(rec, [dimx, dimy,dimx]);        
    end    
    
    % compute rotated projections via FST
    pj_cal = splinterp3(recK, yj, xj, zj);
    pj_cal = real(fftshift(ifft2(ifftshift(pj_cal))));
    pj_cal = croppedOut(pj_cal, [dimx,dimy,Num_pj] );
    
    % compute R factor
    if obj.monitor_R == 1 && mod(iter,monitorR_loopLength) == 0
        Rarr = zeros(Num_pj,1);
        Rarr2 = zeros(Num_pj,1);
        for i=1:Num_pj
            pj = projections(:,:,i); 
            proj_i = pj_cal(:,:,i);
            Rarr(i) = sum(sum( abs(proj_i - pj) ))/ sum(abs(pj(:)));
            Rarr2(i) = norm( proj_i - pj, 'fro' )/ norm( pj, 'fro' );
        end
        errR  = mean(Rarr);
        errR2 = mean(Rarr2);
        fprintf('GraDIRE: Iteration %d. Rfactor=%.4f, R2factor=%.4f \n',iter, errR, errR2);
        errR_arr(iter./monitorR_loopLength) = errR;
        Rarr_record(:,iter./monitorR_loopLength) = Rarr;
        Rarr2_record(:,iter./monitorR_loopLength) = Rarr2;
    else
        fprintf('GraDIRE: Iteration %d \n',iter);
    end    
    
    % compute gradient & apply gradient descent
    grad = -sum_rot_pjs;
    parfor (k = 1:Num_pj,parforArg)
%     for k = 1:Num_pj
        rot_pj_cal  = splinterp2(pj_cal(:,:,k), double(Rot_y(:,:,:,k)), double(Rot_x(:,:,:,k)));
        grad = grad + rot_pj_cal;
    end
    rec = rec - dt*grad;
    rec = max(0,rec);
    if shape_limit == 1
        rec(~shape_constrain) = 0;
    end
    
    if obj.save_temp == 1 && mod(iter,obj.save_loopLength) == 0
        temp_filename = [obj.saveFilename,num2str(iter),'.mat'];
        save(temp_filename, 'rec');
    end
    
    rec_big( ind_V(1,1):ind_V(1,2), ind_V(2,1):ind_V(2,2), ind_V(3,1):ind_V(3,2) ) = rec;
end
if obj.monitor_R == 1
    obj.errR            = errR;
    obj.Rarr_record     = Rarr_record;
    obj.Rarr2_record    = Rarr2_record;
end
obj.reconstruction = rec;
end