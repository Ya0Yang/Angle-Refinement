if ~isdeployed
    addpath ./src/
    addpath ./data/
end

savePreFix = 'CMlineRef_FePt64_10proj_AllRand2_B2';

model_filename   = './models/vesicle.mat';
pj_filename      = sprintf('./data/RealPJ_%s.mat',savePreFix);
angle_filename   = sprintf('./data/Angles_%s.mat',savePreFix);
support_filename = sprintf('./data/Fsupport_%s.mat',savePreFix);

Model = importdata(model_filename);

% % bin the model
% model_b = zeros(32,32,32);
% for i=1:2
%   for j=1:2
%     for k=1:2
%       model_b = model_b + Model(i:2:64,j:2:64,k:2:64);
%     end
%   end
% end
% 
% %Model = model_b;

%% Euler angles for initial positions

% initialize angles
random_seed = 2017;
numProjs = 10;
inpMaxError = 2;

rng(random_seed);
The = linspace(-70,70,numProjs)';
Phi = zeros(size(The));
Psi = zeros(size(The));
% Phi = (rand(numProjs,1) - 0.5)*2*inpMaxError;
% Psi = (rand(numProjs,1) - 0.5)*2*inpMaxError;

Angles = [Phi The Psi];

% made 9th angle phi to 40 for better performance
Angles(9,1) = 40;

% add random jitter
TrueAngles = Angles + (rand(size(Angles))-0.5)*inpMaxError*2;

% add zero degree projection as 1st projection (do this only when
% numProjs is even number)
TrueAngles = cat(1,[0 0 0],TrueAngles); numProjs = numProjs+1;
Angles = cat(1,[0 0 0],Angles); 

% indices which will be used for first three projections for special Euler
% angle calculation using common line
FirstThreeAngleInds_ori = [1 7 9];
%FirstThreeAngleInds_ori = [1 2 3];

% re-order the angle array accordingly
NewAngleInd = FirstThreeAngleInds_ori;
for i=1:size(Angles,1)
  if sum(i==FirstThreeAngleInds_ori)==0
    NewAngleInd(end+1) = i;
  end
end
    
Angles = Angles(NewAngleInd,:);
TrueAngles = TrueAngles(NewAngleInd,:);




FirstThreeAngleInds = [1 2 3];

%%

save(angle_filename,'TrueAngles');
%pause
%% calculate projections
pj = zeros(64,64,numProjs);

parfor pj_num = 1:numProjs
    if mod(pj_num,1)==0
        fprintf('Calculating projection #%d/%d\n',pj_num,numProjs)
    end
%    tmp = calculate3Dprojection_interp(modelK,exp_angles_jittered(pj_num,1),exp_angles_jittered(pj_num,2),exp_angles_jittered(pj_num,3)); 
%    pj(:,:,pj_num) = tmp(nc_padded-n2:nc_padded+n2-1,nc_padded-n2:nc_padded+n2-1);
    tmp = calculate3Dprojection(Model,TrueAngles(pj_num,1),TrueAngles(pj_num,2),TrueAngles(pj_num,3)); 
   pj(:,:,pj_num) = tmp;%(nc_padded-n2:nc_padded+n2-1,nc_padded-n2:nc_padded+n2-1);
end
save(pj_filename,'pj');


%% get best in-plane angles for given first three projections

% search angle range and step
Range = 4;
Step = 0.05;


Pj1 = pj(:,:,FirstThreeAngleInds(1)); Angle1 = Angles(FirstThreeAngleInds(1),:);
Pj2 = pj(:,:,FirstThreeAngleInds(2)); Angle2 = Angles(FirstThreeAngleInds(2),:);
Pj3 = pj(:,:,FirstThreeAngleInds(3)); Angle3 = Angles(FirstThreeAngleInds(3),:);

[minAngle12_1, minAngle12_2] = get_CMline_bestInpAngles_LS_parallel(Pj1,Pj2,Angle1,Angle2,Range,Step);
[minAngle13_1, minAngle13_3] = get_CMline_bestInpAngles_LS_parallel(Pj1,Pj3,Angle1,Angle3,Range,Step);
[minAngle23_2, minAngle23_3] = get_CMline_bestInpAngles_LS_parallel(Pj2,Pj3,Angle2,Angle3,Range,Step);

%%
%opts = optimset('lsqcurvefit');
opts = optimset('Display','on');
% fitting for x axis projection
init_guess = [Angle2 Angle3];
testdata = [minAngle12_1 minAngle12_2 minAngle13_1 minAngle13_3 minAngle23_2 minAngle23_3];
%fun = @(p,xdata) CMline_calc_threelines_error(Angle1,p(1:3),p(4:6),xdata(1:2),xdata(3:4),xdata(5:6));
fun = @(p,xdata) CMline_calc_threelines_error(Angle1,p(1:3),p(4:6))*180/pi;
[p,fminres,~,eflag] = lsqcurvefit(fun,init_guess,testdata,testdata,[],[],opts);


currRefAngles = zeros(size(Angles));
currRefAngles(1,:) = Angle1;
currRefAngles(2,:) = p(1:3);
currRefAngles(3,:) = p(4:6);

%%

% search angle range and step
Range = 4;
Step = 0.5;

% initilize arrays
ImgSize = size(pj,1);
[initX, initY] = ndgrid(0,-ceil((ImgSize-1)/2):floor((ImgSize-1)/2));

AngleSearchRange = -Range:Step:Range;
RangeLength = length(AngleSearchRange);

tic

% loop over all projections

TotErrAr = zeros(1,size(Angles,1));
for i=4:size(Angles,1)  
  fprintf(1,'Running i=%d j=%d',i);toc
  
  Errarr = zeros(RangeLength,RangeLength,RangeLength,length(1:(i-1)));
  parfor PhiInd = 1:RangeLength
    Phi = AngleSearchRange(PhiInd) + Angles(i,1)
    for TheInd = 1:RangeLength
      The = AngleSearchRange(TheInd) + Angles(i,2)
      for PsiInd = 1:RangeLength
        for j=1:(i-1)          
          Psi = AngleSearchRange(PsiInd) + Angles(i,3)
          
          Matrix1 = get_Euler_matrix(Angles(j,:));
          Matrix2 = get_Euler_matrix([Phi The Psi]);

          % get in-plane vectors from the rotation matrices
          [invnormvec1, invnormvec2] = get_commonline_vectors(Matrix1, Matrix2, [0;0;1]);

          % convert vectors to angles
          ComAngle1 = angle(invnormvec1(1)+invnormvec1(2)*1i) - pi/2;
          ComAngle2 = angle(invnormvec2(1)+invnormvec2(2)*1i) - pi/2;


          % get common lines, this implementation is using DFT, but can be
          % done with FFT interpolation or real space rotation
          CMline1F = get_Fouriercommonline_from_angle(pj(:,:,j),ComAngle1, initX, initY);
          CMline2F = get_Fouriercommonline_from_angle(pj(:,:,i),ComAngle2, initX, initY);
          
          % get common line error (can be different metric)
          Errarr(PhiInd,TheInd,PsiInd,j) = sum(abs(CMline1F-CMline2F).^2);
        end
      end
    end
  end
  
  % get Euler angles for minimum error
  ErrarrSum = squeeze(sum(Errarr,4));
  [minErr,minInd] = min(ErrarrSum(:));
  TotErrAr(i) = minErr;
  [minI,minJ,minK] = ind2sub(size(ErrarrSum),minInd);
  currRefAngles(i,:) = [AngleSearchRange(minI) + Angles(i,1), AngleSearchRange(minJ) + Angles(i,2), AngleSearchRange(minK) + Angles(i,3)];

  save(sprintf('./results/R_%s_currRefAngles.mat',savePreFix),'currRefAngles');
  save(sprintf('./results/R_%s_TotErrAr.mat',savePreFix),'TotErrAr');
end
    

    

      
  
  
  
  
  
  
  
  
  
  
  
  
  
  

