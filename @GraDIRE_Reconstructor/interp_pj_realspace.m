function obj = interp_pj_realspace(obj)
tic
projections  = obj.InputProjections;
Num_pj       = obj.NumProjs;
dimx         = obj.Dim1;
dimy         = obj.Dim2;

phiangles    = obj.InputAngles(:,1);
thetaangles  = obj.InputAngles(:,2);
psiangles    = obj.InputAngles(:,3);

vec1         = obj.vector1;
vec2         = obj.vector2;
vec3         = obj.vector3;

dtype        = obj.dtype;
sigma        = obj.sigma;
n1_oversampled = obj.n1_oversampled;
n2_oversampled = obj.n2_oversampled;

if sigma
    [Y_big,X_big,Z_big] = meshgrid(1:n2_oversampled,1:n1_oversampled,1:n2_oversampled,dtype);
    x_cen = floor(obj.n2_oversampled,dtype/2);
    y_cen = floor(obj.n1_oversampled,dtype/2);
    z_cen = floor(obj.n2_oversampled,dtype/2);
    kernel = (X_big-x_cen).^2 + (Y_big-y_cen).^2 + (Z_big-z_cen).^2;
    kernel = exp(-kernel/sigma^2);
    obj.kernel = kernel;
end
clear X_big Y_big Z_big

ncy = round((dimy+1)/2); 
ncx = round((dimx+1)/2); 
k1 = cast((-1*ceil((dimx-1)/2):1:floor((dimx-1)/2)),dtype );
k2 = cast((-1*ceil((dimy-1)/2):1:floor((dimy-1)/2)),dtype );
k3 = k1;
[YY,XX,ZZ] = meshgrid(k1,k2,k3);
XX = XX(:)';
YY = YY(:)';
ZZ = ZZ(:)';
rot_pjs = zeros(dimy,dimx,dimx,Num_pj,dtype);

ncy_big = round((n2_oversampled+1)/2); 
ncx_big = round((n1_oversampled+1)/2); 
[Y, X, Z] = meshgrid((1:n2_oversampled) - ncy_big, (1:n1_oversampled) - ncx_big, 0);
Y = Y(:)';
X = X(:)';
Z = Z(:)';
xj = zeros(n2_oversampled,n1_oversampled,Num_pj);
yj = zeros(n2_oversampled,n1_oversampled,Num_pj);
zj = zeros(n2_oversampled,n1_oversampled,Num_pj);

Rot_x = zeros(dimy, dimx, dimy, Num_pj, dtype);
Rot_y = zeros(dimy, dimx, dimy, Num_pj, dtype);

if obj.use_parallel
  parforArg = Inf;
else
  parforArg = 0;
end

parfor (k = 1:Num_pj, parforArg)
    phi   = phiangles(k);
    theta = thetaangles(k);
    psi   = psiangles(k);
    pj    = projections(:,:,k);    
    
    R1 = MatrixQuaternionRot(vec1,phi);
    R2 = MatrixQuaternionRot(vec2,theta);
    R3 = MatrixQuaternionRot(vec3,psi);
    R =(R1*R2*R3)';
    
    rotCoords = R(1:2,:)*[XX; YY; ZZ];
    rot_x  = double(rotCoords(1,:));
    rot_y  = double(rotCoords(2,:));
    rot_x = reshape(rot_x, [dimy dimx dimy])+ncx;
    rot_y = reshape(rot_y, [dimy dimx dimy])+ncy;
    Rot_x(:,:,:,k) = rot_x;
    Rot_y(:,:,:,k) = rot_y;
    rot_pj = splinterp2(pj, rot_y, rot_x);    
    rot_pjs(:,:,:,k) = rot_pj;        
    
    rotCoords = R'*[X; Y; Z];
    xj(:,:,k) = reshape((rotCoords(1,:)),[n2_oversampled,n1_oversampled])+ncy_big;
    yj(:,:,k) = reshape((rotCoords(2,:)),[n2_oversampled,n1_oversampled])+ncy_big;
    zj(:,:,k) = reshape((rotCoords(3,:)),[n2_oversampled,n1_oversampled])+ncy_big;
end
sum_rot_pjs = sum(rot_pjs,4);

obj.sum_rot_pjs = sum_rot_pjs;
obj.xj = xj;
obj.yj = yj;
obj.zj = zj;
obj.Rot_x = Rot_x;
obj.Rot_y = Rot_y;

timeTakenToInterp = toc;
timeTakenToInterp = round(10*timeTakenToInterp)./10;
fprintf('GraDIRE: projections interpolated in %.12g seconds.\n\n',timeTakenToInterp);

end