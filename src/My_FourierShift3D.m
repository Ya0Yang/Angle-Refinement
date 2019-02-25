function shiftModel = My_FourierShift3D(model,dy,dx,dz)

%dy is first dimension, dx second, dz third

ny = size(model,1); nx = size(model,2);nz = size(model,3);
[X Y Z] = meshgrid(-ceil((nx-1)/2):floor((nx-1)/2),-ceil((ny-1)/2):floor((ny-1)/2),-ceil((nz-1)/2):floor((nz-1)/2));
F = My_IFFTN(model);
Pfactor = exp(2*pi*1i*(dx*X/nx + dy*Y/ny + dz*Z/nz));
shiftModel = My_FFTN(F.*Pfactor);

end