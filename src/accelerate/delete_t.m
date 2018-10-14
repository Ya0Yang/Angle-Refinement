
A = rand(3,3);
B = rand(3,3);
N = 10000;
tic
for i = 1:N
   [invnormvec1, invnormvec2] = get_commonline_vectors(A, B, [0;0;1]);
end
toc

tic
for i = 1:N
    [invnormvec11, invnormvec12] = accel_get_commonline_vectors(A, B);

end
toc


