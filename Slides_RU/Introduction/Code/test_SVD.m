clear; close all;


M = [1 0 1; 0 5 0; 1 0 3];
C = [1 7 2];
y = 1;

tol = 10^(-5);

[P, N] = pinv_null(C, tol);

x = ( eye(3) - N*((N'*M*N) \ (N'*M)) ) * P*y


disp(['Residual of the SVD solution: ', num2str(norm(C*x - y)), '; norm of the solution: ', num2str(norm(x))])

Count = 10000;
tic;
for i = 1:Count
[P, N] = pinv_null(C, tol);
x = ( eye(3) - N*((N'*M*N) \ (N'*M)) ) * P*y;
end
timer = toc;

average_time = timer / Count;
fprintf('Average time: %.2e ms \n', average_time*10^3);
