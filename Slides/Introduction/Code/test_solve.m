clear; close all;


M = [1 0 1; 0 5 0; 1 0 3];
C = [1 7 2];
y = 1;

T1 = M \ C';
t2 = (C*T1) \ y;
x = T1*t2;


disp(['Residual of the backslash solution: ', num2str(norm(C*x - y)), '; norm of the solution: ', num2str(norm(x))])

Count = 10000;
tic;
for i = 1:Count
T1 = M \ C';
t2 = (C*T1) \ y;
x = T1*t2;
end
timer = toc;

average_time = timer / Count;
fprintf('Average time: %.2e ms \n', average_time*10^3);
