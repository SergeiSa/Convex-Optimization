clear; close all;

M = [1 0 1; 0 5 0; 1 0 3];
C = [1 7 2];
y = 1;

x = quadprog(M,[],[],[],C,y)

disp(['Residual of the quadprog solution: ', num2str(norm(C*x - y)), '; norm of the solution: ', num2str(norm(x))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimoptions('quadprog', 'Display', 'none');

Count = 1000;
tic;
for i = 1:Count
x = quadprog(2*M,[],[],[],C,y, [],[],[],options);
end

timer = toc;

average_time = timer / Count;
fprintf('Average time: %.2e ms \n', average_time*10^3);