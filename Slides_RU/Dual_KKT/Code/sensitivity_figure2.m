clear; close all; clc;

%% Quadratic (choose center to the LEFT of feasible set)
Q = [2 0.6; 0.6 1.5];
xc = [-2.0; 0.5];              % center of ellipses (unconstrained minimizer)

% Build q so that xc is the minimizer
q = -2*Q*xc;

f = @(x) x'*Q*x + q'*x;

%% Polytope (smooth pentagon, no right angles)
poly_pts = [
    -0.5  1.8;
     1.2  2.2;
     2.2  0.8;
     1.6 -1.5;
    -0.8 -1.2
];

perturbation_index = 1;

% Convert to Ax <= b form
K = convhull(poly_pts(:,1), poly_pts(:,2));
poly_pts = poly_pts(K(1:end-1),:); % ordered vertices

n = size(poly_pts,1);
A = zeros(n,2); b = zeros(n,1);

for i = 1:n
    p1 = poly_pts(i,:);
    p2 = poly_pts(mod(i,n)+1,:);
    edge = p2 - p1;
    normal = [edge(2), -edge(1)]; % outward normal
    normal = normal / norm(normal);
    
    A(i,:) = normal;
    b(i) = normal * p1';
end

%% Choose SW->NE constraint to perturb (positive slope line)
% We enforce it manually for clarity
% a = [1; -1]; 
% a = a / norm(a);
% 
% b_active = 0.8;
% 
% % Replace one constraint with this one
% A(end,:) = a';
% b(end) = b_active;

u = 0.6;

%% Solve QPs
H = 2*Q;
f_qp = q;

x_star_0 = quadprog(H, f_qp, A, b);

b_pert = b;
b_pert(perturbation_index) = b(perturbation_index) + u;

x_star_u = quadprog(H, f_qp, A, b_pert);

v0 = f(x_star_0);
v1 = f(x_star_u);

% additional levels
v2 = v0 + 3.0;   % infeasible contour
v3 = v1 - 0.8;   % outer contour

%% Function to draw ellipse level set
draw_ellipse = @(v, style) ...
    drawEllipse(Q, xc, v - f(xc), style);

%% Plot
figure('Color', 'w'); hold on; axis equal;

% Polytope
fill(poly_pts(:,1), poly_pts(:,2), [0.9 0.9 0.9], ...
    'EdgeColor','k','LineWidth',1.5,'FaceAlpha',0.4);

% Constraint lines (SW -> NE)
xx = linspace(-2,2,200);
a_active   = A(perturbation_index, :);
b_active   = b(perturbation_index);
b_active_p = b_pert(perturbation_index);

yy0 = (b_active   - a_active(1)*xx)/a_active(2);
yyu = (b_active_p - a_active(1)*xx)/a_active(2);

plot(xx, yy0, 'r-', 'LineWidth',1);
plot(xx, yyu, 'r--', 'LineWidth',2);

% Elliptical contours (analytic)
draw_ellipse(v2, 'k:');   % infeasible (inner)
draw_ellipse(v0, 'b:');   % tangent at x*(0)
draw_ellipse(v1, 'g:');   % tangent at x*(u)
draw_ellipse(v3, 'k:');   % outer

% Optimal points
plot(x_star_0(1), x_star_0(2), 'ko', 'MarkerFaceColor','k', 'MarkerSize',8);
plot(x_star_u(1), x_star_u(2), 'bo', 'MarkerFaceColor','b', 'MarkerSize',8);

% Arrow (movement)
quiver(x_star_0(1), x_star_0(2), ...
       x_star_u(1)-x_star_0(1), ...
       x_star_u(2)-x_star_0(2), ...
       0, 'k', 'LineWidth',2, 'MaxHeadSize',0.6);

% Labels
text(x_star_0(1)+0.1, x_star_0(2), 'x^*(0)');
text(x_star_u(1)+0.1, x_star_u(2), 'x^*(u)');

xlabel('x_1'); ylabel('x_2');
title('QP Sensitivity: geometric illustration');

x_min = -0.5 * (Q \ q);
plot(x_min(1), x_min(2), 'md', 'MarkerFaceColor','m', 'MarkerSize',8);
text(x_min(1)+0.1, x_min(2), 'x_{unconstr}');

grid off;

xlim([-3.5, 2.5])
ylim([-2, 3])


%% ===== Helper function =====
function drawEllipse(Q, xc, level, style)
    if level <= 0, return; end
    
    [V,D] = eig(Q);
    t = linspace(0,2*pi,200);
    circle = [cos(t); sin(t)];
    
    % ellipse transform
    A = V * diag(1./sqrt(diag(D)));
    pts = xc + sqrt(level) * (A * circle);
    
    plot(pts(1,:), pts(2,:), style, 'LineWidth',1.5);
end