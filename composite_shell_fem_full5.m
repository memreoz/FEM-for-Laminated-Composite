function results = composite_shell_fem_full5()
% ============================================================
% Finite Element Analysis of a Variable-Stiffness Composite Plate
% ============================================================
% - 5x5 shell elements (4-node Mindlin-Reissner)
% - Variable fiber angle θ(x,y)
% - User-defined loads and boundary conditions per edge
% - Outputs: element centers, displacements, principal stresses, directions
% ============================================================

% === USER-DEFINED BOUNDARY CONDITIONS AND LOADS ===
bc.fix_bottom = true;
bc.fix_top = false;
bc.fix_left = false;
bc.fix_right = false;

bc.simply_bottom = false;
bc.simply_top = true;
bc.simply_left = true;
bc.simply_right = false;

load_edges.force_bottom = [1e5, 0, 0]; % [Fx, Fy, Fz]
load_edges.force_top = [0, 0, 0];
load_edges.force_left = [0, 0, 0];
load_edges.force_right = [0, 0, 0];

load_edges.moment_bottom = [0, 0, 0]; % [Mx, My, Mz]
load_edges.moment_top = [1000, 0, 0];
load_edges.moment_left = [0, 0, 0];
load_edges.moment_right = [0, 0, 0];

% === GEOMETRY AND MESH ===
nx = 5; ny = 5; Lx = 1; Ly = 1; t = 0.002;
dx = Lx / nx; dy = Ly / ny;
numNodes = (nx+1)*(ny+1); ndof = numNodes * 5;
nodes = zeros(numNodes,2);
count = 1;
for j = 0:ny
    for i = 0:nx
        nodes(count,:) = [i*dx, j*dy];
        count = count + 1;
    end
end
nodeID = reshape(1:numNodes, nx+1, ny+1)';
elements = zeros(nx*ny, 4);
e = 1;
for j = 1:ny
    for i = 1:nx
        n1 = nodeID(j,i); n2 = nodeID(j,i+1);
        n3 = nodeID(j+1,i+1); n4 = nodeID(j+1,i);
        elements(e,:) = [n1 n2 n3 n4];
        e = e + 1;
    end
end
nelem = size(elements,1);

% === MATERIAL ===
E1 = 135e9; E2 = 10e9; G12 = 5e9; nu12 = 0.3;
Q = get_lamina_stiffness(E1,E2,G12,nu12);

% === FIBER ANGLE FIELD ===
theta_func = @(x,y) -45 + 90 * x;
fiber_angle = zeros(nelem,1);
for e = 1:nelem
    xy = nodes(elements(e,:),:);
    xc = mean(xy(:,1)); yc = mean(xy(:,2));
    fiber_angle(e) = theta_func(xc, yc);
end

% === INITIALIZE LOAD VECTORS ===
element_forces = zeros(nelem, 6);
tol = 1e-6;
bottom_nodes = find(abs(nodes(:,2)) < tol);
top_nodes = find(abs(nodes(:,2) - Ly) < tol);
left_nodes = find(abs(nodes(:,1)) < tol);
right_nodes = find(abs(nodes(:,1) - Lx) < tol);

top_elements = find(any(ismember(elements, top_nodes),2));
bottom_elements = find(any(ismember(elements, bottom_nodes),2));
left_elements = find(any(ismember(elements, left_nodes),2));
right_elements = find(any(ismember(elements, right_nodes),2));

for e = bottom_elements'
    element_forces(e,1:3) = element_forces(e,1:3) + load_edges.force_bottom;
    element_forces(e,4:6) = element_forces(e,4:6) + load_edges.moment_bottom;
end
for e = top_elements'
    element_forces(e,1:3) = element_forces(e,1:3) + load_edges.force_top;
    element_forces(e,4:6) = element_forces(e,4:6) + load_edges.moment_top;
end
for e = left_elements'
    element_forces(e,1:3) = element_forces(e,1:3) + load_edges.force_left;
    element_forces(e,4:6) = element_forces(e,4:6) + load_edges.moment_left;
end
for e = right_elements'
    element_forces(e,1:3) = element_forces(e,1:3) + load_edges.force_right;
    element_forces(e,4:6) = element_forces(e,4:6) + load_edges.moment_right;
end

% === BOUNDARY CONDITIONS ===
fixed_dofs = [];
fixed_dofs = [fixed_dofs, get_dofs(bottom_nodes, bc.fix_bottom, bc.simply_bottom)];
fixed_dofs = [fixed_dofs, get_dofs(top_nodes, bc.fix_top, bc.simply_top)];
fixed_dofs = [fixed_dofs, get_dofs(left_nodes, bc.fix_left, bc.simply_left)];
fixed_dofs = [fixed_dofs, get_dofs(right_nodes, bc.fix_right, bc.simply_right)];
free_dofs = setdiff(1:ndof, unique(fixed_dofs));

% === ASSEMBLY ===
K = zeros(ndof); F = zeros(ndof,1);
for e = 1:nelem
    nodes_e = elements(e,:);
    xy = nodes(nodes_e,:);
    theta = fiber_angle(e);
    D = transform_D(Q, theta);
    Ke = shell4node_FE(xy, t, D);
    Fe = element_force_vector(element_forces(e,:), t);
    dofs = zeros(1,20);
    for i = 1:4
        dofs((i-1)*5+1:i*5) = (nodes_e(i)-1)*5 + (1:5);
    end
    K(dofs,dofs) = K(dofs,dofs) + Ke;
    F(dofs) = F(dofs) + Fe;
end

% === SOLVE ===
U = zeros(ndof,1);
U(free_dofs) = K(free_dofs,free_dofs) \ F(free_dofs);
u = U(1:5:end); v = U(2:5:end); w = U(3:5:end);

% === STRESS CALCULATION: PRINCIPAL STRESSES ===
sig1 = zeros(nelem,1); sig2 = zeros(nelem,1); theta_p = zeros(nelem,1);
elem_center = zeros(nelem,2);
for e = 1:nelem
    nodes_e = elements(e,:); xy = nodes(nodes_e,:);
    theta = fiber_angle(e); D = transform_D(Q, theta);
    dofs = zeros(1,20);
    for i = 1:4
        dofs((i-1)*5+1:i*5) = (nodes_e(i)-1)*5 + (1:5);
    end
    Ue = U(dofs);
    [~, dN_dxi] = shapeQuad4(0,0);
    J = dN_dxi * xy; dN_dx = J \ dN_dxi;
    Bm = zeros(3,20);
    for a = 1:4
        Bm(:,(a-1)*5+1:(a-1)*5+5) = [dN_dx(1,a), 0, 0, 0, 0;
                                     0, dN_dx(2,a), 0, 0, 0;
                                     dN_dx(2,a), dN_dx(1,a), 0, 0, 0];
    end
    strain = Bm * Ue;
    stress = D * strain;
    sx = stress(1); sy = stress(2); txy = stress(3);
    R = sqrt((sx - sy)^2 / 4 + txy^2);
    sig1(e) = (sx + sy)/2 + R;
    sig2(e) = (sx + sy)/2 - R;
    theta_p(e) = 0.5 * atan2(2*txy, sx - sy);
    elem_center(e,:) = mean(xy);
end

% === PLOTS ===
figure;
subplot(1,3,1);
trisurf([elements(:,[1 2 3]); elements(:,[1 3 4])], nodes(:,1), nodes(:,2), u, 'FaceColor','interp');
title('X Displacement u [m]'); colorbar; view(45,30);
xlabel('x'); ylabel('y'); zlabel('u');

subplot(1,3,2);
scatter(elem_center(:,1), elem_center(:,2), 80, fiber_angle, 'filled');
title('Fiber Orientation \theta(x,y) [deg]', 'Interpreter','latex');
xlabel('x'); ylabel('y'); colorbar; axis equal;

subplot(1,3,3);
scatter(elem_center(:,1), elem_center(:,2), 80, sig1/1e6, 'filled');
title('Principal Stress \sigma_1 [MPa]', 'Interpreter','latex');
xlabel('x'); ylabel('y'); colorbar; axis equal;

figure;
quiver(elem_center(:,1), elem_center(:,2), cos(theta_p), sin(theta_p), 0.2);
title('Principal Stress Directions');
xlabel('x'); ylabel('y'); axis equal;

% === OUTPUT RESULTS STRUCT ===
results.nodes = nodes;
results.elements = elements;
results.displacement_u = u;
results.displacement_v = v;
results.displacement_w = w;
results.fiber_angle = fiber_angle;
results.element_centers = elem_center;
results.principal_stress_1 = sig1;
results.principal_stress_2 = sig2;
results.principal_direction = theta_p;
end

function dofs = get_dofs(node_list, is_fixed, is_simply)
dofs = [];
for i = 1:length(node_list)
    base = (node_list(i)-1)*5;
    if is_fixed
        dofs = [dofs, base + (1:5)];
    elseif is_simply
        dofs = [dofs, base + [2 3 5]];
    end
end
end
%% === SUPPORTING FUNCTIONS ===
function Q = get_lamina_stiffness(E1,E2,G12,nu12)
Q11 = E1 / (1 - nu12^2 * E2/E1);
Q22 = E2 / (1 - nu12^2 * E2/E1);
Q12 = nu12 * E2 / (1 - nu12^2 * E2/E1);
Q66 = G12;
Q = [Q11 Q12 0; Q12 Q22 0; 0 0 Q66];
end

function Dbar = transform_D(Q, theta)
m = cosd(theta); n = sind(theta);
T1 = [m^2, n^2, 2*m*n;
      n^2, m^2, -2*m*n;
     -m*n, m*n, m^2 - n^2];
T2 = [m^2, n^2, m*n;
      n^2, m^2, -m*n;
     -2*m*n, 2*m*n, m^2 - n^2];
Dbar = T1' * Q * T2;
end

function Ke = shell4node_FE(xy, t, D)
Ke = zeros(20);
D_m = D; D_b = D * t^2 / 12;
D_s = mean(diag(D)) * t/5 * eye(2);
g = [-1 1]/sqrt(3);
[gx, gy] = meshgrid(g, g);
gauss_pts = [gx(:), gy(:)];
for i = 1:4
    xi = gauss_pts(i,1); eta = gauss_pts(i,2);
    [~, dN_dxi] = shapeQuad4(xi, eta);
    J = dN_dxi * xy; detJ = det(J); dN_dx = J \ dN_dxi;
    Bm = zeros(3,20); Bb = zeros(3,20); Bs = zeros(2,20);
    for a = 1:4
        Bm(:,(a-1)*5+1:(a-1)*5+5) = [dN_dx(1,a), 0, 0, 0, 0;
                                     0, dN_dx(2,a), 0, 0, 0;
                                     dN_dx(2,a), dN_dx(1,a), 0, 0, 0];
        Bb(:,(a-1)*5+1:(a-1)*5+5) = [0, 0, 0, dN_dx(1,a), 0;
                                     0, 0, 0, 0, dN_dx(2,a);
                                     0, 0, 0, dN_dx(2,a), dN_dx(1,a)];
        Bs(:,(a-1)*5+1:(a-1)*5+5) = [0, 0, dN_dx(1,a), -1, 0;
                                     0, 0, dN_dx(2,a), 0, -1];
    end
    Ke = Ke + (Bm'*D_m*Bm + Bb'*D_b*Bb + Bs'*D_s*Bs) * detJ * t;
end
end

function Fe = element_force_vector(ef, t)
Fe = zeros(20,1);
for i = 0:3
    Fe(i*5+1) = ef(1)/4;
    Fe(i*5+2) = ef(2)/4;
    Fe(i*5+3) = ef(3)/4;
    Fe(i*5+4) = ef(4)/4;
    Fe(i*5+5) = ef(5)/4;
end
end

function [N, dN_dxi] = shapeQuad4(xi, eta)
N = 0.25 * [(1-xi)*(1-eta);
            (1+xi)*(1-eta);
            (1+xi)*(1+eta);
            (1-xi)*(1+eta)];
dN_dxi = 0.25 * [ -(1-eta),  (1-eta),  (1+eta), -(1+eta);
                 -(1-xi),  -(1+xi),   (1+xi),   (1-xi) ];
end
