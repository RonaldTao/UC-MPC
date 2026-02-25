%% Tube MPC trajectory (nominal plan + ancillary feedback in invariant tube)
% This script:
%   1) Discretizes the open-loop model x^+ = Ad x + Bd u + disturbance
%   2) Designs LQR gain K for tube feedback
%   3) Builds disturbance set W and disturbance invariant set Z (approx. via finite sum)
%   4) Tightens constraints: X_bar = X \ominus Z, U_bar = U \ominus KZ
%   5) Solves nominal MPC on (Ad,Bd) with tightened constraints
%   6) Simulates actual system with u = u_nom + K(x_real - x_nom) and disturbances

%% 0) Environment setup
% add path below for YALMIP and solver
addpath(genpath('C:\Users\xxx\matlab\YALMIP'))
addpath(genpath('C:\Program Files\Mosek\10.1\toolbox\r2017a'))
yalmip('clear')
clear all
clc

figure; hold on

%% 1) System definition and discretization
Ac = [0 0.0067 1.34;
      0 -0.869 43.2;
      0 0.993 -1.34];

Bc = [0.169  0.252;
      -17.3  -1.58;
      -0.169 -0.252];

B_perp = null(Bc');
Bu = B_perp*0.15;

Cc = [1 0 1;
      1 0 0];
Dc = [0 0;
      0 0];

nx = 3;
nu = 2;
ny = 2;

dT = 4e-3;

% Discretize with inputs [Bc Bu]
[Ad, Bd_total, Cd, Dd] = ssdata(c2d(ss(Ac, [Bc Bu], Cc, [Dc [0;0]]), dT));
Bd   = Bd_total(:,1:2);
Bu_d = Bd_total(:,3);

%% 2) Ancillary feedback gain K (LQR)
Q_lqr = diag([1,10,10]);
R_lqr = 0.1;

[K,P,~] = dlqr(Ad, Bd, Q_lqr, R_lqr);
K = -K;                  % use u = Kx (consistent with your original)
Ak = Ad + Bd*K;           % closed-loop tube error dynamics

x0 = [0;0;0]; 

%% 3) Disturbance set W (vertex form -> inequality form)
% Disturbance model in your code:
%   w_k = Bd*f_uncert + Bu_d*w_uncert
% with bounds embedded by the v1,v2,v3 construction
v1 = abs(Bd(1,1)*2.4*1.8) + abs(Bd(1,2)*0.9*1.8) + abs(Bu_d(1)*1);
v2 = abs(Bd(2,1)*2.4*1.8) + abs(Bd(2,2)*0.9*1.8) + abs(Bu_d(2)*1);
v3 = abs(Bd(3,1)*2.4*1.8) + abs(Bd(3,2)*0.9*1.8) + abs(Bu_d(3)*1);

% 8 vertices of a box (then scaled by 1.8 like your original)
W_v = 1.8 * [ v1,  v2,  v3;
             -v1,  v2,  v3;
             -v1, -v2,  v3;
             -v1, -v2, -v3;
             -v1,  v2, -v3;
              v1, -v2,  v3;
              v1, -v2, -v3;
              v1,  v2, -v3];

% Convert vertex to inequality (H z <= h), store as [H h]
[W, w, ~, ~] = vert2lcon(W_v);
W = [W w];

%% 4) Disturbance invariant set Z (finite approximation)
% Z approx = sum_{n=0}^{n_max} Ak^n W
Z_v  = [0,0,0];
n_max = 5;

for n = 0:n_max
    n   % keep your progress print
    FW  = (Ak^n * W_v')';    % push disturbance vertices through Ak^n
    Z_v = Minkowski_sum(Z_v, FW);
end

[Z, z, ~, ~] = vert2lcon(Z_v);
Z = [Z z];

%% 5) Original state/input constraints (X, U) as polytopes
% State constraint vertices (box)
Cx_v = [ 10,  100,   4;
        -10,  100,   4;
        -10, -100,   4;
        -10,  100,  -4;
        -10, -100,  -4;
         10, -100,   4;
         10,  100,  -4;
         10, -100,  -4];

[Cx, cx, ~, ~] = vert2lcon(Cx_v);
Cx = [Cx cx];

% Plot original X in (x1,x2) projection (as in your code)
[k_hull,~] = convhull(Cx_v);
p = patch(Cx_v(k_hull,1)', Cx_v(k_hull,2)', 'm');
p.FaceAlpha = 0.15;

% Input constraint vertices (box)
Cu_v = [ 25,  22;
         25, -22;
        -25,  22;
        -25, -22];

[Cu, cu, ~, ~] = vert2lcon(Cu_v);
Cu = [Cu cu];

%% 6) Tightened constraints: X_bar = X \ominus Z, U_bar = U \ominus KZ
% Tightened state constraint
Cx_bar = P_difference(Cx, Z);

% Convert to vertices for plotting
[Cx_bar_v, ~, ~] = lcon2vert(Cx_bar(:,1:end-1), Cx_bar(:,end), [], []);
[k_hull2,~] = convhull(Cx_bar_v);
p2 = patch(Cx_bar_v(k_hull2,1)', Cx_bar_v(k_hull2,2)', 'r');
p2.FaceAlpha = 0.35;

% Tightened input constraint via KZ
KZ_v = (K * Z_v')';
[KZ, kz, ~, ~] = vert2lcon(KZ_v);
KZ = [KZ kz];

Cu_bar = P_difference(Cu, KZ);

%% 7) Nominal MPC on tightened constraints (tube MPC nominal plan)
Q = eye(nx); 
R = eye(nu); 

N = 0.2/dT;
N = round(N);  % make sure integer

u = sdpvar(repmat(nu,1,N),   repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1), repmat(1,1,N+1));
r = sdpvar(repmat(ny,1,N+1), repmat(1,1,N+1));

constraints = [];
objective   = 0;

for k = 1:N
    % Stage cost (same as your original)
    if k == 1
        objective = objective + 100*norm(Cd*x{k} - r{k}, 2) + 1*norm(u{k}, 2);
    else
        objective = objective + 100*norm(Cd*x{k} - r{k}, 2) + 1*norm(u{k}, 2) ...
                              + 100*norm(u{k} - u{k-1}, 2);
    end

    % Nominal dynamics
    constraints = [constraints, x{k+1} == Ad*x{k} + Bd*u{k}];

    % Tightened bounds (exactly your scalar form)
    constraints = [constraints, ...
        -Cu_bar(1,3) <= u{k}(1) <= Cu_bar(1,3), ...
        -Cu_bar(2,3) <= u{k}(2) <= Cu_bar(2,3), ...
        -Cx_bar(2,4) <= x{k+1}(1) <= Cx_bar(2,4), ...
        -Cx_bar(1,4) <= x{k+1}(2) <= Cx_bar(1,4), ...
        -Cx_bar(3,4) <= x{k+1}(3) <= Cx_bar(3,4)];
end

objective = objective + 100*norm(Cd*x{N+1} - r{N+1}, 2);

parameters_in = {x{1}, [r{:}]};
solutions_out = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective, sdpsettings('solver','mosek'), ...
                       parameters_in, solutions_out);

%% 8) Simulation setup (time + reference)
T      = 15;
steps  = floor(T/dT);
steps2 = steps + 500;

Time = zeros(1, steps+1);
for k = 1:steps+1
    Time(k) = (k-1)*dT;
end

ref = zeros(ny, steps2);
for k = 1:steps2
    if (k-1 < floor(steps/2))
        ref(:,k) = [9; 6.5];
    else
        ref(:,k) = [0; 0];
    end
end

%% 9) Nominal closed-loop simulation (tube MPC nominal trajectory)
xk = [0;0;0];

Xhist  = zeros(nx, steps+1);
Uhist  = zeros(nu, steps+1);

% Store nominal state trajectory X(:,i) (avoid X=[X,x])
X = zeros(nx, steps+1);

count = 0;

for i = 1:steps+1
    future_r   = ref(:, i:i+N);
    Xhist(:,i) = xk;
    X(:,i)     = xk;

    input = {xk, future_r};
    [solutions, diagnostics] = controller{input};

    Uvec = solutions{1};
    Uvec = Uvec(:);               % force column
    uk   = [Uvec(1); Uvec(2)];    % first input
    uk   = uk(:);

    % Propagate nominal system (no disturbance here)
    xk = Ad*xk + Bd*uk;

    Uhist(:,i) = uk;

    count = count + 1;
    if mod(count,100) == 0
        count
    end
end

%% 10) Real system simulation with tube feedback and disturbances
X_real = zeros(nx, steps+1);
U      = zeros(nu, steps+1);
fhist  = zeros(nu, steps+1);

x_real = [0;0;0];

for i = 1:steps+1
    % Log current real state
    X_real(:,i) = x_real;

    % Tube MPC control law: u = u_nom + K (x_real - x_nom)
    u_real = Uhist(:,i) + K*(x_real - X(:,i));
    U(:,i) = u_real;

    % Disturbances (same as your original)
    f_uncert = ([-0.8*sin(0.4*pi*dT*(i-1)) - 0.1*x_real(3)^2;
                  0.1 - 0.2*x_real(3)] * 1.8);
    w_uncert = 1*(sin(0.6*pi*dT*(i-1)));

    fhist(:,i) = f_uncert;

    % Real dynamics
    x_real = Ad*x_real + Bd*u_real + Bd*f_uncert + Bu_d*w_uncert;
end

%% 11) Plots
close all;

figure;
plot(Time, X_real(1,:) + X_real(3,:)); hold on
plot(Time, X(1,:) + X(3,:)); hold on
plot(Time, ref(1,1:steps+1));
title('theta, x1+x3');

figure;
plot(Time, X_real(1,:)); hold on
plot(Time, ref(2,1:steps+1));
title('gamma, x1');

figure;
plot(Time, X_real(1,:)); hold on
plot(Time, X_real(2,:)); hold on
plot(Time, X_real(3,:));
legend('x1','x2','x3')

figure;
plot(Time, X_real(3,:)); hold on
plot(Time, X(3,:));
yline(4,'-','State Constrain');
yline(-4,'-','State Constrain');
title('x3');

figure;
plot(Time, U(1,:));
yline(25,'-','Input Bound');
yline(-25,'-','Input Bound');
title('total de')

figure;
plot(Time, U(2,:));
yline(22,'-','Input Bound');
yline(-22,'-','Input Bound');
title('total df');

figure;
plot(Time, fhist(1,:));
title('f 1');

figure;
plot(Time, fhist(2,:));
title('f 2');

%% 12) Save results
save TMPC_umpc.mat Uhist
save TMPC_x.mat    X_real
save TMPC_u.mat    U

%% =======================================================================
% Helper functions (kept identical logic; only minor formatting)
% =======================================================================

function Cv = Minkowski_sum(Av, Bv) % Vertex representation
% C = A + B in vertex form
nA = size(Av,1);
nB = size(Bv,1);

Cv = zeros(nA*nB, size(Av,2));
id = 0;

for id_A = 1:nA
    for id_B = 1:nB
        id = id + 1;
        Cv(id,:) = Av(id_A,:) + Bv(id_B,:);
    end
end
end

function C = P_difference(A, B) % Inequality representation
% C = A \ominus B (Pontryagin difference)
% A = [U u] meaning U*x <= u
% Uses Kolmanovsky & Gilbert (1998) style support-function computation

U = A(:,1:end-1); u = A(:,end);
V = B(:,1:end-1); v = B(:,end);

options = optimset('Display','none');
for i = 1:length(u)
    % u_i <- u_i - max_{Vx<=v} U_i x
    u(i) = u(i) - U(i,:)*linprog(-U(i,:)', V, v, [], [], [], [], [], options);
end

C = [U u];
end