%% Vanilla MPC trajectory (u = u_opt) + online solve time measurement
% This script:
%   1) Discretizes the open-loop nominal model (Ac, Bc)
%   2) Builds a vanilla MPC (decision variable u{k})
%   3) Simulates the uncertain system dynamics (adds f_uncert and w_uncert)
%   4) Measures average solve time of the optimizer call
%   5) Saves trajectories and plots results

%% 0) Environment setup
% add path below for YALMIP and solver
addpath(genpath('C:\Users\xxx\matlab\YALMIP'))
addpath(genpath('C:\Program Files\Mosek\10.1\toolbox\r2017a'))
yalmip('clear')
clear all
clc
%% 1) Continuous-time model definition (NO state feedback for vanilla MPC)
Ac = [0 0.0067 1.34;
      0 -0.869 43.2;
      0 0.993 -1.34];

Bc = [0.169  0.252;
      -17.3  -1.58;
      -0.169 -0.252];

% Unmatched direction (used only via Bu_d*w_uncert in simulation)
B_perp = null(Bc');
Bu = B_perp*0.15;

% Output model (tracking)
Cc = [1 0 1;
      1 0 0];
Dc = [0 0;
      0 0];

nx = 3;    % number of states
nu = 2;    % number of inputs
ny = 2;    % number of outputs

%% 2) Discretization
dT = 4e-3;  % sampling time

% Discretize open-loop nominal system with inputs [Bc Bu]
[Ad, Bd_total, Cd, Dd] = ssdata(c2d(ss(Ac, [Bc Bu], Cc, [Dc [0;0]]), dT));
Bd   = Bd_total(:,1:2);
Bu_d = Bd_total(:,3);

%% 3) MPC setup
N = 0.2/dT;     % horizon length in steps
N = round(N);   % force integer for sdpvar indexing

% Decision variables over horizon:
%   u{k} : nu x 1
%   x{k} : nx x 1
%   r{k} : ny x 1 (reference)
%   p{k} : nx x 1 (your auxiliary/slack variable)
u = sdpvar(repmat(nu,1,N),   repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1), repmat(1,1,N+1));
r = sdpvar(repmat(ny,1,N+1), repmat(1,1,N+1));
p = sdpvar(repmat(nx,1,N),   repmat(1,1,N));

constraints = [];
objective   = 0;

for k = 1:N
    % Output tracking
    yk = Cd*x{k};

    % Stage cost (same as your original)
    if k == 1
        objective = objective + 100*norm(yk - r{k}, 2) + 1*norm(u{k}, 2) + 100*norm(p{k}, 2);
    else
        objective = objective + 100*norm(yk - r{k}, 2) + 1*norm(u{k}, 2) ...
                              + 100*norm(u{k} - u{k-1}, 2) + 100*norm(p{k}, 2);
    end

    % Nominal dynamics used inside MPC
    constraints = [constraints, x{k+1} == Ad*x{k} + Bd*(u{k})];

    % Input bounds (hard bounds)
    constraints = [constraints, -25 <= u{k}(1) <= 25, ...
                               -22 <= u{k}(2) <= 22];

    % State bounds applied to (x{k+1} - p{k}) (exactly your original form)
    constraints = [constraints, ...
        -10  <= x{k+1}(1) - p{k}(1) <= 10, ...
        -100 <= x{k+1}(2) - p{k}(2) <= 100, ...
        -4   <= x{k+1}(3) - p{k}(3) <= 4];
end

% Terminal tracking cost
objective = objective + 100*norm(Cd*x{N+1} - r{N+1}, 2);

% Build optimizer:
% input: current x and stacked reference r over horizon
% output: stacked u, x, p
parameters_in = {x{1}, [r{:}]};
solutions_out = {[u{:}], [x{:}], [p{:}]};

controller = optimizer(constraints, objective, sdpsettings('solver','mosek'), ...
                       parameters_in, solutions_out);

%% 4) Simulation parameters and reference
T      = 15;
steps  = floor(T/dT);
steps2 = steps + 500;  % ensure ref is long enough for i+N indexing

% Time vector (1 x (steps+1))
Time = zeros(1, steps+1);
for k = 1:steps+1
    Time(k) = (k-1)*dT;
end

% Reference trajectory ref is ny x steps2
ref = zeros(ny, steps2);
for k = 1:steps2
    if (k-1 < floor(steps/2))
        ref(:,k) = [9; 6.5];
    else
        ref(:,k) = [0; 0];
    end
end

%% 5) Closed-loop simulation with uncertainty + solve time measurement
xk = [0;0;0];

% Histories (preallocate)
Xhist  = zeros(nx, steps+1);   % for plotting (same as your original)
Uhist  = zeros(nu, steps+1);   % MPC decision input
Uthist = zeros(nu, steps+1);   % total input (here equals u, but keep for consistency)

fhist  = zeros(nu, steps+1);   % matched uncertainty history
whist  = zeros(1,  steps+1);   % unmatched scalar history

% State trajectory storage (avoid X = [X, x])
X = zeros(nx, steps+1);

solve_time_sum = 0;
solve_count    = 0;

for i = 1:steps+1
    % Log state
    Xhist(:,i) = xk;
    X(:,i)     = xk;

    % Horizon reference for this step
    future_r = ref(:, i:i+N);

    % Controller input
    input = {xk, future_r};

    % Define uncertainties (exactly your original formulas)
    f_uncert = 1.8*[-0.8*sin(0.4*pi*dT*(i-1)) - 0.1*xk(3)^2;
                     0.1 - 0.2*xk(3)];
    w_uncert = 3*sin(0.6*pi*dT*(i-1));

    fhist(:,i) = f_uncert;
    whist(:,i) = w_uncert;

    % ---- TIME ONLY THE CONTROLLER SOLVE ----
    t_start = tic;
    [solutions, diagnostics] = controller{input};
    solve_time_sum = solve_time_sum + toc(t_start);
    solve_count    = solve_count + 1;
    % ---------------------------------------

    % Extract first control move
    U = solutions{1};
    U = U(:);            % force column (prevents row/col implicit expansion issues)
    uk = [U(1); U(2)];
    uk = uk(:);          % ensure 2x1

    % In vanilla MPC, "total input" is just u itself
    U_total = uk;

    % Simulate actual system (nominal + matched + unmatched uncertainty)
    xk = Ad*xk + Bd*uk + Bd*f_uncert + Bu_d*w_uncert;

    % Log inputs
    Uhist(:,i)  = uk;
    Uthist(:,i) = U_total;
end

avg_solve_time = solve_time_sum / solve_count;
fprintf('Average controller solve time: %.6f seconds\n', avg_solve_time);

%% 6) Save results
save MPC_umpc.mat Uhist
save MPC_x.mat    X
save MPC_u.mat    Uthist

%% 7) Plots
close all;

figure;
plot(Time, Xhist(1,:) + Xhist(3,:)); hold on
plot(Time, ref(1,1:steps+1));
title('theta, x1+x3');

figure;
plot(Time, Xhist(1,:)); hold on
plot(Time, ref(2,1:steps+1));
title('gamma, x1');

figure;
plot(Time, X(1,:)); hold on
plot(Time, X(2,:)); hold on
plot(Time, X(3,:));
legend('x1','x2','x3')

figure;
plot(Time, X(3,:));
yline(4,'-','State Constrain');
yline(-4,'-','State Constrain');
title('x3');

figure;
plot(Time, Uhist(1,:));
title('de');

figure;
plot(Time, Uthist(1,:)); hold on
yline(25,'-','Input Bound');
yline(-25,'-','Input Bound');
title('total de')

figure;
plot(Time, Uhist(2,:));
title('df');

figure;
plot(Time, Uthist(2,:)); hold on
yline(22,'-','Input Bound');
yline(-22,'-','Input Bound');
title('total df');

figure;
plot(Time, fhist(1,:));
title('f 1');

figure;
plot(Time, fhist(2,:));
title('f 2');

figure;
plot(Time, whist(1,:));
title('w');