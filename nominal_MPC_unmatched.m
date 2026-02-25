%% Nominal MPC with tightened constraints (YALMIP + MOSEK)
% This script:
%   1) Builds a discrete-time nominal closed-loop model (Am = Ac + Bc*Kx)
%   2) Solves an MPC problem on the nominal model with tightened bounds
%   3) Simulates the closed-loop system and saves u / x histories
%
% IMPORTANT:
%   - The MPC decision variable is u{k}
%   - The "total" control applied for constraint checking is u_total = u{k} + Kx*x{k}

%% 0) Environment setup
% add path below for YALMIP and solver
addpath(genpath('C:\Users\xxx\matlab\YALMIP'))
addpath(genpath('C:\Program Files\Mosek\10.1\toolbox\r2017a'))
yalmip('clear')
clear all
clc

% Load tightened constraint bounds
% Expected:
%   Ucon_nom : input bounds 
%   Xcon_nom : state bounds
load uc.mat;
load xc.mat;

% Make sure bounds are column vectors to avoid dimension mismatch
Ucon_nom = Ucon_nom(:);
Xcon_nom = Xcon_nom(:);

%% 1) Continuous-time model definition
Ac = [0 0.0067 1.34;
      0 -0.869 43.2;
      0 0.993 -1.34];

Bc = [0.169  0.252;
      -17.3  -1.58;
      -0.169 -0.252];

% Unmatched input direction (not used later, but kept for reference)
B_perp = null(Bc');
Bu = B_perp*0.15;

% State feedback gain (used only to define Am and "total input")
Kx = [3.25  0.891  7.12;
     -6.1  -0.898 -10];

% Nominal closed-loop continuous-time A matrix
Am = Ac + Bc*Kx;

% Output model (tracking outputs)
Cc = [1 0 1;
      1 0 0];
Dc = [0 0;
      0 0];

nx = 3;     % number of states
nu = 2;     % number of inputs
ny = 2;     % number of outputs

%% 2) Discretization
dT = 4e-3;  % sampling time (seconds)

% Discretize the closed-loop nominal model with input channels [Bc Bu]
% Resulting discrete system:
%   x_{k+1} = Am_d*x_k + Bd*u_k + Bu_d*w_k
% Here we only use Bd and decision u_k
[Am_d, Bd_total, Cd, Dd] = ssdata(c2d(ss(Am, [Bc Bu], Cc, [Dc [0;0]]), dT));
Bd   = Bd_total(:,1:2);  % matched control channel
Bu_d = Bd_total(:,3);    % unmatched

%% 3) MPC setup
N = 0.2/dT;         % horizon length in steps (should be integer)
N = round(N);       % force integer to avoid sdpvar issues

% Decision variables:
%   u{k} : nu x 1, k=1..N
%   x{k} : nx x 1, k=1..N+1
%   r{k} : ny x 1, k=1..N+1 (reference trajectory over horizon)
u = sdpvar(repmat(nu,1,N),   repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1), repmat(1,1,N+1));
r = sdpvar(repmat(ny,1,N+1), repmat(1,1,N+1));

constraints = [];
objective   = 0;

% Build objective + constraints
for k = 1:N
    % Tracking output at step k
    yk = Cd*x{k};

    % Stage cost
    if k == 1
        objective = objective + 100*norm(yk - r{k}, 2) + 1*norm(u{k}, 2);
    else
        objective = objective + 100*norm(yk - r{k}, 2) + 1*norm(u{k}, 2) ...
                              + 100*norm(u{k} - u{k-1}, 2);
    end

    % Nominal dynamics
    constraints = [constraints, x{k+1} == Am_d*x{k} + Bd*(u{k})];

    % Total input used for constraint checking: u_total = u + Kx*x
    constraints = [constraints, ...
        -Ucon_nom(1) <= Kx(1,1)*x{k}(1) + Kx(1,2)*x{k}(2) + Kx(1,3)*x{k}(3) + u{k}(1) <= Ucon_nom(1), ...
        -Ucon_nom(2) <= Kx(2,1)*x{k}(1) + Kx(2,2)*x{k}(2) + Kx(2,3)*x{k}(3) + u{k}(2) <= Ucon_nom(2), ...
        -Xcon_nom(1) <= x{k+1}(1) <= Xcon_nom(1), ...
        -Xcon_nom(2) <= x{k+1}(2) <= Xcon_nom(2), ...
        -Xcon_nom(3) <= x{k+1}(3) <= Xcon_nom(3)];
end

% Terminal tracking cost
objective = objective + 100*norm(Cd*x{N+1} - r{N+1}, 2);

% Build YALMIP optimizer object:
% Input to optimizer: current state x{1} and stacked reference [r{:}]
% Output: stacked decision inputs [u{:}] and states [x{:}]
parameters_in  = {x{1}, [r{:}]};
solutions_out  = {[u{:}], [x{:}]};

controller = optimizer(constraints, objective, sdpsettings('solver','mosek'), ...
                       parameters_in, solutions_out);

%% 4) Simulation parameters and reference construction
T      = 15;                 % simulation duration (seconds)
steps  = floor(T/dT);         % number of simulation steps
steps2 = steps + 500;         % extra for horizon look-ahead indexing

% Time vector (1 x (steps+1))
Time = zeros(1, steps+1);
for k = 1:steps+1
    Time(k) = (k-1)*dT;
end

% Reference trajectory ref is ny x steps2
% First half: [9; 6.5], second half: [0; 0]
ref = zeros(ny, steps2);
for k = 1:steps2
    if (k-1 < floor(steps/2))
        ref(:,k) = [9; 6.5];
    else
        ref(:,k) = [0; 0];
    end
end

%% 5) Closed-loop simulation
xk = [0;0;0];     % initial state

% Histories
Xhist  = zeros(nx, steps+1);   % state history (for plotting)
Uhist  = zeros(nu, steps+1);   % applied MPC decision input u
Uthist = zeros(nu, steps+1);   % total input u + Kx*x (for constraint checking)

X = zeros(nx, steps+1); % preallocate to avoid growing inside loop

for i = 1:steps+1
    % Save current state
    Xhist(:,i) = xk;
    X(:,i)     = xk;

    % Build horizon reference (needs ref long enough for i+N)
    future_r = ref(:, i:i+N);

    % Call the optimizer
    input = {xk, future_r};
    [solutions, diagnostics] = controller{input};

    % Extract optimizer output
    % NOTE: solutions{1} is stacked [u{1}; u{2}; ...; u{N}] (size 2N x 1 or 1 x 2N)
    U = solutions{1};
    U = U(:);                 % force column to avoid implicit expansion

    uk = [U(1); U(2)];        % first control move
    uk = uk(:);               % force 2x1

    % Total input for constraint checking / logging
    U_total = uk + Kx*xk;     % guaranteed 2x1

    % Propagate nominal discrete-time dynamics
    xk = Am_d*xk + Bd*uk;

    % Log inputs
    Uhist(:,i)  = uk;
    Uthist(:,i) = U_total;
end

%% 6) Save results
save umpchist_nominal.mat   Uhist
save xhist_nominal.mat      X
save utotalhist_nominal.mat Uthist

%% 7) Plots
close all;

figure;
plot(Time, Xhist(1,:) + Xhist(3,:)); hold on
plot(Time, ref(1,1:steps+1));
title('theta, x1 + x3');

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
plot(Time, Uthist(1,:));
yline(25,'-','Input Bound');
yline(-25,'-','Input Bound');
title('total de')

figure;
plot(Time, Uhist(2,:));
title('df');

figure;
plot(Time, Uthist(2,:));
yline(22,'-','Input Bound');
yline(-22,'-','Input Bound');
title('total df');