rho_0 =  0.1;          % bounds for ||x(0)||_inf
v_infnorm = 149.94;  %infinity norm of u_opt
w_infnorm = 1;
considering_x_csts = 1; 
% considering x constraints
state_to_focus = 3;  %1: gamma, 2: q, 3: alpha (default)
%state information: x[1] is gamma; x[2] is q ...
if state_to_focus == 1
    Tx = diag([1 0.01 0.01]); % for gamma
elseif state_to_focus == 2
    Tx = diag([0.01 1 0.01]);   % for q
elseif state_to_focus == 3
    Tx = diag([0.01 0.01 1]); % for alpha
else
    Tx = eye(3);
end

%%
Ae = L1_config.Ae;
Ta = L1_config.Ta; % sampling time
Cs = L1_config.Cs; %low pass filter
gamma1 = L1_config.gamma1;

Am = plant.Am;
B = plant.B;
Bv = plant.Bv; % Bv is Bv
In = plant.In;
Im = plant.Im;
Bu = plant.Bu;
%% coordinate transformation to emphasize the state (angle of attack) with constraints
% apply T^i_x to system for individual bounds
% x_new = Tx*x; u_new = Tu*u; f_new = Tu*f

Am_trans = Tx*Am/Tx; %eig_Am = eig(Am)
B_trans = Tx*B;
Bv_trans = Tx*Bv;
Cout_trans = plant.Cout/Tx;
Bu_trans = Tx*Bu;

% update the properties of uncertainty if necessary, according to f(t,Tx^{-1}x_new)
% note that the uncertainty depends only on angle of attack, which is not
% scaled. Therefore, for this example, there is no need to update the properties of uncertainty

%% determine rho_r in transformed system
s = tf('s');
rho_in = L1norm(s*Tx/(s*In-Am))*rho_0;

% transfer function from the inputs to the states
Hxm = (s*In-Am_trans)\B_trans;
Hxu = (s*In-Am_trans)\Bu_trans;
% transfer function from reference to the closed-loop states
Gxm = Hxm*(Im-Cs);
Hxv = (s*In-Am_trans)\Bv_trans; 

Gxm_l1norm = L1norm(Gxm);
Hxv_l1norm = L1norm(Hxv);
C_l1norm = L1norm(Cs);
Hxu_l1norm = L1norm(Hxu);

% solve rho_r
% update below as needed
 syms rho_r;
 if state_to_focus == 0
     bf_rhor  = (0.8+0.1*rho_r^2)*1.8; 
 end
 if considering_x_csts
    bf_rhor  = (0.8+0.1*(Xcon(3))^2)*1.8;
 end
b_w=1;
% use top one first to see if rho_r > Xcon; if >, then use the bottom one
% to constrain
 rho_r = solve(rho_r == (Gxm_l1norm*bf_rhor + Hxv_l1norm*v_infnorm+Hxu_l1norm*w_infnorm+rho_in)+0.01,rho_r);
 rho_r = vpa(rho_r(1),5) % save as 5 digits
% verify the stability condition
% rho_r = (Gxm_l1norm*bf + Hxv_l1norm*v_infnorm*1+rho_in)
% v_infnorm = (rho_r-rho_in-Gxm_l1norm*bf_rhor)/Hxv_l1norm
% if v_infnorm < 0
%      warning('Some states violate constraints in the presence of the specified uncertainty, bound on x0 and the reference.')
% end
if ~considering_x_csts
    rho = rho_r + gamma1;
    [~,Lf_rho] = cal_f_consts(rho);
else
    [~,Lf_rho] = cal_f_consts(Xcon(3));
end
if state_to_focus == 0 && ~considering_x_csts
    bf_rhor = (0.8+0.1*rho_r^2)*1.8;
    %bf_rhor  = 1.2+0.15*rho_r^2+2;
end

if Gxm_l1norm*Lf_rho>=1 % Lf = 0.2rho 
    disp('Stability condition (2) is not satisifed');
    % inscrease filter bandwidth to make it satisfy
end

%% refine the bounds on u_a
answer = questdlg('Did you run ensure Xr is correct?');
switch answer 
    case 'Yes'
    case 'No'
        error('Please verify the following Xr is correct!');
end 

% formulate Xr
Xr = [-22.898 22.898; -49.104 49.104; -1.8503 1.8503];
% refine Xr based on the constraints on alpha
Xr(1,:)= [max(Xr(1,1), -Xcon(1)); min(Xr(1,2),Xcon(1))];
Xr(2,:)= [max(Xr(2,1), -Xcon(2)); min(Xr(2,2),Xcon(2))];
Xr(3,:)= [max(Xr(3,1), -Xcon(3)); min(Xr(3,2),Xcon(3))];

bf_Xr = uncert_config.bf_vec_fcn(Xr(3,2));
rho_ur =  C_l1norm*bf_Xr;

fprintf('rho_in = %.3e, ||Gxm||_L1 =  %.3e, ||Hxv||_L1=%.3e, ||v||_inf = %.3e, rho_r = %.2e\n', ...
        rho_in,Gxm_l1norm, Hxv_l1norm, v_infnorm, rho_r);
fprintf('bf_rhor = %.3e,  rho_ur =[%.3e, %.3e] \n', ...
    bf_rhor, rho_ur(1), rho_ur(2));
% return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% assume the constraint on the sample time is satisfied; later, need to
% verify this. For verification, one needs to compute the constants \bar \alpha_i(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute \alpha_i_bar and evaluate whether sample time is small enough then find the error bound
bf_rho = bf_rhor; % when computing bf_rhor, we already consider the constraint set for X 
I = eye(3);
Phi = Ae\(expm(Ae*Ta)-I);
Phi_inv = inv(Phi);
syms tau
disp('start computing alpha_i_bar...');
tvec = linspace(0,Ta,10); N = length(tvec);
% f0 = @(tau) norm(expm(Ae.*(Ta-tau))*B,'inf');
% f2 = @(tau) norm(expm(Ae.*(Ta-tau))*Phi_inv*expm(Ae*Ta),'inf');
% f3 = @(tau) norm(expm(Ae.*(Ta-tau))*B,'inf')

f0 = @(tau) f0_fcn(tau,Ae,Ta,B);
f3 = @(tau) f3_fcn(tau,Ae,Ta,Bu);
f2 = @(tau) f2_fcn(tau,Ae,Ta,Phi_inv); 

% f3 = @(tau) f3_fcn(tau,Ae,Ta,B);
a0_bar = integral(f0,0,Ta);
a3_bar = integral(f3,0,Ta);
a1 = nan(1,N);
% a2 = nan(1,N);

for i = 1:N
    t = tvec(i);
    a1(i) = norm(expm(Ae*t),'inf');
%     a2(i) = 
%     a3(i) = integral(f3,0,t);
%     a2(i) = int(norm(expm(Ae*(Ta-tau))*Phi_inv*expm(Ae*Ta),'inf'),0,t);
%     a3(i) = int(norm(expm(Ae*(Ta-tau))*B,'inf'),0,t);
end
a1_bar = max(a1);
a2_bar = integral(f2,0,Ta); %max(a2);

gamma0_T = (bf_rho*a0_bar+a3_bar*b_w)*(a1_bar+a2_bar+1);

tmp  = L1norm(Hxm*Cs*pinv(B)*(s*I-Ae))/(1-Gxm_l1norm*Lf_rho)*gamma0_T-gamma1; % condition (42)
if tmp>=0
    tmp
    warning('Estimation sample time is not small enough')    
end
gamma2 = C_l1norm*Lf_rho*gamma1 + L1norm(Cs*pinv(B)*(s*In-Ae))*gamma0_T;

% difference between the nominal system and the adaptive system
rho_til_x = gamma1 + Gxm_l1norm*bf_rhor+Hxu_l1norm*b_w; % bound on ||xn-x||_inf
%rho_til_x = gamma1 + Gxm_l1norm*bf_Xr;
% rho_til_y = norm(Cout_trans,'inf')*rho_til_x;
rho_u = gamma2 + rho_ur;

fprintf('alpha_0_bar = %.2e, alpha_1_bar = %.2e, alpha_2_bar =%.2e\n',a0_bar,a1_bar,a2_bar);
fprintf('gamma0_T = %.2e, gamma1 = %.2e, gamma2 =%.2e, \nrho_til_x =  %.2e, rho_u1= %.2e, rho_u2= %.2e\n',gamma0_T,gamma1,gamma2,rho_til_x,rho_u(1),rho_u(2)); % rho_til_y = %.2e
%% Compute maximum v(0) such that the actual system is guaranteed to stay in X 
if state_to_focus == 3
    v0_infnorm = (Xcon-rho_til_x - rho_in)/Hxv_l1norm
end

% \alpha_bar_0
function rst = f0_fcn(tau,Ae,Ta,B) 
N = length(tau);
rst = zeros(1,N);
for i =1:N
    rst(i) = norm(expm(Ae.*(Ta-tau(i)))*B,'inf');
end
end
% \alpha_bar_2
function rst = f2_fcn(tau,Ae,Ta,Phi_inv) 
N = length(tau);
rst = zeros(1,N);
for i =1:N
    rst(i) = norm(expm(Ae.*(Ta-tau(i)))*Phi_inv*expm(Ae*Ta),'inf');
end
end
% \alpha_bar_3
function rst = f3_fcn(tau,Ae,Ta,Bu) 
N = length(tau);
rst = zeros(1,N);
for i =1:N
    rst(i) = norm(expm(Ae.*(Ta-tau(i)))*Bu,'inf');
end
end
% 
% function rst = f3_fcn(tau,Ae,Ta,B) 
% N = length(tau);
% rst = zeros(1,N);
% for i =1:N
%     rst(i) = norm(expm(Ae.*(Ta-tau(i)))*B,'inf');
% end
% end
function [bf,Lf]= cal_f_consts(rho)
Lf = max(0.2*rho,0.2)*1.8;
%Lf = max(0.3*rho,0.3);
bf = (0.8+0.1*rho^2)*1.8; 
%bf = 1.2+0.15*rho^2+2;
end

function bf_rho= bf_fcn(rho)
bf_rho = (0.8+0.1*rho^2)*1.8;
%bf_rho = 1.2+0.15*rho^2+2;
end
