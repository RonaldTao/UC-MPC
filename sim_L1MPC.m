function simRst = sim_L1MPC(Ae_d,BI_xhat_d,Uhist,X,x0, rHist,plant,uncert_config,baseline_config,L1_config,sim_config,Xcon_nom,Ucon_nom)
%% --- Simulate reference governor with an L1 adaptive controller in the lopp---
Am_d = plant.Am_d;
Am = plant.Am;
% Ap_d = plant.Ap_d;
Bv_d = plant.Bv_d;
B_d = plant.B_d;
Bu = plant.Bu;
Bu_d = plant.Bu_d;
Bv = plant.Bv;
B = plant.B;
B_pinv = pinv(B);           % pseudo-inverse of B
B_perp_pinv = pinv(plant.B_perp); % pseudo-inverse of B_perp
BI_d = plant.BI_d;
dT = sim_config.dT; %1e-4

ntimes = floor(plant.dT_sys/dT);
ntimes_adapt = floor(L1_config.Ta/dT);
f_uncert = uncert_config.f_uncert;
w_uncert = uncert_config.w_uncert;
add_uncert = uncert_config.add_uncert;

Kx = baseline_config.Kx;
Kv = baseline_config.Kv;

Ae = L1_config.Ae;
adapt_gain = L1_config.adapt_gain;

nLen = size(rHist,2);
nx = length(x0); % 3
m = size(B_d,2); % 2

umpcHist   = Uhist;
xHist   = zeros(nx,nLen);
tHist   = (0:1:nLen-1)*dT;
xcHist = xHist;  
uHist = zeros(m,nLen);
xCur  = x0;             % states of the actual system (discrete-time implementation)
xcCur = x0;             % states of the actual system (continuous-time implementation)
%xnomCur = x0;           % states of the nominal system (continuous-time implementation)
xrefCur = x0;           % states of the reference system (continuous-time implementation)
umpcCur  = umpcHist(:,1);

xHist(:,1) = x0;
%umpcHist(:,1) = umpc0;
% for L1 implementation
xhatHist = zeros(nx,nLen); 
xhatCur = x0;

filter = ss(L1_config.Cs);
filter_d = c2d(filter,dT);

% xfilterHist = zeros(m,nLen);
xfilterCur = zeros(m,1);
xfilterRefCur = zeros(m,1);

BsigmahatHist = zeros(nx,nLen);
% sigmahatHist = zeros(m,nLen);

sigmahat2Hist = zeros(nx-m,nLen);
sigmahatHist = zeros(m,nLen);

uL1Hist = zeros(m,nLen);

% for uncertainty
fHist = zeros(m,nLen);
wHist = zeros(1,nLen);
count = 1;
solve_time_sum = 0;   % accumulate controller solve time
solve_count    = 0;
for i=1:1:(nLen)
    %ref=rHist(:,i:i+1/plant.dT_sys);
    xHist(:,i) = xCur;
    xcHist(:,i) = xcCur;
    %xnomHist(:,i) = xnomCur;
    
    xhatHist(:,i) = xhatCur;
    % --------------- should use xnomCur instead of xcCur ----------------
%     Kp    =  findK(umpcCur, xcCur, rCur,  Am_d,  Bv_d, APd, bPd); % one
    % --------------------------------------------------------------------
    if mod(i-1,ntimes)==0
        umpcCur = umpcHist(:,count);
        count=count+1;
    end

    %umpcCur = umpcHist(:,i);
    %umpcHist(:,i) = umpcCur;
    t_start1 =tic;
    % ----------adaptive law ---------------
    if i==1 || mod(i-1,ntimes_adapt)==0
        xtilCur = xhatCur-xcCur;
        %***xtilCur = xhatCur-xCur;
        BsigmahatCur = -adapt_gain*xtilCur;        
    end
    BsigmahatHist(:,i) = BsigmahatCur;
    
    % -----------control law----------------
    sigmahatCur = B_pinv*BsigmahatCur;
    sigmahatHist(:,i) = sigmahatCur; 
    
    sigmahat2Cur = B_perp_pinv*BsigmahatCur;
    sigmahat2Hist(:,i) = sigmahat2Cur; 
    
%     u_bl = Kfb*xCur + Kff*umpcCur;
    if L1_config.include_L1 
        % compute L1 control signal
        uL1 = -(filter_d.c*xfilterCur + filter_d.d*sigmahatCur);      
    else
        uL1 = [0;0];
    end
    solve_time_sum = solve_time_sum + toc(t_start1);

    uL1Hist(:,i) = uL1;
    
%     %% for testing
%     uL1 = [0;0];    
    uBL = Kx*xcCur + Kv*umpcCur;
    %***uBL = Kx*xCur + Kv*umpcCur;
    uHist(:,i) = uL1+uBL;

    t = (i-1)*dT;
    if add_uncert == 1 
        fHist(:,i) = f_uncert(t,xcCur); 
        %***fHist(:,i) = f_uncert(t,xCur);  
    else
        fHist(:,i) = zeros(m,1);
    end
    wHist(:,i) = w_uncert(t);
    t_start2 = tic;
     % update state predictor
    %xhatCur = Ae_d*xhatCur + BI_xhat_d*(Bv*umpcCur+ Bv*uL1+ BsigmahatCur+(Am-Ae)*xCur);
    % continuous-time formulation + Euler integration
    xhatCur_dot = Ae*xhatCur + Bv*umpcCur+ B*uL1+ BsigmahatCur+(Am-Ae)*xcCur;
    %xhatCur_dot = Ae*xhatCur + Bv*umpcCur+ B*uL1+ BsigmahatCur+(Am-Ae)*xCur;
    xhatCur = xhatCur + xhatCur_dot*dT;
    solve_time_sum = solve_time_sum + toc(t_start2);
    solve_count = solve_count + 1;
    %%%%%%%%%%%%%% update system states %%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %There is a problem when using the discrete-time formulation for
    %simulation: probably because the state predictor uses a continous-time
    %formulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % discrete-time formulation
    if i==1 || mod(i-1,ntimes)==0
        xCur =  Am_d*xCur + BI_d*(Bv*umpcCur+ B*(uL1+fHist(:,i)));
    end
    
    % continuous-time formulation + Euler integration
    % actual states
    xcCur = xcCur+(Am*xcCur+Bv*umpcCur+ B*(uL1+fHist(:,i))+Bu*wHist(:,i))*dT;
    %xcCur = xcCur+(Am*xcCur+Bv*umpcCur+ B*(uL1+fHist(:,i)))*dT;   
    %xnomCur = xnomCur + (Am*xnomCur + Bv*umpcCur)*dT;
    % reference states
    %xrefCur = xrefCur + (Am*xrefCur + Bv*umpcCur + B*(fHist(:,i)-(filter_d.c*xfilterRefCur + filter_d.d*fHist(:,i)))+Bu*wHist(:,i))*dT;
    %***xrefCur = Am_d*xrefCur+BI_d*(Bv*umpcCur+B*(fHist(:,i)-(filter_d.c*xfilterRefCur + filter_d.d*fHist(:,i))));
    % update filter state
    xfilterCur = filter_d.a*xfilterCur + filter_d.b*sigmahatCur;    
    xfilterRefCur = filter_d.a*xfilterRefCur + filter_d.b*fHist(:,i);
end
avg_solve_time = solve_time_sum / solve_count;

fprintf('Average controller solve time: %.6f seconds\n', avg_solve_time);
% simRst.tHist = tHist;
simRst.umpcHist = umpcHist;
simRst.xHist = xHist;
simRst.BsigmahatHist = BsigmahatHist;
simRst.sigmahatHist = sigmahatHist;
simRst.uL1Hist = uL1Hist;
simRst.fHist = fHist;
simRst.xhatHist = xhatHist;
simRst.xcHist = xcHist;
simRst.uHist = uHist;
simRst.wHist = wHist;
disp('Simulation complete!');