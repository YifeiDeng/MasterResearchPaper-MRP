% Apply Adaptive DGBS method on the Barrier down and out call Option using
% rQMC
%--------------------------------------------------------------------------
% INPUTS:
%
%   M:            Maximum Depth
%   run:          inner simulations runs for rQMC
%   sobol_num:    #s of copies of Sobol sequence with a rotation shift
%   B_l:          the lower barrier for the down-and-out option
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% V: the estimated price of the down-and-out barriar call option under
% the VG model within the Adpative DGBS framework using rQMC
% se: the estimated standard error of the down-and-out barriar call option
% under the VG model within the Adpative DGBS framework using rQMC

function [V,se] = Barrier_Option_ADGBS_Call_downout_RQMC(M, run, sobol_num, B_l)
T_ = 0.46575; % set maturity
v = 0.49083; % the volatility parameter for gamma process
theta = -0.28113; % the drift parameter for gamma process
sigma = 0.19071; % the volatility for VG
r = 0.0549;
q = 0.011;
S_0 = 100;
K = 100;
% B_u = 105;

T = 2^M; % set index T for generating the process
t = zeros(T,1);
t(1) = T_; % t_1 = T
t = [0;t]; % t_0 = 0

mu = r + 1/v * log(1-theta*v-1/2*sigma^2*v)-q;
up = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) + theta/2;
un = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) - theta/2;
vp = up^2 * v;
vn = un^2 * v;

mean_B = 0;
squares_B = 0;

mu_p = max(0,mu);
mu_n = min(0,mu);

po_B = zeros(run,1);

G_p = zeros(T+1,1);
G_n = zeros(T+1,1);
G_p(1) = 0;
G_n(1) = 0;

call_M = zeros(sobol_num,1);

for k = 1:sobol_num
    qmc_set_orignial = net(sobolset(2*T+2),run);
    r_direct = unifrnd(0,1,[1,2*T+2]);
    shifted_qmc = mod(qmc_set_orignial+repelem(r_direct,run,1),1);
    call = zeros(run,1);
    for N = 1:run
        G_p(2) = gaminv(shifted_qmc(N,1),up^2*T_/vp, vp/up);
        G_n(2) = gaminv(shifted_qmc(N,2),un^2*T_/vn, vn/un);
        X = G_p(2) - G_n(2) + T_*mu;
        bl_max = max(0,X); su_min = min(0,X);
        lint(1) = 1;rint(1) = 2;nint = 1;pos = 2;
        in = false;
        % criterias for down-out call option
        if S_0*exp(X)>=K && S_0*exp(-G_n(2)-T_*mu_n) >= B_l
            in = true;  done = false;
            for m = 1:M
                lint_c = lint; rint_c = rint; nint_c = nint; nint = 0;
                for j = 1:nint_c
                    pos = pos + 2;
                    i_l = lint_c(j); i_m = pos-1; i_r = rint_c(j); t(i_m) = (t(i_r)+t(i_l))/2;delta = (t(i_r)-t(i_l))/(2*v);
                    % simulating the beta random variate
                    Yp = betainv(shifted_qmc(N,pos-1),delta, delta);
                    G_p(i_m) = G_p(i_l) + (G_p(i_r) - ...
                        G_p(i_l))*Yp;
                    % simulating the beta random variate
                    Yn = betainv(shifted_qmc(N,pos),delta, delta);
                    G_n(i_m) =  G_n(i_l) + (G_n(i_r) - ...
                        G_n(i_l))*Yn;
                    su_min = min([su_min, G_p(i_m) - G_n(i_m) + t(i_m)*mu]);
                    bl_max = max([bl_max, G_p(i_m) - G_n(i_m) + t(i_m)*mu]);
                    i_n = [i_l,i_m];
                    i_p = [i_m,i_r];
                    for w = 1:2
                        lb = G_p(i_n(w))-G_n(i_p(w))+ t(i_n(w))*mu+(t(i_p(w))-t(i_n(w)))*mu_n;
                        ub = G_p(i_p(w))-G_n(i_n(w))+ t(i_n(w))*mu+(t(i_p(w))-t(i_n(w)))*mu_p;
                        % If exceeding the lower barrier, option worthless,
                        % change in = false
                        if S_0*exp(su_min) < B_l
                            in = false;
                            done = true;
                            break
                        end
                        
                        % if above does not meet, then no yet definite and continue the iterations
                        if S_0*exp(lb) < B_l
                            nint = nint + 1; lint(nint) = i_n(w); rint(nint) = i_p(w);
                        end
                    end
                end
                % if done = true (meet either knock out or knock in conditions)
                % or nint == 0, then break the loop
                if done == true || nint == 0
                    break
                end
            end
        end
        
        % if the knock-out condition meet, then worthless
        if S_0*exp(su_min) < B_l
            in = false;
        end
        
        if in == true
            call(N) = exp(-r*T_)*max([0,S_0*exp(X)-K]); % call
        elseif in == false
            call(N) = 0;
        end
        
        if mod(N,10000) == 0
            fprintf('%d th rounds.\n',N);
        end
    end
    if mod(k,5) == 0
        fprintf('%d th iterations.\n',k);
    end
    call_M(k) = sum(call)/run;
end

se = sqrt(var(call_M)/(size(call_M,1)));
V = 1/size(call_M,1) * sum(call_M);
end
