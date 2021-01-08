% Apply Adaptive DGBS method on the up-and-in call Barrier Option
%--------------------------------------------------------------------------
% INPUTS:
%
%   M:            Maximum Depth
%   run:          simulations run
%   B_u:          the upper barrier for the up-and-in option
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% V: the estimated price of the up-and-in barriar call option under
% the VG model within the Adpative DGBS framework
% se: the estimated standard error of the up-and-in barriar call option
% under the VG model within the Adpative DGBS framework
%--------------------------------------------------------------------------
function [V,se] = Barrier_Option_ADGBS_Call_upin(M, run, B_u) 
T_ = 0.46575; % set maturity
v = 0.49083; % the volatility parameter for gamma process
theta = -0.28113; % the drift parameter for gamma process
sigma = 0.19071; % the volatility for VG
r = 0.0549;
q = 0.011;
S_0 = 100;
K = 100;
% B_l = 80;

T = 2^M; % set index T for generating the process
t = zeros(T,1);
t(1) = T_; % t_1 = T
t = [0;t]; % t_0 = 0

mu = r + 1/v * log(1-theta*v-1/2*sigma^2*v)-q;
up = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) + theta/2;
un = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) - theta/2;
vp = up^2 * v;
vn = un^2 * v;

mu_p = max(0,mu);
mu_n = min(0,mu);

po_B = zeros(run,1);
path = zeros(1,run);

G_p = zeros(T+1,1);
G_n = zeros(T+1,1);
G_p(1) = 0;
G_n(1) = 0;

rand_unif_variate = rand(run,2);

for N = 1:run
    G_p(2) = gaminv(rand_unif_variate(N,1),up^2*T_/vp, vp/up);
    G_n(2) = gaminv(rand_unif_variate(N,2),un^2*T_/vn, vn/un);
    X = G_p(2) - G_n(2) + T_*mu;
    bl_max = max(0,X); su_min = min(0,X);
    lint(1) = 1;rint(1) = 2;nint = 1;pos = 2;
    % criterias for up-in(po = (S_T-K)+)/up-out(po = 0)
    if S_0*exp(X)>=K && S_0*exp(G_p(2) + T_*mu_p) > B_u
        in = true;  done = false;
        for m = 1:M
            lint_c = lint; rint_c = rint; nint_c = nint; nint = 0;
            path(N) = path(N) + 1;
            for j = 1:nint_c
                pos = pos + 1;
                i_l = lint_c(j); i_m = pos; i_r = rint_c(j); t(i_m) = (t(i_r)+t(i_l))/2;delta = (t(i_r)-t(i_l))/(2*v);
                % simulating the beta random variate
                if delta > 0 && delta <= 1/2
                    U = rand;
                    V = rand;
                    if rand<1/2
                        S = 1;
                    else
                        S = 0;
                    end
                    Yp = 1/2 + S/(2*sqrt(1+1/(U^((-1/delta)-1)*(cos(2*pi*V))^2)));
                elseif delta > 1/2
                    S = 2;
                    while S > 1
                        U2 = unifrnd(0,1);
                        V2 = unifrnd(-1,1);
                        S = U2^2 + V2^2;
                    end
                    Yp = 1/2+U2*V2/S*sqrt(1-S^(2/(2*delta-1)));
                end
                G_p(i_m) = G_p(i_l) + (G_p(i_r) - ...
                    G_p(i_l))*Yp;
                % simulating the beta random variate
                if delta > 0 && delta <= 1/2
                    U = rand;
                    V = rand;
                    if rand<1/2
                        S = 1;
                    else
                        S = 0;
                    end
                    Yn = 1/2 + S/(2*sqrt(1+1/(U^((-1/delta)-1)*(cos(2*pi*V))^2)));
                elseif delta > 1/2
                    S = 2;
                    while S > 1
                        U2 = unifrnd(0,1);
                        V2 = unifrnd(-1,1);
                        S = U2^2 + V2^2;
                    end
                    Yn = 1/2+U2*V2/S*sqrt(1-S^(2/(2*delta-1)));
                end
                G_n(i_m) =  G_n(i_l) + (G_n(i_r) - ...
                    G_n(i_l))*Yn;
                su_min = min([su_min, G_p(i_m) - G_n(i_m) + t(i_m)*mu]);
                bl_max = max([bl_max, G_p(i_m) - G_n(i_m) + t(i_m)*mu]);
                i_n = [i_l,i_m];
                i_p = [i_m,i_r];
                for w = 1:2
                    lb = G_p(i_n(w))-G_n(i_p(w))+ t(i_n(w))*mu+(t(i_p(w))-t(i_n(w)))*mu_n;
                    ub = G_p(i_p(w))-G_n(i_n(w))+ t(i_n(w))*mu+(t(i_p(w))-t(i_n(w)))*mu_p;
                    % up-in, has value, do not change in = true
                    if S_0*exp(bl_max) > B_u
                        in = true;
                        done = true;
                        break
                    end
                    % if above does not meet, then continue the iterations
                    if S_0*exp(ub) > B_u
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
    % if the stock price does not exceed upper barrier, the contract is worthless, set in = false
    if S_0*exp(bl_max) <= B_u
        in = false;
    end
    
    if in == true
        po_B(N) = exp(-r*T_)*max([0,S_0*exp(X)-K]); % call
    else
        po_B(N) = 0;
    end
    if mod(N,10000) == 0
        fprintf('%d th rounds.\n',N);
    end    
end

se = sqrt(var(po_B)/size(po_B,1));
V = 1/size(po_B,1) * sum(po_B);

end
