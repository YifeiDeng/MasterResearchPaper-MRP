% Look back Option using Adaptive DGBS method: dynamically exclude those intervals (i.e., all enclosed subsequent grid points) from the partition that definitely neither contribute to the minimal nor the maximal values of the process path
%--------------------------------------------------------------------------
% INPUTS:
%
%   M:            Maximum Depth
%   run:          simulations run
%   epsilon:      the tolerance level
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% V: the estimated price of the floating strike lookback call option under
% the VG model within the Adpative DGBS framework
% se: the estimated standard error of the floating strike lookback call
% option under the VG model within the Adpative DGBS framework
%--------------------------------------------------------------------------
function [V,se] = Lookback_Option_ADGBS(M, run, epsilon)
T = 2^M; % set index T for generating the process
T_ = 0.40504; % set maturity
v = 0.2505; % the volatility parameter for gamma process
theta = -0.2859; % the drift parameter for gamma process
sigma = 0.1927; % the volatility for VG
r = 0.0548; % risk free rate
q = 0; % dividend paying
mu = r + 1/v * log(1-theta*v-1/2*sigma^2*v)-q;
S_0 = 100;

up = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) + theta/2;
un = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) - theta/2;
vp = up^2 * v;
vn = un^2 * v;
% mu_p = max(0,mu);
mu_n = min(0,mu);

G_p = zeros(T+1,1);
G_n = zeros(T+1,1);
G_p(1) = 0;
G_n(1) = 0;

t = zeros(T,1);
t(1) = T_; % t_1 = T
t = [0;t]; % t_0 = 0

wow = zeros(run,1);
path = zeros(1,run);
po = zeros(1,run);
means = 0;
squares = 0;

for N = 1:run
    G_p(2) = gaminv(rand,up^2*T_/vp, vp/up);
    G_n(2) = gaminv(rand,un^2*T_/vn, vn/un);
    X = G_p(2) - G_n(2) + T_*mu;
    go = false;
    bl_max = max(0,X); su_min = min(0,X);
    lint(1) = 1;rint(1) = 2;nint = 1;pos = 2;
    for m = 1:M
        lint_c = lint; rint_c = rint; nint_c = nint; nint = 0;
        sl = bl_max;
        path(N) = path(N) + 1;
        if go == true
            break;
        else
            for j = 1:nint_c
                pos = pos + 1;
                i_l = lint_c(j); i_m = pos; i_r = rint_c(j); t(i_m) = (t(i_r)+t(i_l))/2;delta = (t(i_r)-t(i_l))/(2*v);
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
                i_n = [i_l,i_m];
                i_p = [i_m,i_r];
                for w = 1:2
                    lb = G_p(i_n(w))-G_n(i_p(w))+ t(i_n(w))*mu+(t(i_p(w))-t(i_n(w)))*mu_n;
                    if lb <= su_min
                        nint = nint + 1; lint(nint) = i_n(w); rint(nint) = i_p(w);
                        sl = min(sl,lb);
                    end
                end
            end
        end
        if exp(-r*T_)*S_0*(exp(su_min)-exp(sl)) <= 2*epsilon
            go = true;
            break
        end
    end
    if mod(N,10000) == 0
        fprintf('%d th rounds.\n',N);
    end
%     index = find(~t);
%     wow(N) = index(2)-1;
    po(N) = exp(-r*T_)*S_0*(exp(X)-1/2*(exp(su_min)+exp(sl)));
    delta = po(N) - means;
    means = means + delta/N;
    squares = squares + delta*(po(N)-means);
end
V = 1/run*sum(po);
se = sqrt(squares/(run*(run-1)));

end

% tic
% Lookback_Option_ADGBS
% toc
