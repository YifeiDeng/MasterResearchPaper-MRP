% Adaptive DGBS method: dynamically exclude those intervals (i.e., all enclosed subsequent grid points) from the partition that definitely neither contribute to the minimal nor the maximal values of the process path
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
% x_l: the lower bound for the VG processes for run simulations
% x_u: the upper bound for the VG processes for run simulations
% point: the average points simulated in each path
%--------------------------------------------------------------------------
function [x_l, x_u, point] = Adaptive_DGBS(M, run, epsilon)
T = 2^M; % set index T for generating the process
T_ = 0.40504;
mu = 0.31356; % the drift for VG
v = 0.2505; % the volatility parameter for gamma process
theta = -0.2859; % the drift parameter for gamma process
sigma = 0.1927; % the volatility for VG
up = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) + theta/2;
un = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) - theta/2;
vp = up^2 * v;
vn = un^2 * v;
mu_p = max(0,mu);
mu_n = min(0,mu);

G_p = zeros(T+1,1);
G_n = zeros(T+1,1);

t = zeros(T,1);
t(1) = T_; % t_1 = T
t = [0;t]; % t_0 = 0

x_l = zeros(run,1);
x_u = zeros(run,1);

X_ = zeros(T+1,1);
wow = zeros(run,1);
point = zeros(1,run);

G_p(1) = 0;
G_n(1) = 0;

for N = 1:run
    G_p(2) = gaminv(rand,up^2*T_/vp, vp/up);
    G_n(2) = gaminv(rand,un^2*T_/vn, vn/un);
    X = G_p(2) - G_n(2) + T_*mu;
    bl_max = max(0,X); su_min = min(0,X);
    lint(1) = 1;rint(1) = 2;nint = 1;pos = 2;
    go = false;
    for m = 1:M
        lint_c = lint; rint_c = rint; nint_c = nint; nint = 0;
        sl = bl_max; bu = su_min;
        if go == true
            return
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
                X_(i_m) = G_p(i_m) - G_n(i_m);
                su_min = min([su_min, G_p(i_m) - G_n(i_m) + t(i_m)*mu]);
                bl_max = max([bl_max, G_p(i_m) - G_n(i_m) + t(i_m)*mu]);
                i_n = [i_l,i_m];
                i_p = [i_m,i_r];
                for w = 1:2
                    lb = G_p(i_n(w))-G_n(i_p(w))+ t(i_n(w))*mu+(t(i_p(w))-t(i_n(w)))*mu_n;
                    ub = G_p(i_p(w))-G_n(i_n(w))+ t(i_n(w))*mu+(t(i_p(w))-t(i_n(w)))*mu_p;
                    if lb <= su_min || ub >= bl_max
                        nint = nint + 1; lint(nint) = i_n(w); rint(nint) = i_p(w);
                        sl = min(sl,lb);
                        bu = max(bu,ub);
                    end
                end
            end
            if max(su_min-sl, bu-bl_max) < (2*epsilon)
                go = true;
                return
            end
        end
    end
    index = find(~t);
    if mod(N,10000) == 0
        fprintf('%d th rounds.\n',N);
    end
    point(N) = index(2)-1;
    wow(N) = index(2)-1;
    X_(:) = X_(:) + (mu*t);
    X_(2) = G_p(2) - G_n(2) + T_*mu;
    x_l(N) = 1/2*(su_min + sl);
    x_u(N) = 1/2*(bl_max + bu);
end
% whole = sort([t(1:wow(10)),X_(1:wow(10),1)]);
% hold on
% plot(whole(:,1),whole(:,2))
% plot(whole(:,1),repelem(x_l(10),size(whole(:,1),1)),'b-')
% plot(whole(:,1),repelem(x_u(10),size(whole(:,1),1)),'r-')
% hold off
end

