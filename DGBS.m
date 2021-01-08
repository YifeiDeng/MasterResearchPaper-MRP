% Truncated Double Gamma Bridge Sampling (DGBS) of a VG Process X(t) = B(G(t; mu, v), theta , sigma) for a 2k -Point Equal-Length Partition of [0, T ]
%--------------------------------------------------------------------------
% INPUTS:
%
%   M:            Maximum Depth
%   run:          simulations run
%--------------------------------------------------------------------------
%
% OUTPUT:
%
% X: realizations of the VG processes for run simulations
%--------------------------------------------------------------------------
function [X] = DGBS(M, run)
T = 2^M; % set index T for generating the process
T_ = 5; % set maturity
mu = .31; % the drift for VG
v = .25; % the volatility parameter for gamma process
theta = -.29; % the drift parameter for gamma process
sigma = .19; % the volatility for VG

up = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) + theta/2;
un = 1/2*(sqrt(theta^2 + 2*sigma^2/v)) - theta/2;
vp = up^2 * v;
vn = un^2 * v;

t = linspace(0, T_, T+1);

G_p = zeros(T+1,1);
G_n = zeros(T+1,1);

X = zeros(T+1,run);
path = zeros(1,run);

G_p(1) = 0;
G_n(1) = 0;


for j = 1:run
    G_p(end) = gaminv(rand,up^2*T_/vp, vp/up);
    G_n(end) = gaminv(rand,un^2*T_/vn, vn/un);
    for l = 1:M
        path(j) = path(j)+1;
        for m = 1:2^(l-1)
            i = 2*m - 1;
            delta_1 = up^2*T_/(vp*2^l);
            if delta_1 > 0 && delta_1 <= 1/2
                U = rand;
                V = rand;
                if rand<1/2
                    S = 1;
                else
                    S = 0;
                end
                Yp = 1/2 + S/(2*sqrt(1+1/(U^((-1/delta_1)-1)*(cos(2*pi*V))^2)));
            elseif delta_1 > 1/2
                S = 2;
                while S > 1
                    U2 = unifrnd(0,1);
                    V2 = unifrnd(-1,1);
                    S = U2^2 + V2^2;
                end
                Yp = 1/2+U2*V2/S*sqrt(1-S^(2/(2*delta_1-1)));
            end
            G_p(i*T/2^l+1) = G_p((i-1)*T/2^l+1) + (G_p((i+1)*T/2^l+1) - ...
                G_p((i-1)*T/2^l+1))*Yp;
            delta_2 = un^2*T_/(vn*2^l);
            if delta_2 > 0 && delta_2 <= 1/2
                U = rand;
                V = rand;
                if rand<1/2
                    S = 1;
                else
                    S = 0;
                end
                Yn = 1/2 + S/(2*sqrt(1+1/(U^((-1/delta_2)-1)*(cos(2*pi*V))^2)));
            elseif delta_2 > 1/2
                S = 2;
                while S > 1
                    U2 = unifrnd(0,1);
                    V2 = unifrnd(-1,1);
                    S = U2^2 + V2^2;
                end
                Yn = 1/2+U2*V2/S*sqrt(1-S^(2/(2*delta_2-1)));
            end
            G_n(i*T/2^l+1) =  G_n((i-1)*T/2^l+1) + (G_n((i+1)*T/2^l+1) - ...
                G_n((i-1)*T/2^l+1))*Yn;
            X(i*T/2^l+1,j) = G_p(i*T/2^l+1) - G_n(i*T/2^l+1);
        end
    end
    if mod(j,10000) == 0
        fprintf('%d th rounds.\n',j);
    end
    X(:,j) = X(:,j) + (mu*t)';
    X(end,j) = G_p(end) - G_n(end) + T_*mu;
%     hold on
%     plot(t,X(:,j))
%     hold off
end


% title('100 Simus of a VG process: (mu, sigma, v, theta)=(.31,.19,.25,-.29)(DGBS)')
% axis auto
% xlabel('Times (years)')
% ylabel('Variance Gamma Process X(t)')
end

