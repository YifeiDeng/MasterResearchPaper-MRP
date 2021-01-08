% Brownian-Gamma Bridge Sampling (BGBS) of a VG process X(t) = B(G(t; mu, v), theta , sigma) for a 2k -Point Equal-Length Partition of [0, T]
k = 10; % the maximum depth
T = 2^k; % set index T for generating the process
T_ = 0.40504; % set maturity to .40504 years
mu = 0.31356; % the drift for VG
v = 0.2505; % the volatility parameter for gamma process
theta = -0.2859; % the drift parameter for gamma process
sigma = 0.1927; % the volatility for VG
t = linspace(0, T_, T+1);
G = zeros(T+1,100);
X = zeros(T+1,100);
G(1,:) = 0;
G(end,:) = gamrnd(T_/v, v);
X(1,:) = 0;
X(end,:) = normrnd(theta*G(end,1), sqrt(sigma^2*G(end,1)));

for j = 1:100
    for l = 1:k
        for m = 1:2^(l-1)
            i = 2*m - 1;
            Y = betarnd(T_*T/(v*2^l), T_*T/(v*2^l));
            G(i*T/2^l+1,j) = G((i-1)*T/2^l+1,j) + (G((i+1)*T/2^l+1,j) - ...
                G((i-1)*T/2^l+1,j))*Y;
            b = G((i+1)*T/2^l+1,j) - G(i*T/2^l+1,j);
            z = normrnd(0, sqrt(b*sigma^2*Y));
            X(i*T/2^l+1,j) = Y * X((i+1)*T/2^l+1,j) + ...
                (1-Y)*X((i-1)*T/2^l+1,j) + z;
        end
    end
    X(:,j) = X(:,j) + (mu*t)';
    hold on
    plot(linspace(0, T_, T+1),X(:,j))
    hold off
end

title('100 Simus of a VG process: (mu, sigma, v, theta)=(.31,.19,.25,-.29)(BGBS)')
axis auto
xlabel('Times (years)')
ylabel('Variance Gamma Process X(t)')
