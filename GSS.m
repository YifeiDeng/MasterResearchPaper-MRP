% GSS Brownian-Gamma Beidge Sampling  of a gamma process (i.e. G(t; u, v)) 
% for a $2^{-k}$-point Equal-Length Partition of [0,T]

k = 10; % 2^10 equal length partitions of [0,T]
T = 2^k;
T_= 5; % set maturity to 5 years
h = 2^(-k)*T_;
G = zeros(T+1,100);
G(1,:) = 0;
t = linspace(0, T_, T+1);

mu = 1; % the drift parameter
v = 0.656; % the volatility parameter for gamma process

for j = 1:100
    for i = 2:(2^k+1)
        Q = gamrnd(mu^2*h/v, v/mu);
        G(i,j) = G(i-1,j) + Q;
    end
    hold on
    plot(t,G(:,j))
    hold off
end
title('100 Simulations of a gamma process with mu = 1, v = 0.656 (GSS)')
axis auto
xlabel('Time (years)')
ylabel('Time Transfromed (years)')


