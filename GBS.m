% Gamma Bridge Sampling of a Process G(t; µ, ?) for a 2k -Point 
% Equal-Length Partition of [0, T]

k = 10; % the maximum depth
T = 2^k; % set index T for generating the process
T_ = 5; % set maturity to 5 years
mu = 1; % the drift parameter
v =  0.656; % the volatility parameter for gamma process
t = linspace(0, T_, T+1);
G = zeros(T+1,100);
G(1,:) = 0;


for j = 1:100
    G(end,:) = gamrnd(mu^2*T_/v, v/mu);
    for l = 1:k
        for m = 1:2^(l-1)
            i = 2*m - 1;
            Y = betarnd(mu^2*T_/(v*2^l), mu^2*T_/(v*2^l));
            G(i*T/2^l+1,j) = G((i-1)*T/2^l+1,j) + (G((i+1)*T/2^l+1,j) - ...
                G((i-1)*T/2^l+1,j))*Y;
        end
    end
    hold on
    plot(t,G(:,j))
    hold off
end

title('100 Simulations of a gamma process with mu = 1, v = 0.656 (GBS)')
axis auto
xlabel('T (maturity)')
ylabel('Time Transfromed (years)')