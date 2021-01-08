% code for producing Figure 6.2

% 1. fixed the simulation runs for MC as N = 25*2^(16) and inner simulation runs for rQMC as n = 2^(16) 
% 2. fixed the toleranece level to 10^(-6)
% 3. initialized # of copies of sobol sequence with a rotation shift
% 4. initialized the vector of depth m is chosen to be 8, sufficiently large

sobol_num = 25;
run = sobol_num * 2^(16);
run_rqmc = 2^(16);
epsilon = 10^(-6);
M = 2:1:9;

% Initialized vector for storing
V = zeros(size(M,2),1);
V_rQMC = zeros(size(M,2),1);
se = zeros(size(M,2),1);
se_rQMC = zeros(size(M,2),1);

for i = 1:size(M,2)
    [V(i),se(i)] = Lookback_Option_ADGBS(M(i), run, epsilon);
    [V_rQMC(i),se_rQMC(i)] = Lookback_Option_ADGBS_RQMC(M(i), run_rqmc, sobol_num, epsilon);
end

% Gives the plot
figure(1)
plot(M, V, 'b-')
xlim([min(M)-0.1, max(M)+0.1])
xlabel('Depth(m)')
ylabel('Estimated Price')
title('Monte-Carlo Method')

figure(2)
plot(M, V_rQMC, 'r-')
xlim([min(M)-0.1, max(M)+0.1])
xlabel('Depth(m)')
ylabel('Estimated Price')
title('Randomized Quasi-Monte-Carlo Method')