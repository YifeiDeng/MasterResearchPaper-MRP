% code for producing Table 6.3

% 1. the depth m is chosen to be 8, sufficiently large
% 2. initialized vector of toleranece
% 3. initialized simulation runs for both MC and rQMC 
% 4. initialized # of copies of sobol sequence with a rotation shift
M = 8;
epsilon = [10^(0),10^(-2),10^(-4),10^(-6),10^(-8),10^(-10)];
run = 25*2^(16);
run_rqmc = 2^(16);
sobol_num = 25;

% Initialized vector for storing
V = zeros(size(epsilon,2),1);
V_rQMC = zeros(size(epsilon,2),1);
se = zeros(size(epsilon,2),1);
se_rQMC = zeros(size(epsilon,2),1);

for i = 1:size(epsilon,2)
    [V(i),se(i)] = Lookback_Option_ADGBS(M, run, epsilon(i));
%     [V_rQMC(i),se_rQMC(i)] = Lookback_Option_ADGBS_RQMC(M, run_rqmc, sobol_num, epsilon(i));
end