% code for producing Figure 6.1

% 1. fixed the depth m is chosen to be 8, sufficiently large
% 2. fixed the toleranece level to 10^(-6)
% 3. initialized # of copies of sobol sequence with a rotation shift
% 4. initialized vectors of simulation runs for both MC and rQMC 
M = 8;
epsilon = 10^(-6);
sobol_num = 25;
run = sobol_num * [2^(8), 2^(9), 2^(10), 2^(11), 2^(12), 2^(13), 2^(14), 2^(15),2^(16),2^(17)];
run_rQMC = [2^(8), 2^(9), 2^(10), 2^(11), 2^(12), 2^(13), 2^(14), 2^(15),2^(16),2^(17)];

% Initialized vector for storing
V = zeros(size(run,2),1);
V_rQMC = zeros(size(run,2),1);
se = zeros(size(run,2),1);
se_rQMC = zeros(size(run,2),1);

for i = 1:size(run,2)
    [V(i),se(i)] = Lookback_Option_ADGBS(M, run(i), epsilon);
    [V_rQMC(i),se_rQMC(i)] = Lookback_Option_ADGBS_RQMC(M, run_rQMC(i), sobol_num, epsilon);
end

% Gives the plot
figure(1)
plot(run, se, 'b-')
xlim([min(run)-10^4, max(run)+10^4])
xlabel('Total Simulation Runs')
ylabel('Standard Error')
title('Monte-Carlo Method')

figure(2)
plot(run, se_rQMC, 'r-')
xlim([min(run)-10^4, max(run)+10^4])
xlabel('Total Simulation Runs')
ylabel('Standard Error')
title('Randomized Quasi-Monte-Carlo Method')


