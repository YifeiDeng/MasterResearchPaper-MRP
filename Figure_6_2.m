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

% Gives the plot(MC)
figure(1)
plot(M, V, 'b-')
xlim([min(M)-0.1, max(M)+0.1])
xlabel('Depth(m)')
ylabel('Estimated Price')
title('Monte-Carlo Method')

% log-log scale
figure(2)
plot(log(M), log(V), 'b-')
xlim([min(M)-0.1, max(M)+0.1])
xlabel('log(Depth(m))')
ylabel('log(Estimated Price)')
title('Monte-Carlo Method (log-log scale)')
axis auto

% Gives the plot(rQMC)
figure(3)
plot(M, V_rQMC, 'r-')
xlim([min(M)-0.1, max(M)+0.1])
xlabel('Depth(m)')
ylabel('Estimated Price')
title('Randomized Quasi-Monte-Carlo Method')

% log-log scale
figure(4)
plot(log(M), log(V_rQMC), 'r-')
xlim([min(M)-0.1, max(M)+0.1])
xlabel('log(Depth(m))')
ylabel('log(Estimated Price)')
title('Randomized Quasi-Monte-Carlo Method (log-log scale)')
axis auto

% Conbined the plot
coef_MC = polyfit(log(M), log(V)',1);
coef_rQMC = polyfit(log(M), log(V_rQMC)',1);
hold on
plot(log(M), log(V), 'b-')
plot(log(M), log(M)*coef_MC(1)+coef_MC(2), 'b--')
plot(log(M), log(V_rQMC), 'r-')
plot(log(M), log(M)*coef_rQMC(1)+coef_rQMC(2), 'r--')
hold off
xlabel('log(Depth(m))')
ylabel('log(Estimated Price)')
title('MC and rQMC (log-log scale)')
axis auto
legend('MC','Linear Fit for MC: y = -0.0352x+2.3045','rQMC','Linear Fit for rQMC: y = -0.0220x+2.2826','FontSize',13)
