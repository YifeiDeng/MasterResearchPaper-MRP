% code for producing Table 6.4

% 1. the depth m is chosen to be 8, sufficiently large
% 2. initialized the vector of upper and lower barrier
% 3. initialized simulation runs for both MC and rQMC 
% 4. initialized # of copies of sobol sequence with a rotation shift
M = 8;
B_u = [105, 110, 120];
B_l = [80, 85, 90];
sobol_num = 25;
run = sobol_num*2^(16);
run_rqmc = 2^(16);

% Initialized vector for storing
V_upin = zeros(size(B_u,2),1);
V_downout = zeros(size(B_u,2),1);

V_upin_rQMC = zeros(size(B_u,2),1);
V_downout_rQMC = zeros(size(B_u,2),1);

se_upin = zeros(size(B_u,2),1);
se_downout = zeros(size(B_u,2),1);

se_upin_rQMC = zeros(size(B_u,2),1);
se_downout_rQMC = zeros(size(B_u,2),1);

for i = 1:size(B_u,2)
    tic
    [V_upin(i),se_upin(i)] = Barrier_Option_ADGBS_Call_upin(M, run, B_u(i));
    toc
    
    tic
    [V_downout(i),se_downout(i)] = Barrier_Option_ADGBS_Call_downout(M, run, B_l(i));
    toc
    
    tic
    [V_upin_rQMC(i),se_upin_rQMC(i)] = Barrier_Option_ADGBS_Call_upin_RQMC(M, run, sobol_num, B_u(i));
    toc
    
    tic
    [V_downout_rQMC,se_downout_rQMC] = Barrier_Option_ADGBS_Call_downout_RQMC(M, run, sobol_num, B_l(i));
    toc
end