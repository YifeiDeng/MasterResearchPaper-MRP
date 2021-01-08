% code for producing Table 6.2

% 1. the depth m is chosen to be 10, sufficiently large
% 2. initialized vector of toleranece
% 3. simulation runs sets to 10^6, i.e 1 millon
M = 10;
epsilon = [10^(-2), 10^(-6), 10^(-10), 10^(-14)];
run = 10^6;

% Initialized vector for storing
point_ADGBS = zeros(size(epsilon,2),run);

for i = 1:size(epsilon,2)
    [X_, x_l, x_u, point_ADGBS(i,:),index] = Adaptive_DGBS(M, run, epsilon(i));
end

Aver_point_path = mean(point_ADGBS,2);