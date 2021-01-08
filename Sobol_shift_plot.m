% Illustration of Cranley-Patterson Rotation on Sobol Sequence with shift
clear;
point = 2^6;
qmc_set_orignial = net(sobolset(2*point),point);
qmc_set_orignial_x = qmc_set_orignial(:,1:point);
qmc_set_orignial_y = qmc_set_orignial(:,(point+1):(2*point));
r_direct_x = unifrnd(0,1,[1,point]);
r_direct_y = unifrnd(0,1,[1,point]);
shifted_qmc_x = mod(qmc_set_orignial_x+repelem(r_direct_x,point,1),1);
shifted_qmc_y = mod(qmc_set_orignial_y+repelem(r_direct_y,point,1),1);

hold on
plot(qmc_set_orignial_x(:,1),qmc_set_orignial_y(:,1),'bo')
plot(shifted_qmc_x(:,1),shifted_qmc_y(:,1),'ro')
hold off
legend('Sobol Sequence','Shifted Sobol Sequence')