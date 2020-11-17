%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculate the kinectic parameters in gasification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[E,lnA,R2,k,m,gof,r_sample, fd] = kin_per_gas(ind)
% file1 = 'coal750.mat';
% file2 = 'coal800.mat';
% file3 = 'coal850.mat';
R = 8.31447;
alpha = (0.01:0.025:0.99)';

T = [700; 750; 800; 850];
T = T(ind);
N = length(ind);

k = zeros(N,1);
m= k;
gof = cell(N,1);
f = gof;
r_sample = zeros(length(alpha), N);

figure()
hold on
for i =1 : N
    file = sprintf('na2fe%d.mat',int16(T(i)));
    [k(i),m(i),f{i},gof{i},r_sample(:,i)] = gas_kinetic(file);
    plot(alpha, r_sample(:,i)/60,'.', 'DisplayName', sprintf('coal%d-data',int16(T(i))));
    plot(alpha, f{i}(alpha),'-', 'DisplayName', sprintf('coal%d-fit',int16(T(i))));
    disp(f{i}.r)
end

lm = fitlm(1./(T+273.15),log(k));
E = lm.Coefficients.Estimate(2)*(-R);
lnA = lm.Coefficients.Estimate(1);
R2 = lm.Rsquared.Adjusted;

hold off
figure()
plot(1./T, log(k),'^')

%% for extracting data purpose
r_sample = r_sample/60;
fd = r_sample;
for i = 1: N
    fd(:,i) = f{i}(alpha);
end
