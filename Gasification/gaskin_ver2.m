
%% loading parameters
filename = 'na.mat';
data = importdata(filename);
time = 0:0.05:600;
scale = 1e2;
gap = 10;
time = time(1:gap:end);
data=data(1:gap:end,:);
dt = 0.05*gap;
dim = size(data);
alpha = (0:0.05:1)';
T = 273.15+[600; 700; 800];
R = 8.31447;

y = zeros(dim(1), dim(2));
sample = zeros(size(alpha,1),dim(2));

%% select fit type
ft = fittype('model_integrated(x, k, m)');
% ft = fittype('model_integrated(x, k, n, m)');
% ft = fittype('model_grain(x,k)');

f = cell(dim(2),1);
k = zeros(3,1);
m = k;
t = zeros(dim(2),1);

for i = 1:dim(2)
    f0 = fit(time',data(:,i),'spline');
    y(:,i) = differentiate(f0,time);
    [x, ind] = unique(data(:,i));
    sample(:,i) = interp1(x,y(ind,i),alpha,'spline');           % sampling
    f{i} = fit(alpha,sample(:,i),ft,'StartPoint',[2, 0.4]);  % model fit
    
%     f{i} = fit(data(:,i),y(:,i),ft,'StartPoint',[2, 0.5]);
%     sample(:,i) = f{i}(alpha);
%     
    k(i) = (f{i}.k)/scale;
    m(i) = f{i}.m;
    figure(i)
    plot(data(:,i), y(:,i));
    hold on
    plot(f{i},alpha,sample(:,i));
    hold off
    
    t(i) = interp1(data(ind,i),time(ind),0.5);
end

%% regression method I
lmI = fitlm(-1000./T/R, log(k));

%% regression method II
E = zeros(size(alpha,1),1);
lnA = E;
R2 = E;
for i = 1 : size(alpha,1)
    lm = fitlm(1./T*1000, log(sample(i,:)'));
    E(i) = lm.Coefficients.Estimate(2)*(-R);
    lnA(i) = lm.Coefficients.Estimate(1);
    R2(i) = lm.Rsquared.Adjusted;
end

%% regression method III

lmIII=fitlm(1./T*1000/R, log(-log(0.5)./t));

%% regression method IV

lmIV=fitlm(1./T*1000, log(-log(0.5).^(1-m)./t./(1-m)));

%% regression method V



