
%% loading parameters
filename = 'coal.mat';
data = importdata(filename);
time = 0:0.05:600;
scale = 1e2;
gap = 10;
% precision = 1;
time = time(1:gap:end);
data=data(1:gap:end,:);
dt = 0.05*gap;
dim = size(data);
alpha = (0:0.05:0.99)';   % lesson learned! carefully debug
T = 273.15+[600; 700; 800];
R = 8.31447;


% data = floor(data*precision)/precision;
%% calculate first derivative
rate = zeros(dim(1),dim(2));
rate(2:end-1,:) = (data(3:end,:)-data(1:end-2,:))/2/dt; % central difference
% rate(2:end,:) = (data(2:end,:)-data(1:end-1,:))/dt;
rate(1,:) = (data(2,:)-data(1,:))/dt;
rate(end,:) = (data(end,:)-data(end-1,:))/dt;

y = zeros(dim(1), dim(2));
sample = zeros(size(alpha,1),dim(2));

%% select fit type
% ft = fittype('model_integrated(x, k, m)');
    ft = fittype('model_general(x, k, n, m)');
% ft = fittype('model_grain(x,k)');
    
% ft = fittype('model_RPM(x, k, m)');

f = cell(dim(2),1);
k = zeros(3,1);
m = k;
t = zeros(dim(2),1);

    
for i = 1:dim(2)
    data(data(:,i)>=1,i) = 1;
    y(:,i) = smooth(data(:,i),rate(:,i)*scale,0.02,'rloess','Robust'); % smoothing
    [x, ind] = unique(data(:,i));
    sample(:,i) = interp1(x,y(ind,i),alpha,'spline');           % sampling
%     f{i} = fit(alpha,sample(:,i),ft,'Lower',[1e-6, 2.0],...
%     'StartPoint',[1, 2.2]);  % model fit -- RPM

%     f{i} = fit(alpha,sample(:,i),ft,...
%     'StartPoint',[2, 0.5]);      % model fit -- integrated
        
     f{i} = fit(alpha,sample(:,i),ft,...
    'StartPoint',[1, 2.2, 2.2]);  % model fit -- RPM

%     f{i} = fit(data(:,i),y(:,i),ft,'StartPoint',[2, 0.5]);
%     sample(:,i) = f{i}(alpha);
%     
    k(i) = (f{i}.k)/scale;
    m(i) = f{i}.m;
    
    figure(i)
    plot(data(:,i), y(:,i));
    hold on
    plot(data(:,i), rate(:,i)*scale);
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



