%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function tries to fit the gasification to three models: grain,
%%% integrated and RPM. data is two column: time alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% loading parameters
filename = 'FE700.mat';
data = importdata(filename);
time = data(:,1);
% scale = 1e2;
% gap = 10;
% precision = 1;
% time = time(1:gap:end);
% data=data(1:gap:end,:);
% dt = 0.05;
dim = size(data);
alpha = (0:0.05:0.99)';   % lesson learned! carefully debug
T = 273.15+700;
R = 8.31447;


% data = floor(data*precision)/precision;
%% calculate first derivative
rate = zeros(dim(1), 1);
for i = 2 : dim(1)-1
    rate(i) = (data(i+1,2)-data(i-1,2))/(data(i+1,1)-data(i-1,1)); % central difference
end
% rate(2:end,:) = (data(2:end,:)-data(1:end-1,:))/dt;
rate(1) = (data(2,2)-data(1,2))/(data(2,1)-data(1,1));
rate(end) = (data(end,2)-data(end-1,2))/(data(end,1)-data(end,1));


%% select fit type
ft = fittype('model_integrated(x, k, m)');
%     ft = fittype('model_general(x, k, n, m)');
% ft = fittype('model_grain(x,k)');
    
% ft = fittype('model_RPM(x, k, m)');    

r_smooth = smooth(data(:,1),rate,0.1,'rloess','Robust'); % smoothing
[x, ind] = unique(data(:,2));
r_sample = interp1(x,r_smooth(ind,:),alpha,'spline');           % sampling
%     f{i} = fit(alpha,sample(:,i),ft,'Lower',[1e-6, 2.0],...
%     'StartPoint',[1, 2.2]);  % model fit -- RPM

%     f{i} = fit(alpha,sample(:,i),ft,...
%     'StartPoint',[2, 0.5]);      % model fit -- integrated
        
[f,gof] = fit(alpha,r_sample,ft,'StartPoint',[1, 2.2]);  % model fit -- RPM


k = (f.k);
m = f.m;
    
figure()
plot(data(:,2), r_smooth);
hold on
% plot(data(:,2), rate);
plot(f,alpha,r_sample);
hold off




