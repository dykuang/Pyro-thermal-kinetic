%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function preprocess the gasfication data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data,t0,tf] = preprocess(filename)
% filename = 'FE750.mat';  % two cols: first col is time, the other is weight
data_r = importdata(filename);
% tf = data_r(end,1);
t0 = data_r(end,1)-600; % isothermal in the last 600 min, important
% t0 = 145; % NA750
% t0 = 148.6; % NA800
% t0 = 144.9; % NA850
% t0 = 132;
% tf = 737.8 ;
data = data_r(data_r(:,1)>t0,:);
data = data(1:10:end,:);
data(:,1)=data(:,1)-data(1,1);
data = [data zeros(size(data,1),3)]; %expand to contain other info: alpha, derivative and its smooth
w0 = data(1,2);
wf = data(end,2);

%% calculate degree of conversion
% data(:,2)=smooth(data(:,1),data(:,2),0.10,'rloess');
data(:,3) = (w0-data(:,2))./(w0-wf);
% data(:,3)=smooth(data(:,1),data(:,3),0.15,'rloess');

%% calculate derivative and smooth it
data(1,4) = (data(2,3)-data(1,3))/(data(2,1)-data(1,1));
data(end,4)= (data(end,3)-data(end-1,3))/(data(end,1)-data(end-1,1));
data(2:end-1,4) = (data(3:end,3)-data(1:end-2,3))./(data(3:end,1)-data(1:end-2,1));

data(:,5)=smooth(data(:, 1),data(:,4),0.07); % smooth with respect to time? alpha? turns out the same, choose the appropriate level
% data(:,5) = data(:,4);
figure()
plot(data(:,1),data(:,5));
figure()
plot(data(:,3),data(:,5));
hold on
plot(data(:,3),data(:,4));
figure()
plot(data(:,1),data(:,3));
end