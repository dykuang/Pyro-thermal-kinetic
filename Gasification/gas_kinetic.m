%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs the first step gasification kinetic analysis
% obtaining coefficient k for different model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k,m,f,gof, r_sample] = gas_kinetic(filename)

alpha = (0.01:0.025:0.99)';
datar = importdata(filename);
data = datar(:,[1 3 5]);  % tim, alpha and its derivative-smooth

%% select fit type
% ft = fittype('model_integrated(x, k, m)');
ft = fittype('model_NDM(x, k, m, r)');
% ft = fittype('model_grain(x,k)');    
% ft = fittype('model_RPM(x, k, m)');  


[x, ind] = unique(data(:,2));
r_sample = interp1(x,data(ind,3),alpha,'spline');           % sampling

scale = 60;     
exclude = (alpha>1 |alpha <0.34);
[f,gof] = fit(alpha,r_sample/scale,ft,...
              'StartPoint',[1e-3, 0.55, 0.4],...
              'Exclude',exclude,...
              'Robust','Bisquare');  % model fit -- RPM
% [f,gof] = fit(alpha,r_sample*scale,ft,'StartPoint',1); # modify initial
% condition, seems to have multiple minimum

k = (f.k);
m = f.m;
% m = 2/3;
% disp(gof)    
% figure()
% plot(data(:,2), data(:,3));
% hold on
% plot(data(:,2), rate);
% plot(f,alpha,r_sample/scale);
% hold off







end