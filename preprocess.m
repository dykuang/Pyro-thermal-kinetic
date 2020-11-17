%%%%% This function does a preprocess of the original experiment data
%%%%% input: time, temperature, weight
%%%%% output: time, temperature, weight, weight loss, dw/dt, x and dx/dt.
%%%%% Be careful to change ad value for different materials
function [data,aT,ar, Tm] = preprocess(file_name, begin, final, ad)
data_original = importdata(file_name);
% data_original = file_name;
% begin = 110;
% final = 900;

data = data_original(data_original(:,2)>=begin & data_original(:,2)<=final, :); % start at 100C
% data = data.data;
% data = data(:,1:3);
% data = biomass5;

%% determine the end temperature of analysis
% hr = diff(data(:,2))./diff(data(:,1));
% for final = 1:size(hr,1)
%     if hr(final) < 25-2
%         break
%     end
% end
% data = data(1:final,:);

%% compute weight loss, ad value, coal:8.24%  pine:0.53% corn:5.92%.
% ad = 0.0824;
w0 = data(1,3);
wf = data(end,3);
% wf = w0*ad;

% data(:,3) = smooth(data(:,1),data(:,3),0.010,'rloess');

data(:,4) = (data(:,3)-w0*ad)/(w0*(1-ad));  %need to subtract ash content

%% compute first derivative of weight loss(DTG)
% using central difference
data(1,5)=(data(2,3)-data(1,3))/(data(2,1)-data(1,1));
data(end,5)= (data(end,3)-data(end-1,3))/(data(end,1)-data(end-1,1)); 
data(2:end-1,5) = (data(3:end,3)-data(1:end-2,3))./(data(3:end,1)-data(1:end-2,1));
data(:,5)= -data(:,5)/w0*(1-ad);

% scale=10e3;
DTG_smooth = smooth(data(:,1),data(:,5),0.015,'rloess');
% f = fit(data(:,1),data(:,5),'smoothingspline');
% for i= 1 : 5
%  DTG_smooth = smooth(data(:,1),DTG_smooth,0.015,'rloess');
% end
% DTG_smooth = f(data(:,1));
%% compute conversion degree
% ash = w0*ad;
data(:,6) = (w0-data(:,3))/(w0-wf);

%% compute rate of conversion degree with respect to time/T
% using central difference.
% data(1,7)=(data(2,6)-data(1,6))/(data(2,1)-data(1,1));
data(1,7)=(data(2,6)-data(1,6))/(data(2,2)-data(1,2));
data(end,7)= (data(end,6)-data(end-1,6))/(data(end,2)-data(end-1,2)); 
data(2:end-1,7) = (data(3:end,6)-data(1:end-2,6))./(data(3:end,2)-data(1:end-2,2));

DC_smooth = smooth(data(:,2),data(:,7),0.015,'rloess');
% Itarative smoothing. It Does not seem to improve much.
% for i= 1: 5
%  DC_smooth = smooth(data(:,1),DC_smooth,0.015,'rloess');
% end

%% visualization
figure(1)
plot(data(:,2), data(:,4))  % weight and weight loss
title('weight loss')

figure(2)
plot(data(:,2), data(:,5), data(:,2), DTG_smooth) % weight loss rate
title('DTG')

figure(3)
plot(data(:,2), data(:,6)) % conversion degree
title('conversion degree')

figure(4)
plot(data(:,2), data(:,7), data(:,2), DC_smooth) % rate of conversion degree
title('rate of conversion degree')

figure(5)
plot(data(:,6), data(:,7), data(:,6), DC_smooth)
title('DC V.S Alpha')

figure(6)
plot(data(:,1),data(:,2))
title('analysis window')

%% for convenient copy to excel
data(:,5) = DTG_smooth;
data(:,7) = DC_smooth;

%% obtaining sampled points for kinetic analysis
alpha = (0.1:0.05:0.95)';
% alpha = linspace(0.01, 0.99, 101)';
% [x, ind] = unique(data(:,6));
% raw_rate = interp1(data(ind,6), DC_smooth(ind), alpha);
% raw_T = 273.15+interp1(data(ind,6), data(ind,2), alpha);
f = fit(data(:,6),data(:,2),'smoothingspline');
aT = f(alpha)+273.15;
fr = fit(data(:,6),data(:,7),'smoothingspline');
ar = fr(alpha);

%% Obtaining peak temperature
[val, ind] = max(data(:,5));
Tm = data(ind, 2) + 273.15;

end