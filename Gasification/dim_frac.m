%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses empirical regression formula for fractor dimension
% estimation of the sample surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,D,R2] = dim_frac(filename)
%% importdata
data = importdata(filename);
x = log(-log(data(20:37,1)));
y = log(data(20:37,2));

lm = fitlm(x, y);

C = lm.Coefficients.Estimate(1);
D = lm.Coefficients.Estimate(2)+3;
R2 = lm.Rsquared.Adjusted;
end
