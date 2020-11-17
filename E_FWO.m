function [E, Ap, Rsquare, SE] =  E_FWO(csv_file, sample_ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function applies FWO method for deciding activation engergy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


par = set_parameter();
R = par.R;
beta = par.beta;
alpha = par.alpha;


data = importdata(csv_file);
% data = csv_file;
% data = data.data;

sample_size = size(data,1);
E = zeros(sample_size,1);
Rsquare = E;
Ap = E;
SE = E;


for i = 1:size(alpha,1)
      lm_fwo = fitlm(1./data(i,sample_ind)', log(beta(sample_ind)'));
      E(i) = -R*lm_fwo.Coefficients.Estimate(2)/1.052;
      SE(i) = R*lm_fwo.Coefficients.SE(2)/1.052;
      Ap(i) = lm_fwo.Coefficients.Estimate(1);
      Rsquare(i,2) = lm_fwo.Rsquared.Adjusted;
end