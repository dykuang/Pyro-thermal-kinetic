function [E, Ap, Rsquare] =  E_KAS(csv_file, sample_ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function applies KAS method for deciding activation engergy
% sample_ind selects arrays for regression[
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
      lm_kas = fitlm(1./data(i,sample_ind)', log(beta(sample_ind)./(data(i,sample_ind).^2))');
      E(i,1) = -R*lm_kas.Coefficients.Estimate(2);
      SE(i,1) = R*lm_kas.Coefficients.SE(2);
      Ap(i,1) = lm_kas.Coefficients.Estimate(1);
      Rsquare(i,1) = lm_kas.Rsquared.Adjusted;
end

