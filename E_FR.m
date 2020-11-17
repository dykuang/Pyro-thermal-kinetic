function[E,Ap,Rsquare] = E_FR(csv_file_T, csv_file_a, sample_ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function utilizes FR method in calculating activation energy
% csv_file_T contains alpha vs T. Each column corresponds to a fixed
% heating rate.
% csv_file_a contains alpha vs ln(beta*d alpha/dT).Each column corresponds to a fixed
% heating rate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


par = set_parameter();
R = par.R;
beta = par.beta;
alpha = par.alpha;


Tdata = importdata(csv_file_T);
% Tdata = Tdata.data;

Y = importdata(csv_file_a);
% Y = log([beta(1)*Y(:,1) beta(2)*Y(:,2) beta(3)*Y(:,3) beta(4)*Y(:,4) beta(5)*Y(:,5)]);
% Y = Y.data;
Y = log(Y);

sample_size = size(Tdata,1);
E = zeros(sample_size,1);
Rsquare = E;
Ap = E;
SE = E;

for i = 1:size(alpha,1)
          lm_fr = fitlm(1./Tdata(i,sample_ind)', Y(i,sample_ind)');
          E(i) = -R*lm_fr.Coefficients.Estimate(2);
          SE(i) = R*lm_fr.Coefficients.SE(2);
          Ap(i) = lm_fr.Coefficients.Estimate(1);
          Rsquare(i) = lm_fr.Rsquared.Adjusted;
end
