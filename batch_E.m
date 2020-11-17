%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This script is used to benchmark different Activation Energy algorithm%
%  file_T: A .mat file containing temperatures during reaction whose rows correspond to conversion degree and
%  columns diferent heating rate
%  filt_rate: A .mat file containing the rate of conversion degree
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_T = 'cc_T.mat';
file_rate = 'cc_r.mat';
ind = [1 2 3 4 5];
alpha =0.10:0.05:0.95;
R = 8.3144598;
beta = [5;10;15;20;25];

[E(:,1),Ap_KAS, Rsquare_KAS] = E_KAS(file_T, ind);
[E(:,2),Ap_FR, Rsquare_FR] = E_FR(file_T,file_rate,ind);
E(:,3) = E_VA(file_T, ind);
[E(:,4),Ap_FWO, Rsquare_FWO] = E_FWO(file_T,ind);
E(:,5) = E_VO(file_T,ind);
% E = E/1000;
% E(end-3:end,:) = E(end-4,3);
figure
plot(alpha', E(:,1), alpha', E(:,2), alpha', E(:,3), alpha', E(:,4), alpha', E(:,5));
legend('KAS','FR','VA','FWO','VO')



%% calculate A with model from masterplot
% n = 4; %remember to change this for different model!
% % g = @(x) (1-(1-x).^(1/3)).^2;
% % f = @(x) 1.5*(1-x).^(2/3).*(1-(1-x).^(1/3)).^(-1);
% if n>1
%     g = @(x) ((1-x).^(1-n)-1)/(n-1);
% else 
%     g = @(x) -log(1-x);
% end
% 
% f = @(x) (1-x).^n;
% 
% 
% lnA_KAS = (Ap)+log(g(alpha'))+log(E(:,1))-log(R); % intercept from KAS/(DAEM?)
% lnA_FR = Ap_FR - log(f(alpha')); % intercept from FR
% 
% % Direct from rate equation
% j = 3; % index for method
% aT = importdata(file_T);
% dadt = importdata(file_rate);
% lnA_r  = log(dadt)+[E(:,j) E(:,j) E(:,j) E(:,j) E(:,j)]./R./aT - log(f(alpha'))*ones(1,5);
% 
% % Direct from model function g
% Pint = lnA_r;
% for s = 1 : 5
% Pint(1,s)=integral(@(x) exp(-E(1,s)/R./x), 0, aT(1,s), 'AbsTol',1e-6,'RelTol',1e-3);
% for i = 2 : size(aT,1)
%     Pint(i,s) = Pint(i-1,s)+...
%     integral(@(x) exp((-E(i,s)-E(i-1,s))/R./x/2), aT(i-1,s), aT(i,s), 'AbsTol',1e-6,'RelTol',1e-3); 
% end
% end
% 
% P = 0.00484*exp(-1.0516*[E(:,3) E(:,3) E(:,3) E(:,3) E(:,3)]./R./aT);
% lnA_VA = log(g(alpha'))*ones(1,5)+ ones(18,1)*log(beta)'-log([E(:,3) E(:,3) E(:,3) E(:,3) E(:,3)]) - log(P) + log(R);
% lnA_VA_int = log(g(alpha'))*ones(1,5)+ ones(18,1)*log(beta)' - log(Pint);
% % From a minimization
% lnA_min = lnA_KAS;
% % lnA_min_s = lnA_min;
% for i = 1 : length(alpha)
%     error = @(x) 0.9*norm(log(dadt(i,3))+(E(i,j)/R/aT(i,3)) - log(f(alpha(i))) -(x))...
%      +norm(log(g(alpha(i)))+ log(beta(3)) + log(R) - log(P(i,3)) -log(E(i,j)) - (x));       
% % +norm(log(g(alpha(i)))+ log(beta(3)) - log(Pint(i,3)) - (x));
%                 
%     lnA_min(i) = fminbnd(error, 25,70);
% %     errorI(i) = norm(log(dadt(i,3))+(E(i,j)/R/aT(i,3)) - log(f(alpha(i))) -(lnA_min(i)));
% %     errorII(i) = norm(log(g(alpha(i)))+ log(beta(3)) + log(R) - log(P(i,3)) -log(E(i,j)) - (lnA_min(i)));
% 
%     % reverse solve for E
%     
% end
%% Recover dadt v.s T curve
% data_exp = preprocess('rawcoal15.csv'); % be aware of the final temperature.
% dadt_sim_KAS = exp(lnA_KAS).*exp(-E(:,1)./R./aT(:,3)).*f(alpha');
% dadt_sim_VA = exp(lnA_VA_int(:,3)).*exp(-E(:,3)./R./aT(:,3)).*f(alpha');
% dadt_sim_R = exp(lnA_r(:,3)).*exp(-E(:,j)./R./aT(:,3)).*f(alpha');
% dadt_sim_min = exp(lnA_min).*exp(-E(:,j)./R./aT(:,3)).*f(alpha');
% dadt_sim_ave = exp(0.5*(lnA_r(:,3)+lnA_VA(:,3))).*exp(-E(:,j)./R./aT(:,3)).*f(alpha');
% figure()
% plot(data_exp(:,2),data_exp(:,7),...
%     aT(:,3)-273.15,dadt_sim_R,'+',...
%     aT(:,3)-273.15,dadt_sim_VA,...
%     aT(:,3)-273.15,dadt_sim_min/2,'*');
% legend('exp','eqn_r','Pint','min')

%% Recover 1-alpha v.s. curve
% eqn_g_r = exp(lnA_min).*E(:,j)./beta(3)./R.*P(:,3);
% eqn_g_r_int = exp(lnA_VA_int(:,3))./beta(3).*Pint(:,3);
% eqn_g_KAS = exp(lnA_KAS)*R./E(:,1)./exp(Ap);
% 
% g_inv = ((n-1).*eqn_g_r + 1).^(1/(1-n)); 
% figure()
% % plot(data_exp(:,2),1-data_exp(:,6),...
% %     aT(:,3)-273.15, (1-sqrt(eqn_g_r_int)).^3,'*',...
% %     aT(:,3)-273.15, (1-sqrt(eqn_g_KAS)).^3,'+',...
% %     aT(:,3)-273.15, (1-sqrt(eqn_g_r)).^3);
% % legend('exp','eqn_int','KAS','min')
% 
% plot(data_exp(:,2),1-data_exp(:,6),...
%     aT(:,3)-273.15, g_inv,'+');
% legend('exp','simulation')

