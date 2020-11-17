function [x, w, E_alpha_val, E_sigma_val, dH, dG, dS, dE, A] ...
         = analysis_plots(name, beta_ind, k0, ini, exclude_lb)
% 
% name = 'coal';
% beta_ind = 3;
% k0 = 1.0e14;
% ini = [200, 20];
% exculde_lb = 1100, 900 , 1100 [coal, corn, cc]

step = [30, 25, 20, 15, 10];
order = 3;
beta = [5 10, 15, 20, 25];
beta = beta(beta_ind);
step = step(beta_ind);
demo_ind = 3;

% coef = zeros(length(beta), 3*order);
w = zeros(length(beta), order);
R2 = zeros(length(beta), order);
x = cell(length(beta),1);
alpha = 0.2:0.05:0.8;
E_alpha_val = zeros(length(beta), order, length(alpha));
E_sigma_val = zeros(length(beta), order, length(alpha));

dH = zeros(length(beta), order, length(alpha));
dG = zeros(length(beta), order, length(alpha));
dS = zeros(length(beta), order, length(alpha));
dE = zeros(length(beta), order, length(alpha));
A = zeros(length(beta), order, length(alpha));

T0 = importdata(strcat(name, '_T0.mat'));
Tm = importdata(strcat(name, '_Tm.mat'));

T0 = T0(beta_ind);
Tm = Tm(beta_ind);



for b = 1: length(beta)
   data_name = strcat(name, '_', sprintf('%d', beta(b)), '_', 'tga.mat');
   data = importdata(data_name);
   B_label(b) = sprintf("%d K/min", beta(b));
   %% DG & DTG plots
   figure(1)
   hold on
   plot(data(:,2)+275.13, data(:,6))

   figure(2)
   hold on
   plot(data(:,2)+275.13, data(:,7))

   %% Multi-Gaussian DAEM fit
%    [coef, components, w, gof] = gaussian_decompose(data(1:step:end,7), data(1:step:end,2)+275.13, order)
   [x_b, coef_b, w_b, R2_b] = daem_point_fit(data(1:step(b):end,7), data(1:step(b):end,2)+275.13, ...
                                             k0, beta(b), order, ini, exclude_lb);
%    coef(b, :) = coef_b;
   x{b} = x_b;
   w(b,:) = w_b;
   R2(b) = R2_b;
   
   %% Thermo dynamic parameters
   for o = 1:order
       E_alpha_func = fit(data(1:step(b):end,6), x{b}(o,:,1)', 'smoothingspline');
       E_alpha_val(b,o,:) = E_alpha_func(alpha);  % Beta, components, pts selected
       E_sigma_func = fit(data(1:step(b):end,6), x{b}(o,:,2)', 'smoothingspline');
       E_sigma_val(b,o,:) = E_sigma_func(alpha);
       T_func = fit(data(1:step(b):end,6), data(1:step(b):end,2)+275.13, 'smoothingspline');
       T_val = T_func(alpha);
%        disp(size(T_val))
       [m_dH, m_dG, m_dS, m_dE, m_A] = get_thermo_dyn_pars( 1e3*squeeze(E_alpha_val(b,o,:)),...
                                                            T_val',...
                                                            T0(b), Tm(b), beta(b), -1, false);
                                                        
       dH(b,o,:) = m_dH;
       dG(b,o,:) = m_dG;
       dS(b,o,:) = m_dS;
       dE(b,o,:) = m_dE;
       A(b,o,:) = m_A;
   end
 
end
    
figure(1)
ylim([0 1.1])
legend(B_label, 'Location', 'southeast')
xlabel('T (K)')
ylabel('$\alpha$', 'Interpreter','latex')
hold off

figure(2)
legend(B_label, 'Location', 'northeast')
xlabel('T (K)')
ylabel('$\frac{d\alpha}{dT}$', 'Interpreter','latex')

hold off

% for i = 2:2*length(beta)+2
%     figure(i)
%     xlabel('T (K)')
%     ylabel('$\frac{d\alpha}{dT}$', 'Interpreter','latex')
% end


%% plot E
E_plot(alpha, squeeze(E_alpha_val(demo_ind,:,:)), squeeze(E_sigma_val(demo_ind,:,:)));
bar_err_plot(squeeze(E_alpha_val(demo_ind,:,1:2:end))', squeeze(E_sigma_val(demo_ind,:,1:2:end))') % bar plot
ylabel('E (kJ/mol)')
%% plot dH
% type = ["b--", "k--", "r--", 'g--', 'y--', 'c--'];
% figure()
% hold on
% for i = 1:order
%     label(i) = sprintf("component %d", i);
%     plot(alpha, 1e-3*squeeze(dH(demo_ind,i,:)), type(i))
% end
% xlabel('$\alpha$', 'Interpreter','latex')
% ylabel('$\Delta H$ (kJ/mol)', 'Interpreter','latex')
% xlim([0.2 0.8])
% legend(label, 'Location','north')
% hold off
bar_err_plot(1e-3*squeeze(dH(demo_ind,:,1:2:end))', 1e-3*squeeze(E_sigma_val(demo_ind,:,1:2:end))')
ylabel('$\Delta H$ (kJ/mol)', 'Interpreter','latex')
%% plot dG
% figure()
% hold on
% for i = 1:order
%     plot(alpha, 1e-3*squeeze(dG(demo_ind,i,:)), type(i))
% end
% xlabel('$\alpha$', 'Interpreter','latex')
% ylabel('$\Delta G$ (kJ/mol)', 'Interpreter','latex')
% xlim([0.2 0.8])
% legend(label, 'Location','north')
% hold off
bar_err_plot(1e-3*squeeze(dG(demo_ind,:,1:2:end))', 1e-3*squeeze(E_sigma_val(demo_ind,:,1:2:end))')
ylabel('$\Delta G$ (kJ/mol)', 'Interpreter','latex')
%% plot dS
% figure()
% hold on
% for i = 1:order
%     plot(alpha, 1e-3*squeeze(dS(demo_ind,i,:)), type(i))
% end
% xlabel('$\alpha$', 'Interpreter','latex')
% ylabel('$\Delta S$ (J/mol)', 'Interpreter','latex')
% xlim([0.2 0.8])
% legend(label, 'Location','north')
% hold off
bar_err_plot(1e-3*squeeze(dS(demo_ind,:,1:2:end))', 0*squeeze(dS(demo_ind,:,1:2:end))')
ylabel('$\Delta S$ (kJ/mol)', 'Interpreter','latex')
%% plot dE
% figure()
% hold on
% for i = 1:order
%     plot(alpha, 1e-3*squeeze(dE(demo_ind,i,:)), type(i))
% end
% xlabel('$\alpha$', 'Interpreter','latex')
% ylabel('$\Delta E$ (kJ/mol)', 'Interpreter','latex')
% xlim([0.2 0.8])
% legend(label, 'Location','north')
% hold off
bar_err_plot(1e-3*squeeze(dE(demo_ind,:,1:2:end))', 1e-3*squeeze(E_sigma_val(demo_ind,:,1:2:end))')
ylabel('$\Delta E$ (kJ/mol)', 'Interpreter','latex')
ylim([100 200])
%% plot A
% figure()
% hold on
% for i = 1:order
%     plot(alpha, squeeze(A(demo_ind,i,:)), type(i))
% end
% xlabel('$\alpha$', 'Interpreter','latex')
% ylabel('$A$ ', 'Interpreter','latex')
% xlim([0.2 0.8])
% ylim([0.6 1.2])
% legend(label, 'Location', 'north')
% hold off
A_m = 1 - T0(demo_ind).*dS.*(1./dH + 1e6*E_sigma_val.^2./dH.^3);
A_err = (T0(demo_ind).*dS).*1e3*E_sigma_val./dH.^2;
% bar_err_plot(squeeze(A(demo_ind,:,1:2:end))', 0*squeeze(A(demo_ind,:,1:2:end))')
bar_err_plot(squeeze(A_m(1,:,1:2:end))', squeeze(A_err(1,:,1:2:end))')
ylabel('$A$ ', 'Interpreter','latex')
ylim([0.5, 1.5])
%% f(E)
daem_distribution( w(demo_ind,:), squeeze(E_alpha_val(demo_ind,:,7)), squeeze(E_sigma_val(demo_ind,:,7)) );
end