% This script for making bar plots from the pipe line
% batchE -> thermo-dynamics
figure()
alpha = 0.2:0.1:0.8;
plot(alpha', 1e-3*E(3:2:end-3,4), 'bo-', alpha', 1e-3*E(3:2:end-3,3), 'ko-');
legend('FWO', 'VA', 'Location', 'northwest')
set(gca, 'xticklabel',{'0.2'; '0.3'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'});
xlabel('$\alpha$', 'Interpreter','latex')
grid on
ylim([min(min(1e-3*E(3:2:end-3,4))) - 10, max(max(1e-3*E(3:2:end-3,3))) + 10])
ylabel('E (kJ/mol)')
xlim([0.2 0.8])


k0  = log10(exp(m_lnk));
bar_err_plot(k0(3:2:end-3,:), 0*k0(3:2:end-3,:))
ylabel('$\log_{10} k_0 (s^{-1})$', 'Interpreter','latex')
ylim([min(min(k0(3:2:end-3,:))) - 1, max(max(k0(3:2:end-3,:))) + 1])
% saveas(a, 'path\to\file\abc1.png','png');

bar_err_plot(1e-3*m_dH(3:2:end-3,:), 0*m_dH(3:2:end-3,:))
ylabel('$\Delta H$ (kJ/mol)', 'Interpreter','latex')
% ylim([min(min(1e-3*m_dH(3:2:end-3,:))) - 2, max(max(1e-3*m_dH(3:2:end-3,:))) + 2])

bar_err_plot(1e-3*m_dG(3:2:end-3,:), 0*m_dG(3:2:end-3,:))
ylabel('$\Delta G$ (kJ/mol)', 'Interpreter','latex')
% ylim([min(min(1e-3*m_dH(3:2:end-3,:))) - 2, max(max(1e-3*m_dH(3:2:end-3,:))) + 2])


bar_err_plot(1e-3*m_dS(3:2:end-3,:), 0*m_dS(3:2:end-3,:))
ylabel('$\Delta S$ (kJ/mol)', 'Interpreter','latex')
ylim([min(min(1e-3*m_dS(3:2:end-3,:))) - 0.02, max(max(1e-3*m_dS(3:2:end-3,:))) + 0.02])

bar_err_plot(1e-3*m_dE(3:2:end-3,:), 0*m_dE(3:2:end-3,:))
ylabel('$\Delta E$ (kJ/mol)', 'Interpreter','latex')
% ylim([160, 210])

bar_err_plot(m_A(3:2:end-3,:), 0*m_dH(3:2:end-3,:))
ylabel('$A$', 'Interpreter','latex')
ylim([min(min(m_A(3:2:end-3,:))) - 0.05, max(max(m_A(3:2:end-3,:))) + 0.05])
legend('Location', 'northeast')