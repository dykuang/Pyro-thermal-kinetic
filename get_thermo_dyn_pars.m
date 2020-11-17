function [m_dH, m_dG, m_dS, m_dE, m_A, m_lnk] = get_thermo_dyn_pars(E, T, T0, Tm, beta,k0, return_mean)
% function [dH, dG, dS, dE, A, lnk] = get_thermo_dyn_pars(E, T, T0, Tm, beta, k0)
    n = length(beta);
    par = set_parameter();
    if k0 > 0
        lnk = k0.*ones(length(E), n);
    else
        lnk = zeros(length(E), n);
    end
    dH = zeros(length(E), n);
    dG = zeros(length(E), n);
    dS = zeros(length(E), n);
    dE = zeros(length(E), n);
    A = zeros(length(E), n);
    for i = 1:n
        if k0<=0
            lnk(:,i) = log(beta(i)) + log(E) + E/par.R/Tm(i) - log(par.R) - 2*log(Tm(i));
        end
        dH(:, i) = E - par.R*T(:,i);
        dG(:, i) = E + par.R*Tm(i).*( log(par.B) + log(Tm(i)) - log(par.h) - lnk(:,i) );
%         dG(:, i) = E + par.R*T(:,i).*( log(par.B) + log(T(:,i)) - log(par.h) - lnk(:,i) );
        dS(:, i) = ( dH(:, i) - dG(:, i) )./Tm(i);
%         dS(:, i) = ( dH(:, i) - dG(:, i) )./T(:,i);
        dE(:, i) = dH(:,i) - T0(i)*dS(:, i);
        A(:, i)  = dE(:,i)./dH(:,i);
    end
    
    if return_mean
        m_lnk = mean(lnk, 2);
        m_dH = mean(dH, 2);
        m_dG = mean(dG, 2);
        m_dS = mean(dS, 2);
        m_dE = mean(dE, 2);
        m_A = mean(A, 2);
    else
        m_lnk = lnk;
        m_dH = dH;
        m_dG = dG;
        m_dS = dS;
        m_dE = dE;
        m_A = A;
    end

end