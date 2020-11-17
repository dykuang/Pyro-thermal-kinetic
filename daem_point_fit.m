function [x, coef, w, R2] = daem_point_fit(r, T, k0, beta, order, ini, exclude_lb)
    [coef, components, w, gof] = gaussian_decompose(r, T, order, exclude_lb);
    
    x = zeros(order, length(T), 2);
    pred = zeros(order, length(T));
    for c = 1:order
        for i = 1:length(T)
            x(c,i,:) = daem_fit(components(i,c), T(i), ini, beta, k0, 0, 1);
            pred(c, i) = daem_rate_func(T(i), x(c,i,:), beta, k0);
        end
    end
    
    type = ["b--", "k--", "r--", 'g--', 'y--', 'c--'];
    figure()
    hold on
    for i = 1:order
        label(i) = sprintf("Component %d", i);
        plot(T, pred(i,:), type(i))
    end
    plot(T, r, 'k.')
    label(order+1) = "exp data";
    pred_sum = sum(pred, 1);
    plot(T, pred_sum, 'm-')
    label(order+2) = "DAEM fit";
    dim = [0.7 0.2 0.3 0.3];
    str = sprintf('R^2: %.4f', gof.rsquare);
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    legend(label)
    xlabel('T (K)')
    ylabel('$\frac{d\alpha}{dT}$', 'Interpreter','latex')

    Bbar = mean(r);
    SStot = sum((r - Bbar).^2);
%     SSreg = sum((pred_sum - Bbar).^2);
    SSres = sum((r - pred_sum').^2);
    R2 = 1 - SSres/SStot;
    
%     disp(SStot)
%     disp(SSres)
%     disp(R2)

end