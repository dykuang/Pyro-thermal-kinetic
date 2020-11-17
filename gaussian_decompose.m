function [coef_sorted, components, w, gof] = gaussian_decompose(y, x, order, exclude_lb)
% inputs:
%        y: 1d array
%        x: 1d grid where y are sampled
%    order: number of components to decompose
% exclude_lb:  points to disgard, 
%
% output:
%     coef_sorted: parameters for components gaussian
%      components: decomposed gaussians
%               w: weights for each gaussian
%             gof: goodness of fit
% StartPoint matters!
    g_type = sprintf('gauss%d', order);
    fo = fitoptions('Method','NonlinearLeastSquares',...
                    'Robust', 'LAR', ...
                    'TolFun', 1e-7, ...
                    'TolX', 1e-7, ...
                    'Lower',zeros(1,3*order),...
                    'Upper',[],...
                    'StartPoint',[1, 450, 20, 1, 600, 50, 1, 700, 100], ...
                    'Exclude', x > exclude_lb);
%   'coal': [],  ...  
%   'corn': [1, 450, 20, 1, 600, 40, 1, 700, 200],  ...  
%   'cc': [1, 450, 20, 1, 600, 50, 1, 700, 100],  ...  
    y_scale = max(y);
    [gfit, gof] = fit(x, y/y_scale, g_type, fo);
    coef = coeffvalues(gfit);
    coef(1:3:end) = coef(1:3:end)*y_scale;
    coef_sorted = zeros(1, 3*order);
    figure()
    type = ["b--", "k--", "r--", 'g--', 'y--', 'c--'];
    
    components = zeros(length(x), order);
    
    w = coef(1:3:end).*coef(3:3:end);
    w=w/sum(w);
    
    [w, w_ind] = sort(w, 'descend');
    hold on
    for i = 1:order
        label(i) = sprintf("peak %d", i);
        components(:, i) = gaussian(x, coef(3*w_ind(i)-2),coef(3*w_ind(i)-1),coef(3*w_ind(i)) );
        coef_sorted(3*i-2 : 3*i) = coef(3*w_ind(i)-2:3*w_ind(i));
        plot(x, components(:,i), type(i) )
    end
    plot(x, y, 'k.')
    label(order+1) = "exp data";
    plot(x, sum(components, 2), 'm-')
    label(order+2) = "fit";
    dim = [0.7 0.2 0.3 0.3];
    str = sprintf('R^2: %.4f', gof.rsquare);
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    legend(label)
    xlabel('T (K)')
    ylabel('$\frac{d\alpha}{dT}$', 'Interpreter','latex')
    disp(gof)
       
    function val = gaussian(x, a,b,c)
        val = a*exp(-((x-b)./c).^2);
    end

end