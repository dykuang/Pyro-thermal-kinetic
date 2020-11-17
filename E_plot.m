function [] = E_plot(x, E,sigma)
    color = ["b", "k", "r", 'g', 'y', 'c'];
    
    figure()
    hold on
    yu = E+sigma;
    yl = E-sigma;
    dim = size(E);
%     disp(size(yu))
    for i = 1:dim(1)
%         disp(size(yu(i,:)))
        label(i) = sprintf("component %d", i);
        p = plot(x, E(i,:), color(i), 'linewidth', 2);
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        patch([x fliplr(x)], [yu(i,:) fliplr(yl(i,:))], color(i), 'linestyle', '--', 'FaceAlpha',.3)
    end
    
    
    xlabel('$\alpha$', 'Interpreter','latex')
    ylabel('E (kJ/mol)')
    xlim([0.2 0.8])
    legend(label, 'Location','north')
    hold off
end