function [] = daem_distribution(w, E, sigma)
    figure()
    x = linspace(50, 500, 200);
    f = cell(1, length(w));
    val_sum = 0;
    type = ["b--", "k--", "r--", 'g--', 'y--', 'c--'];
    color = ["b", "k", "r", 'g', 'y', 'c'];
    hold on
    for i =1:length(w)
        label(i) = sprintf("component %d", i);
        f{i} = @(x) w(i)./sqrt(2*pi)./sigma(i).*exp( - 0.5*(x-E(i)).^2./sigma(i).^2);
        val_i = f{i}(x);
%         plot(x, val_i, type(i));
        patch([x fliplr(x)], [val_i fliplr(zeros(size(val_i)))], color(i), 'linestyle', '--', 'FaceAlpha',.7)
        val_sum = val_sum + val_i;
    end
    label(length(w)+1) = sprintf("sum");
%     plot(x, val_sum, type(length(w)+1));
    patch([x fliplr(x)], [val_sum fliplr(zeros(size(val_sum)))], 'g', 'linestyle', '--', 'FaceAlpha',.3)
    xlabel('$ E$ (kJ/mol)', 'Interpreter','latex')
    ylabel('f(E)')
    legend(label)

end