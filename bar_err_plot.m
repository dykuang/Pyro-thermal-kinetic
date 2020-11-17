function [] = bar_err_plot(y, err)
figure(); 
hb = bar(y); % get the bar handles
hold on;
l = cell(1,size(y,2));
colors = {[0 0 1], [0.3,0.3,0.3], [1,0,0], [0,1,0], [0.2,0.2,0.6]};
for k = 1:size(y,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,k), err(:,k), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
%     l{k} = sprintf('Component %d', k);
    l{k} = sprintf('%d K/min', k*5);
    hb(k).FaceColor = colors{k};
end
grid on
   
legend(hb,l, 'Location', 'northwest');

% Set Axis properties
set(gca, 'xticklabel',{'0.2'; '0.3'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'});
xlabel('$\alpha$', 'Interpreter','latex')
ylim([min(min(y)) - 2, max(max(y)) + 2])



end