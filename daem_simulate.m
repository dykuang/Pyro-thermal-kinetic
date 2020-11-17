function pred = daem_simulate(T, x, k0, beta, compare, y)

num_pts = length(T);
pred = zeros(1,num_pts);
order = (length(x)+1)/3;
if length(x)<3
    for i = 1:num_pts
        pred(i) = daem_alpha_func(T(i), x, beta, k0);
    end
else
    w = x(end-order+2:end);
    wn = 1 - sum(w);
    for i = 1:num_pts
        for j = 1:order-1
            pred(i) = pred(i) + w(j)*daem_alpha_func(T(i), x( (2*j-1):(2*j) ), beta, k0);
        end
        pred(i) = pred(i) + wn*daem_alpha_func(T(i), x(2*order-1:2*order), beta, k0);
    end
    
end
    
if compare
    figure()
    plot(T, y, 'b-')
    hold on
    plot(T, pred, 'k--')
    legend('Experiment', 'DAEM model')

end
