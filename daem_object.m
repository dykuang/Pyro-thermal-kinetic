function asum = daem_object(x, y, T, beta, k0, option, order)
% x: E1, sigma1, E2, sigma2, ..., w1, wn-1
    asum = 0;
    if option == 0
        for i =1:length(y)
            asum = asum+ (y(i)-daem_rate_func(T(i), x, beta, k0))^2;
%             asum = asum+ abs(y(i)-daem_rate_func(T(i), x, beta, k0, T0));
        end
        
    else
        for i =1:length(y)
            asum = asum+ ( y(i)-daem_alpha_func(T(i), x, beta, k0) )^2;
        end
        
        if order > 1
            asum = asum*(x(2*order+1));
            for n = 2:order
                comp_sum = 0;
                for i =1:length(y)
                    comp_sum = comp_sum + daem_alpha_func(T(i), x((2*n-1):2*n), beta, k0);
                end
                
                if n<order
                    w_ind = 2*order+n-1;
                    comp_sum = x(w_ind)*comp_sum;
                elseif n==order
                    w_ind = (2*order+1):(3*order-1);
                    wn = 1 - sum(x(w_ind));
                    comp_sum = wn*comp_sum;                    
                end
                asum = asum+comp_sum;
            end
        end
    end
    asum = asum/length(y);
end