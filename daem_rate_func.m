function rate = daem_rate_func(T, x0, beta,k0)
    format long
    E0 = x0(1);
    sigma0 = x0(2);
    R = 8.314;
    
%     integrand_1 = @(x, p) exp(-p./R./x);
    
%     inter_inte = @(E) integral( @(TT)integrand_1(TT, E), T0, T, 'ArrayValued',true);
    inter_inte = @(E) (R.*T.^2./E*1e-3).*exp(-1e3*E./R./T);
    
    integrand_2 = @(E) k0/beta/sqrt(2*pi)./sigma0*1e-3 .*...
                     exp(-1e3*E./R./T ...
                         - k0/beta.*inter_inte(E) ...
                         - 0.5*(E-E0).^2./(sigma0.^2) ...                     
                         );
                 
%     rate = integral(integrand_2, 0, inf, 'ArrayValued', true);
    rate = 1e3*integral(integrand_2, 5e1, 3.5e2, 'ArrayValued', true);
end