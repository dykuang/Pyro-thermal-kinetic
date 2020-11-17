function alpha = daem_alpha_func(T, x0, beta, k0)
    format long
    E0 = x0(1);
    sigma0 = x0(2);
    R = 8.314;
    
    integrand = @(E) 1./sqrt(2*pi)./sigma0.*...
                     exp( (-k0.*R.*T.^2./beta./E).*exp(-E./R./T) ...
                          -0.5*(E-E0).^2./(sigma0.*2) );
                 
%     alpha = integral(integrand, 0, inf, 'ArrayValued', true);
    alpha = integral(integrand, 5e4, 3.5e5, 'ArrayValued', true);   
        
   
end