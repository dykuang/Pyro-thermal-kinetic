function [x,fval,exitflag,output] = daem_fit(y, T, x0, beta, k0, option, order)
    lb = zeros(1,3*order-1);
    ub = ones(1, 3*order-1);
    for i = 1:order
        lb(2*i-1) = 50;
        lb(2*i) = 0.001;
        ub(2*i-1) = 350;
        ub(2*i) = 50;
    end
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    options = optimoptions('patternsearch','MaxIterations', 1e6, ...
                           'FunctionTolerance', 1e-8, 'MeshTolerance',1e-8);
    if option == 0
        [x,fval,exitflag,output] = patternsearch(@(x) ...
                            daem_object(x, y, T, beta, k0, option) , ...
                                   x0,A,b,Aeq,beq,lb,ub, nonlcon, options);
    else
        [x,fval,exitflag,output] = patternsearch(@(x) ...
                            daem_object(x, y, T, beta, k0, option,order) , ...
                                   x0, A,b,Aeq,beq,lb,ub, nonlcon, options);
    end


end