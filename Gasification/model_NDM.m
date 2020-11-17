function fval = model_NDM(x, k, m, r)

     fval = k*exp(-(x-m).^2/r.^2/2);

end