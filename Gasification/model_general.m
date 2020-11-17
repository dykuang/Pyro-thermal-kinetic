function f = model_general(x, k, n, m)
 
  f = k*(1-x).^n.*(1-log(1-x)).^m;
  
end