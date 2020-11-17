function f = model_RPM(x, k, m)
 
  f = k*(1-x).*sqrt(1-m*log(1-x));
  
end