function f = objfunction(x)

f1 = x(1) + x(2) -2*x(1)^2 - x(2)^2 + x(1)*x(2);
f = 1/(1 + f1) ;                              %% to change it to maximization