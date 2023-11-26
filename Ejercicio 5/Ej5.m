%% Primero expresamos la ecuacion que rige el sistema 
Co = 5;
Ce = 12;
c = @(t) Ce.*(1-exp(-0.04*t)) + Co.*exp(-0.06*t);
cReq = 0.85*Ce;
error = 1e-6;

%% A) 
%Para esto solo es requerito plotear f(t) = c(t) - cReq.
f = @(t)c(t)-cReq ;
x = linspace(0,85,1000);
plot(x,f(x))
%Notamos que la raiz se encuentra entre 42 y 43.

%% B)
%Propongo como g(t) = t - f(t)
g = @(t) t - f(t);
x = linspace(42,44,1000);
plot(x,g(x))

%Al plotear la función g(t), notamos que cumple las tres hipotesis
%necesarias para realizar el método de punto fijo.
figure
syms n
dg = diff(g(n));
plot(x,subs(dg,n,x));

%Finalmente hacemos aplicamos el método
xi =43;

while(1)
    xi = g(xi);
    if (abs(f(xi)) < error)
        break
    end
end
display('El valor de la raiz es: ')
display(xi)
display('El error es:')
display(f(xi))
%% C
% Ahora hay que redefinir g(t) = t - g(t)/g'(t)
syms n
df = diff(f(n));
df = @(t) subs(df,n,t)
g = @(t) t - f(t)./df(t)

xi =43;

while(1)
    xi = g(xi)
    if (abs(f(xi)) < error)
        break
    end
end
xi = vpa(xi)
e = f(xi)

display('El valor de la raiz es: ')
display(xi)
display('El error es:')
display(e)