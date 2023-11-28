%% Dado que me piden error menor a 10^-3, tomo 1000 puntos por cada intervalo
%Es decir que  h = 1/1000
N = 1000;
h = 1/N;
n = linspace(h,10-h,10*N);
Ts = @(f,n,k) h/3 * (f(n(1)) + f(n(k)) +  4 * sum(f(n(2:2:k-1))) + 2 * sum(f(n(3:2:k-1))) );


f = @(t) t.^3 ./ (exp(t)-1);

phi = @(x) Ts(f,n,x*N);

%creo la tabla de valores
phis = [phi(1),phi(2),phi(3),phi(4),phi(5),phi(6),phi(7),phi(8),phi(9),phi(10)]';
errors = [integral(f,0,1) - phi(1),
    integral(f,0,2) - phi(2),
    integral(f,0,3) - phi(3),
    integral(f,0,4) - phi(4),
    integral(f,0,5) - phi(5),
    integral(f,0,6) - phi(6),
    integral(f,0,7) - phi(7),
    integral(f,0,8) - phi(8),
    integral(f,0,9) - phi(9),
    integral(f,0,10) - phi(10)];
tb = table(phis,errors);
disp('Tabla:');
disp(tb);

%% 2 
x = linspace(h,30,100000);
plot(x,f(x));
grid();
