%% Primero defino el intervalo, los nodos a usar y la funcion a integrar.
u = 10;
N = 200;
n = linspace(0,u,N);
h = n(2) - n(1); %Al ser equiespaciadas, puedo consiguer h asi
c = @(t) cos(pi/2 * t.^2);
s = @(t) sin(pi/2 * t.^2);

%% Aplico la T(f) de trapceios

Tt = @(f,n) h/2 * (f(n(1)) + f(n(N)) + 2 * sum(f(n(2:N-1))));

C = Tt(c,n)
S = Tt(s,n)
Cmatlab = integral(c,0,10);
Smatlab = integral(s,0,10);
errorC = Cmatlab - C
errorS = Smatlab - S

%% Aplico T(f) de Simpson
% Es anti-intuitivo, pero los pares son los impares en matlab por arrancar
% a contar desde el 1
Ts = @(f,n) h/3 * (f(n(1)) + f(n(N)) +  4 * sum(f(n(2:2:N-1))) + 2 * sum(f(n(3:2:N-1))) );

C = Ts(c,n)
S = Ts(s,n)

errorC = Cmatlab - C
errorS = Smatlab - S