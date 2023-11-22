clc
clear

% Define los extremos del intervalo [a, b] y la cota de error
a = 0;
b = 1;
cotaerror = 1e-6;  % Por ejemplo, puedes usar 1e-6 para una tolerancia de 10^-6
% Número máximo de iteraciones
maxiter = 1000;

disp('Biseccion')
% Llama a la función de bisección y almacena los resultados en las variables de salida
[sol, niter, error] = biseccion(a, b, cotaerror, maxiter);
