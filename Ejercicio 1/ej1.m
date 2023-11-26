clc

% Ejercicio 1
A = [3 -1 0 0 0; 0 5 0 -1 2; 1 -1 9 0 0; -2 3 -1 12 1; 0 0 -1 1 15];
B = [-1; 0; 1; 1; 0]; 

% a) Realizar un programa para resolverlo por eliminacion de Gauss y Sustitucion Regresiva
fprintf("Ejercicio a\n");
[AB, Y] = eliminacion_gauss(A, B);
fprintf("Matriz Triangular: \n");
disp(AB(:, 1:5))
fprintf("\nSolución al Sistema: \n");
disp(Y)

% b) Realizar un programa para resolverlo por descomposición LU y Sustitucion Regresiva
fprintf("\nEjercicio b\n")
[L, U, Y] = LU(A, B);
fprintf("Matriz L: \n");
disp(L);
fprintf("\nMatriz U: \n");
disp(U);
fprintf("\nSolución al Sistema: \n");
disp(Y);

% c) Realizar un programa para resolverlo por factorizacion QR y Sustitucion Regresiva y Progresiva según corresponda
fprintf("\nEjercicio c\n")
[Q, R, Y] = QR(A, B);
fprintf("Matriz Q: \n");
disp(Q);
fprintf("\nMatriz R: \n");
disp(R);
fprintf("\nSolución al Sistema: \n");
disp(Y);

% d) Realizar un programa para resolverlo por el metodo de Jacobi con un error menor que 10^(-6) 
fprintf("\nEjercicio d\n")
x0 = zeros(size(B));
cotaerror = 1e-6;
[Y, n, error] = metodo_jacobi(A, B, x0, cotaerror);
fprintf("Cantidad de Iteraciones: \n");
disp(n);
fprintf("\nError: \n");
disp(error);
fprintf("\nSolución al Sistema: \n");
disp(Y);

% e) Realizar un programa para resolverlo por el metodo de Gauss-Seidel con un error menor que 10^(-6) 
fprintf("\nEjercicio e\n")
x0 = zeros(size(B));
cotaerror = 1e-6;
[Y, n, error] = metodo_gauss_seidel(A, B, x0, cotaerror);
fprintf("Cantidad de Iteraciones: \n");
disp(n);
fprintf("\nError: \n");
disp(error);
fprintf("\nSolución al Sistema: \n");
disp(Y);

% % f) Realizar un programa para resolverlo por el metodo de SOR con un error menor que 10^(-6) 
fprintf("\nEjercicio f\n")
x0 = zeros(size(B));
tol = 1e-6;
maxIter = 100;
fprintf("Omega = 0.9 \n")
omega = 0.9;
[Y, iter, error] = sor(A, B, x0, omega, tol, maxIter);
fprintf("Cantidad de Iteraciones: \n");
disp(iter);
fprintf("\nError: \n");
disp(error);
fprintf("\nSolución al Sistema: \n");
disp(Y);

% x0 = zeros(size(B));
% tol = 1e-30;
% maxIter = 25;
% 
% % Inicializar vectores para almacenar resultados
% omegas = 0.1:0.1:1.9;
% errors = zeros(size(omegas));
% 
% % Iterar sobre los valores de omega
% for i = 1:length(omegas)
%     omega = omegas(i);
%     [x, ~, error] = sor(A, b, x0, omega, tol, maxIter);
%     errors(i) = error;
% end
% 
% % Graficar el error en función de omega
% figure;
% plot(omegas, errors, 'o-', 'LineWidth', 2);
% title('Error vs. Omega en el método SOR');
% xlabel('Omega');
% ylabel('Error');
% set(gca, 'YScale', 'log');
% grid on;

function x = sust_reg_sup(T,B)
    s = size(T);
    nfilas = s(1);
    ncol = s(2);
    x(nfilas) = B(nfilas)/T(nfilas,ncol);
    for i=(nfilas-1):-1:1
        for j=(i+1):ncol
            x(i) = (B(i)-dot(T(i,(i+1):ncol),x((i+1):nfilas)))/T(i,i);
        end
    end
    x=x';
end

function x = sust_reg_inf(T, B)
    s = size(T);
    nfilas = s(1);
    
    x = zeros(nfilas, 1); 

    x(1) = B(1) / T(1, 1);

    for i = 2:nfilas
        suma = B(i);
        for j = 1:(i-1)
            suma = suma - T(i, j) * x(j);
        end
        x(i) = suma / T(i, i);
    end
end

function [MA,sol]=eliminacion_gauss(A,B)
    MA=[A,B];
    s=size(MA);
    nfilas=s(1);
    ncolumnas=s(2);
    columna=1;
    for fija=1:(nfilas-1)
        for fila=(fija+1):nfilas
            MA(fila,:)=MA(fila,:)-(MA(fila,columna)/MA(fija,fija))*MA(fija,:);
        end
        columna=columna+1;
    end
    sol=sust_reg_sup(MA(1:nfilas,1:ncolumnas-1),MA(:,ncolumnas));
end

function [L, U, Y] = LU(A, B)
    [m, n] = size(A);
    
    % Inicializar matrices L, U
    L = eye(m);
    U = A;
    P = eye(m);
    
    for k = 1:min(m-1, n)
        % Buscar el índice de fila con el máximo valor en la columna k
        [~, pivot_row] = max(abs(U(k:m, k)));
        pivot_row = pivot_row + k - 1;

        % Intercambiar filas en U y P
        U([k, pivot_row], :) = U([pivot_row, k], :);
        P([k, pivot_row], :) = P([pivot_row, k], :);

        % Eliminación gaussiana para calcular L y U
        for i = k+1:m
            factor = U(i, k) / U(k, k);
            L(i, k) = factor;
            U(i, :) = U(i, :) - factor * U(k, :);
        end
    end
    
    LB=[L,B];
    s=size(LB);
    nfilas=s(1);
    ncolumnas=s(2);
    Y_aux=sust_reg_inf(LB(1:nfilas,1:ncolumnas-1),LB(:,ncolumnas));
    
    ULB=[U,Y_aux];
    s=size(ULB);
    nfilas=s(1);
    ncolumnas=s(2);
    Y=sust_reg_sup(ULB(1:nfilas,1:ncolumnas-1),ULB(:,ncolumnas));
end

function [Q, R, Y] = QR(A, B)
    [m, n] = size(A);
    Q = zeros(m, n);
    R = zeros(n);

    for j = 1:n
        % Proyección ortogonal sobre los vectores ya calculados
        v = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end

        % Normalizar y almacenar el vector en Q
        R(j, j) = norm(v);
        Q(:, j) = v / R(j, j);
    end
    
    QB = Q' * B;
    RQB=[R, QB];
    s=size(RQB);
    nfilas=s(1);
    ncolumnas=s(2);
    Y = sust_reg_sup(RQB(1:nfilas,1:ncolumnas-1), RQB(:,ncolumnas));
end

function [sol,niteraciones,error]=metodo_jacobi(A,B,x0,cotaerror)
    D=diag(diag(A));
    U=triu(A)-D;
    L=tril(A)-D;
    M=-inv(D)*(L+U);
    N=inv(D)*B;
    cont=1;
    xant=x0;
    xsig=M*xant+N;
    while norm(xsig-xant,inf)>cotaerror
        cont=cont+1;
        xant=xsig;
        xsig=M*xant+N;
    end
    sol=xsig;
    niteraciones=cont;
    error=norm(xsig-xant,inf);
end

function [sol,niteraciones,error]=metodo_gauss_seidel(A,B,x0,cotaerror)
    D=diag(diag(A));
    U=triu(A)-D;
    L=tril(A)-D;
    M=-inv(D+L)*U;
    N=inv(D+L)*B;
    cont=1;
    xant=x0;
    xsig=M*xant+N;
    while norm(xsig-xant,inf)>cotaerror
        cont=cont+1;
        xant=xsig;
        xsig=M*xant+N;
    end
    sol=xsig;
    niteraciones=cont;
    error=norm(xsig-xant,inf);
end

function [x, iter, error] = sor(A, b, x0, omega, tol, maxIter)
    n = size(A, 1);
    x = x0;
    iter = 0;
    while iter < maxIter
        for i = 1:n
            x(i) = (1 - omega) * x0(i) + (omega / A(i,i)) * (b(i) - A(i,1:i-1)*x(1:i-1) - A(i,i+1:n)*x0(i+1:n));
        end
        error = norm(x - x0, inf);
        if error < tol
            return;
        end
        x0 = x;
        %x0'
        pause(0.1)
        iter = iter + 1;
    end
end