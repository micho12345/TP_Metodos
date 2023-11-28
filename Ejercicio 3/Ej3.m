%Creo una funcion para calcular y
function[y]= transformacion (x)
    y= 1./ (x.^2);
end
%Con los valores obtenidos se arma la matriz A
function[matriz]= armamatriz(y)
    n=length(y);
    matriz=[y,ones(n,1)];
end

matrizA = armamatriz(y)

%Hago una funcion para trasponer la matriz
function[A_transpose]= traspuesta(A)
    A_transpose= transpose(A)
end

matrizAT= traspuesta(matrizA)

%Calculo la matriz B
matrizB= matrizAT*matrizA

function[resultado]=producto_matricial(matriz,vector)
    resultado=matriz*vector;
end

%Obtengo el vectorD
vectorD= producto_matricial(matrizAT,k)

%Programa para factorizar una matriz por el metodo de Cholesky

function L = factorizacionCholesky(A)
    % Dimensiones de la matriz A
    n = size(A, 1);
    % Inicializar la matriz L
    L = zeros(n);
    % Realizar la factorización de Cholesky
    for j = 1:n
        for i = 1:j
            if i == j
                L(j, j) = sqrt(A(j, j) - sum(L(j,1:j-1).^2));
            else
                L(j, i) = (A(j, i) - sum(L(i, 1:i-1) * L(j, 1:i-1)')) / L(i, i);
            end
        end
    end
end

L=factorizacionCholesky(matrizB)

function x = despejo_cholesky (L, b)
    % Dimensiones de la matriz L y el vector b
    n = size(L, 1);
    % Resolución de Ly = b mediante sustitución hacia adelante
    y = zeros(n, 1);
    for i = 1:n
        y(i) = (b(i) - L(i, 1:i-1) * y(1:i-1)) / L(i, i);
    end
    % Resolución de L'x = y mediante sustitución hacia atrás
    x = zeros(n, 1);
    for i = n:-1:1
    x(i) = (y(i) - L(i+1:n, i)' * x(i+1:n))/ L(i,i);
    end
end

x=despejo_cholesky(L, vectorD)
%Para pronosticar la tasa de crecimiento para x=2 mg/L
function tasa_crecimiento(2,kmax,c)
    tasa=kmax * x^2 / (x^2+c)
end

tasa=tasa_crecimiento (2,kmax,c)