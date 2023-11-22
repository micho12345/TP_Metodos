%Programa para el metodo de biseccion
%a,b extremos del intervalo
function [sol,niter,error]=biseccion(a,b,cotaerror, maxiter)
cont=0;
xant=a;
xsig=b;
while (abs(xsig-xant)>cotaerror) && (cont < maxiter)
    cont=cont+1;
    if fun(xant)*fun(xsig)<0
        xmed=(xant+xsig)/2;
        if fun(xmed)*fun(xant)<0
            xsig=xmed;
        elseif fun(xsig)*fun(xmed)<0
            xant=xmed;
        else
            sol=xmed;
            error=0;
            niter=cont;
            break
        end
    elseif fun(xant)*fun(xsig)>0
        disp('Intervalo Incorrecto')
        break
    elseif fun(xant)==0
       sol=xant;
       error=0;
       niter=cont;
       break
   else
       sol=xsig;
       error=0;
       niter=cont;
       break
   end
end
sol = xmed;
niter = cont;
error = abs(xsig - xant);
% Muestra los resultados
fprintf('La solución aproximada es: %f\n', sol);
fprintf('Número de iteraciones: %d\n', niter);
fprintf('Error estimado: %e\n', error);
end

function y = fun(x)     %funcion f(x)
    y = exp(-(x^2)) - 1 + (1/((x^2)+1));
end


