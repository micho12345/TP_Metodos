%Funcion que genera y grafica el polinomio interpolador.
%x:vector de primeras coordenadas de los datos debe ser columna
%y:vector de segundas coordenadas de los datos debe ser columna
function p=interpolador(x,y)
V=vander(x);
sol=V\y';
p=sol';
%graficas
plot(x,y,'*k','linewidth',9)
grid minor
hold on
xlabel('x')
ylabel('y')
title('Polinomio Interpolador de datos')
a=min(x);
b=max(x);
t=a:0.01:b;
gp=polyval(p,t);
%z=exp(t);
plot(t,gp,'r','linewidth',3)
%plot(t,z,'g')