%Ejercicio 9

%a)
function y=metodo_euler(a,b,y0,N)
h=(b-a)/N;
t=a:h:b;
lt=length(t);
y=zeros(1,lt);
y(1)=y0;
for i=1:N
   y(i+1)=y(i)+h*Fed(t(i),y(i));
end
plot(t,y)
xlabel('Tiempo')
ylabel('Velocidad')
title('Grafico de velocidad en funcion del tiempo')
end

function v=Fed(t,v)
v=9.8-0.18*(v+8.3*(v/46)^2.2);
end

a=0
b=15
y0=0
N=150


x=metodo_euler(0,15,0,150);
y=x(151)

%b)
function y=metodo_runge_kutta_4(a,b,y0,N)
h=(b-a)/N;
t=a:h:b;
lt=length(t);
y=zeros(1,lt);
y(1)=y0;
for i=1:N
    k1=Fed(t(i),y(i));
    k2=Fed(t(i)+0.5*h,y(i)+0.5*k1*h);
    k3=Fed(t(i)+0.5*h,y(i)+0.5*k2*h);
    k4=Fed(t(i)+h,y(i)+k3*h);
    y(i+1)=y(i)+(h/6)*(k1+2*k2+2*k3+k4);
end
plot(t,y)
xlabel('Tiempo')
ylabel('Velocidad')
title('Grafico de velocidad en funcion del tiempo')
end


x=metodo_runge_kutta_4(0,15,0,150);
y=x(151)

%c1)
function dydt = odefunc(t, y)
    dydt = 9.8 - 0.18 * (y + 8.3 * (y / 46) ^ 2.2);
end

% Condiciones iniciales y tiempo de integración
tspan = [0 20];  

% Resolviendo la ecuación diferencial
[t, y] = ode45(@odefunc, tspan, y0);

% Grafico la solución
plot(t, y)
xlabel('Tiempo (s)')
ylabel('Velocidad (m/s)')
title('Solución de la ecuación diferencial usando ode45')

x= y(length(y))

%c2)
function dydt = odefunc(t, y)
    dydt = 9.8 - 0.18 * (y + 8.3 * (y / 46) ^ 2.2);
end

% Condiciones iniciales y tiempo de integración
tspan = [0 20];  

% Resolviendo la ecuación diferencial
[t, y] = ode23(@odefunc, tspan, y0);

% Grafico la solución
plot(t, y)
xlabel('Tiempo (s)')
ylabel('Velocidad (m/s)')
title('Solución de la ecuación diferencial usando ode45')

x= y(length(y))
