function Practica_Angeli

close all
clear all
clc
 
%%%%% Ejemplo de un sistema bi-estable%%%%%%%%
%  Angeli D, Ferrell JE, Sontag ED. Detection of multistability, bifurcations, and hysteresis in a large class of biological positive-feedback systems. PNAS [Internet] 2004;101:1822–7. Available from: http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=357011&tool=pmcentrez&rendertype=abstract

%% Pregunta A)	Usando dos condiciones iniciales para cada condición (y0_a=[x1(0), y1(0)]=[0,0] y y0_b=[0,9]): 
% i) Generar tres diagramas de tiempo vs- variable y1, para tres valores del parámetro de bifurcación v=0.75, 1 y 1.9.  
% ii) Graficar esta misma información en un diagrama-fase; i.e. x1(t) vs y1(t). 
%(opcional: añadir a este diagrama el campo vectorial, dado por [dx1(t)/dt, dy1(t)/dt]
%(?6 figuras).

% constant parameter values:
% Parametros:
alpha1=1; alpha2=1; beta1=200; beta2=10;gamma1=4;
gamma2=4; K1=30; K2=1;%v=1;

% Primero, un análisis dinámico - integración numérica del sistema, 
tspan = [0 10];
y0A = [0  0];
y0B = [0  0.9];

% for the vector plot
[xP,yP] = meshgrid(0:0.1:1,0:0.1:1);

figure_num=1;
ii=1;
for v=[0.75, 1, 1.9];%1.9; % 1.9 only the high steady state %1; bistable 0.75 only the low one
 
[tl,yl] = ode45(@(t,y)Angeli(t,y,v),tspan,y0A);
[th,yh] = ode45(@(t,y)Angeli(t,y,v),tspan,y0B);

subplot(2,3,figure_num)
plot(tl,yl(:,2), 'LineWidth', ii,'Color', 'k')
hold on
plot(th,yh(:,2), 'LineWidth', ii,'Color', 'r')
legend('y0 = [0  0]', 'y0 = [0  9]');
scatter(tl(end), yl(end,2), 'k', 'filled');
scatter(th(end),yh(end,2), 'r', 'filled');
xlabel('t');
ylabel('y(t)');
title(['Dinámica de y, v=' num2str(v)])
 %%
subplot(2,3,3+figure_num)
plot(yl(:,1),yl(:,2), 'LineWidth', ii,'Color', 'k')
hold on
plot(yh(:,1),yh(:,2), 'LineWidth', ii,'Color', 'r')
scatter(yl(end,1), yl(end,2), 'k', 'filled');
scatter(yh(end,1), yh(end,2), 'r', 'filled');
dx1dt=alpha1.*(1-xP)-beta1.*xP.*(v.*yP).^gamma1./(K1+(v.*yP).^gamma1);  
dy1dt=alpha2.*(1-yP)-beta2.*yP.*xP.^gamma2./(K2+xP.^gamma2);
quiver(xP, yP,dx1dt,dy1dt,4, 'b')
xlabel('x1(t) ');
ylabel('y1(t)');
title(['diagrama de fase, v=' num2str(v)])
axis square
xlim([0,1]);
ylim([0,1]);


figure_num=figure_num+1;
end

%% Pregunta B: Caracterizar las separatrices
% Regresemos al valor del parámetro de bifucación v=1.
% probando varias condiciones iniciales, grafica, en una misma figura, las
% trayectorias x1(t) vs y1(t). ¿te puedes dar una idea del tamaño de las
% cuencas de atracción de cada uno de los atractores? ¿qué sucede si
% aumentas el valor de v a 1.6? 

jj=1;
figure;
for v=[1, 1.6]
subplot(1,2,jj); 
hold on

for ii=1:1:100
    
    % generate 3 random numbers
    r1=rand;
    r2=rand;
    r3=rand;
    
    if r1<0.5    
    
    y0 = [r2<0.5 r3];
    else
    y0 = [r3 r2<0.5];
    end;
 
[~,yBi] = ode45(@(t,y)Angeli(t,y,v),tspan,y0);
 
if abs(yBi(end,2)-0.15)<0.2
 plot(yBi(:,1),yBi(:,2),'k')
 
else
 plot(yBi(:,1),yBi(:,2),'r')
   
end
 
end;
 
xlabel('x(0)');
ylabel('y(0)');
title(['diagrama de fase, v=' num2str(v)])

[~,yll] = ode45(@(t,y)Angeli(t,y,v),tspan,y0A);
[~,yhh] = ode45(@(t,y)Angeli(t,y,v),tspan,y0B);
scatter(yll(end,1), yll(end,2), 'k', 'filled');
scatter(yhh(end,1), yhh(end,2), 'r', 'filled');
dx1dt=alpha1.*(1-xP)-beta1.*xP.*(v.*yP).^gamma1./(K1+(v.*yP).^gamma1);  
dy1dt=alpha2.*(1-yP)-beta2.*yP.*xP.^gamma2./(K2+xP.^gamma2);
quiver(xP, yP,dx1dt,dy1dt,4, 'b')
axis square
xlim([0,1]);
ylim([0,1]);
jj=jj+1;
end;
%% Pregunta (C)
% Los sistemas bifurcantes presentan un comportamiento característico,
% llamado "alentamiento crítico", cuando se acercan al punto de
% bifurcación. ¿podemos ver algo así en este sistema? para comprobarlo,
% grafica, en la misma figura, t vs y1(t) para v=0.2:0.1:1. ¿qué observas? 
LineWidth=0.1;
LineStyle=[':', '--',  ':','--', ':', '--', ':','--', ':','--'];
figure;
ii=1;
col=['k', 'k',  'b','b', 'g', 'g', 'm','m', 'r', 'r'];
for v=0.2:0.1:1;
legendInfo{ii} = ['v = ' num2str(v)];
    
%[tl,yl] = ode45(@(t,y)Angeli(t,y,v),tspan,y0A);
[th,yh] = ode45(@(t,y)Angeli(t,y,v),tspan,y0B);

%plot(tl,yl(:,2), 'LineWidth', LineWidth,'Color', 'k')
hold on
plot(th,yh(:,2), 'LineWidth', LineWidth,'Color', col(ii), 'LineStyle', LineStyle(ii));
LineWidth=LineWidth+1;
ii=ii+1;
end;
legend(legendInfo)
xlabel('t');
ylabel('y(t)');
title('Critical slowing down');
text(0, 0.15, 'Recovery time from a transient perturbation increases')
text(0, 0.1,'as system approaches the bifurcation point'); 
text(0, 0.9, 'I.C: Transient perturbation'); 
text(10, 1, ':( Unhealthy stable state'); 
text(10, 0.18, ':) Healthy stable state'); 

axis square

function dydt = Angeli(~,y,v)

% Parametros:
alpha1=1; alpha2=1; beta1=200; beta2=10;gamma1=4;
gamma2=4; K1=30; K2=1;%v=1;

dydt = zeros(2,1);

dydt(1) = alpha1*(1-y(1))-beta1*y(1)*(v*y(2))^gamma1/(K1+(v*y(2))^gamma1);  
dydt(2) = alpha2*(1-y(2))-beta2*y(2)*y(1)^gamma2/(K2+y(1)^gamma2);
