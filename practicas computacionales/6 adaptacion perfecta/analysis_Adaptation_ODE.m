close all
clear all
clc

y0=[0 1]; % R(0), X(0)
tspan = [0 10];

k1=6.0903; k2=3.3667; k3=3.5561; k4=1.1144;

for S=0:1:5
[t_A,y_A] = ode45(@(t,y)adaptation_network(t,y, k1, k2, k3, k4,S),tspan,y0);

Rmax=max(y_A(:,1));

% los puntos de equilibrio de y1 y y2, notar cómo el de y1 es independiente
% de S
y1ss=@(k1, k2, k3, k4,S)k1*k4/(k2*k3);
y2ss=@(k1, k2, k3, k4,S)k3*S/k4;

subplot(1,2,1)
line([0 tspan(end)], [y1ss(k1, k2, k3, k4,S) y1ss(k1, k2, k3, k4,S)], 'Color', 'k');
hold on
plot(t_A, y_A(:,1),'b')
xlabel('Time t');
ylabel('y_1(t)') %, punto de equilibrio independiente de S (Adaptación perfecta)');


subplot(1,2,2)
line([0 tspan(end)], [y2ss(k1, k2, k3, k4,S) y2ss(k1, k2, k3, k4,S)], 'Color', 'r');
hold on
plot(t_A, y_A(:,2),'g')
xlabel('Time t');
ylabel('y_2(t)')%, controlador');

end

hold on 

Fig3A= ...
[0	0.012590557	0.118586476
0.5	1.60587993	1.465118316
1	0.819641849	1.941786253
2	0.732285575	2.15416639
4	0.533866523	1.562079164
6	0.561275345	1.156997778  ];
%24   [0.224603914	2.03716536]-[0.034795764	0.803038463]]; 

subplot(1,2,1)
plot(Fig3A(:,1),Fig3A(:,2), 'rs','LineWidth',2,'MarkerEdgeColor','r',...
                       'MarkerFaceColor','r',...
                       'MarkerSize',10)
hold on
title('Adaptación perfecta')




