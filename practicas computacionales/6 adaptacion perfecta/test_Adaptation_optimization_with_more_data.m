close all
clear all
clc

%% now we test if this works with the new experiment

% lo que sacó la optimización

costo_final= 0.0136;
% parámtros obtenidos por medio de optimización paramétrica

% % 
%  k1= 6.0903;
%  k2=3.3667;
%  k3=3.5561;
%  k4=1.1144;
%  

k1=4.1076;
k2 =6.5134;
k3= 1.2829;
k4= 1.1219 ;


% ahora hacemos otro experimento.


%% Experimental data
% Input experimental data (if you have large amounts of data, it might be
% better to store it in a sepparate .mat file and load it here, using the
% command load('name_of_data_file.mat').

t_exp_1=[   0    0.5000    1.0000    2.0000    4.0000    6.0000 ];%  24.0000]; % Time of measurement, in hours
y_exp_1=[0.0126    1.6059    0.8196    0.7323    0.5339    0.5613 ];%   0.1898]; % pStatStat3


S=2.5;


%% Visualize if the optimization is ok
%Run the model with the parameter we found
tspan=[0 6];
y0=[.1 .1]; % R(0), X(0)


[t_1,y_1] = ode45(@(t,y)adaptation_network(t,y, k1, k2, k3, k4,S),tspan,y0);


%figure();
plot(t_1,y_1(:,1), 'r',t_exp_1, y_exp_1, 'ro','MarkerSize', 10,'LineWidth',2)
xlabel('Time [h]', 'FontSize', 12);
ylabel('pStat3/Stat3','FontSize', 12);

set(gcf, 'Position', [100 100 500 500]); 
axis square
%xlim([0, t_1(end)]);

%ylim([0, 1]);


%% experiment 2

%% Figure 3B: pSTAT3, and STAT3 were analyzed in cell lysates of SW480 cells treated with
%rhIFN?. SW480 cells were platted in 6 well plates and incubated with cytokines for 30
%minutes, and then media removed and replaced with fresh media. Cells were collected
%at 0.5, 1, 2, 4, 6 h after new media was added. GAPDH was used as loading control.
%n=3.


%time post treatment;	pStat3/Stat3
% Fig3B= ...
% [0	0.076738192
% 0.5	1.286765985
% 1	1.347379532
% 2	1.250921412
% 4	0.396501979];


Fig3B= ...
[0 1.286765985
0.5 1.347379532
1 1.250921412
2 0.396501979
4 0.473782901];


%%
%t_exp_2=[0    0.5000    1.0000    2.0000    4.0000]; % Time of measurement, in hours
%y_exp_2=[ 0.0767    1.2868    1.3474    1.2509    0.3965];%   0.1898]; % pStatStat3

%t_exp_2=[    0.5000    1.0000    2.0000    4.0000]; % Time of measurement, in hours
%y_exp_2=[    1.2868    1.3474    1.2509    0.3965];%   0.1898]; % pStatStat3


t_exp_2=[ 0   0.5000    1.0000    2.0000    4.0000]; % Time of measurement, in hours
y_exp_2=[1.286765985 1.347379532 1.250921412 0.396501979 0.473782901];


hold on
plot(t_exp_2, y_exp_2, 'bd','MarkerSize', 10,'LineWidth',2)
% so we have to run the experiment in parts:
% for 0.5 hours assuming high IFNg, i.e. 2.5
%%
tspan_2A=[0 .5];
y0_2A=[.1 .1]; % R(0), X(0)

[t_2A,y_2A] = ode45(@(t,y)adaptation_network(t,y, k1, k2, k3, k4,S),tspan_2A,y0_2A);
%%
tspan_2B=[0 6];
y0_2B=y_2A(end,:); % R(0), X(0)
S=.1;

[t_2B,y_2B] = ode45(@(t,y)adaptation_network(t,y, k1, k2, k3, k4,S),tspan_2B,y0_2B);



plot(t_2A-0.5,y_2A(:,1), '--g','LineWidth',2)

plot(t_2B,y_2B(:,1), '-b','LineWidth',2)

%%
title('Adaptation is not enough','FontSize', 12);

legend('optimized adaptation model (to data from experiment 1)',...
       'Experiment 1','Experiment 2','Adaptation model simulation of exp 2 (pre)',...
       'Adaptation model simulation of exp 2 (exp)');
   
   
   %%%%%% Necesitamos ALGO, algun mecanismo, para que el equilibrio sea
   %%%%%% sea diferente de 0 aun cuando S=0
    
