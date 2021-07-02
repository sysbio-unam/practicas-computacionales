function x=Optimization_adaptation_network
 tic 
close all
clear all
clc

%I wrote this example inside a function, so that all
% the necessary m-files are in the same place.

%% Experimental data
% Input experimental data (if you have large amounts of data, it might be
% better to store it in a sepparate .mat file and load it here, using the
% command load('name_of_data_file.mat').

t_exp=[   0    0.5000    1.0000    2.0000    4.0000    6.0000 ];%  24.0000]; % Time of measurement, in hours
y_exp=[0.0126    1.6059    0.8196    0.7323    0.5339    0.5613 ];%   0.1898]; % pStatStat3


%% Mathematical model (its the same as the one inside the CostFunction)
%B_t_function=@(B_0, k_bp, t)  1+  exp(-k_bp*t)*(B_0 -1);

%% Optimization

% Give an initial guess for the parameter

k10=.5;
k20=.5;
k30=.5;
k40=.5;


S=2.5;

xinit= [k10 k20 k30 k40];
% Set some options

%options= optimset('Algorithm','active-set','TolFun',1e-4,'TolX',1e-6,'MaxIter',1e4,'MaxFunEvals',1e6);
options= optimset('TolFun',1e-4,'TolX',1e-6,'MaxIter',1e4,'MaxFunEvals',1e6);

% Run the optimization!
[x,J,flag]=fminsearch(@(x)CostFunction(x,S, t_exp, y_exp),xinit,options);


%% Visualize if the optimization is ok
%Run the model with the parameter we found
tspan=t_exp(1):0.1:t_exp(end)+1;
y0=[0 0 0]; % R(0), X(0)

k1=x(1);
k2=x(2);
k3=x(3);
k4=x(4);



[t_opt,y_opt] = ode45(@(t,y)adaptation_network(t,y, k1, k2, k3, k4,S),tspan,y0);

Cost=CostFunction(x,S, t_exp, y_exp)

%figure();
plot(t_opt,y_opt(:,1), 'b',t_exp, y_exp, 'bo','MarkerSize', 10,'LineWidth',2)
xlabel('Time [h]', 'FontSize', 12);
ylabel('pStat3/Stat3','FontSize', 12);
title(['Ks=' num2str(x) ' cost=' num2str(Cost)],'FontSize', 12);
set(gcf, 'Position', [100 100 500 500]); 
axis square
xlim([0, t_opt(end)]);
%ylim([0, 1]);


toc

end


function Cost=CostFunction(x,S, t_exp, y_exp)


k1=x(1);
k2=x(2);
k3=x(3);
k4=x(4);


%% Solve the ODE - with these parameters

%initial conditions

X_0=0; Y_0=0;
y0=[X_0 Y_0]; %initial condition [X, Y]
%integration interval
tspan = [0 t_exp(end)]; %integration interval
% Call the ODE solver 
[t,y] = ode45(@(t,y)adaptation_network(t,y, k1, k2, k3, k4, S),tspan,y0);

y_pre=y(:,1);

% interpolate those values corresponding to the measurements
Y_pre_Interpol= interp1(t,y_pre,t_exp);

%% Calculate the cost of the predition vs. the experimental data
Cost=(sum((Y_pre_Interpol-y_exp).^2));

end


function dydt = adaptation_network(t,y, k1, k2, k3, k4,S)

%R=y(1)
%X=y(2)
dydt =[k1*S-k2*y(2)*y(1); % dR/dt=f1(R, X)
       k3*S-k4*y(2)]     % %dX/dt=f2(R,X)

   
end  % note that second term is uncoupled
