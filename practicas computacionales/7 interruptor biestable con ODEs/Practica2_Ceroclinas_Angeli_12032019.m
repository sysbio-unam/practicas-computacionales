% close all
% clear all
% clc

%% encontrar las ceroclinas de la función de angeli


%dx/dt = alpha1*(1-x) -beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
% dx/dt=fx(x,y) --> fx(x, y) =0 --> dx/dt=0 --> x está en equilibrio.
% xss=g(y)

%dy= alpha1*(1-x) - beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
% análogo
% ...
% yss=h(x)

% g(y) y h(x) se llaman ceroclinas (nullclines) de x y de y,
% respectivamente.

alpha1=1; alpha2=1; beta1=200; beta2=10; gamma1=4; gamma2=4; K1=30; K2=1; %v=3;



for  v=[0.1, 1, 3]
    
%%


%%
%%% ceroclina de x: dx=0
%0 = alpha1*(1-x)-beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
%alpha1*(1-x)=beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
%(1-x)/x=beta1*(v*y)^gamma1/(K1+(v*y)^gamma1)/alpha1
%1/x-1=beta1*(v*y)^gamma1/(K1+(v*y)^gamma1)/alpha1
%1/x=(beta1*(v*y)^gamma1/(K1+(v*y)^gamma1)/alpha1)+1
%x=((beta1*(v*y)^gamma1/(K1+(v*y)^gamma1)/alpha1)+1)^(-1)

xss=@(y)1./((beta1*(v.*y).^gamma1./(K1+(v.*y).^gamma1)./alpha1)+1);
%%
  
%%% encontrar la ceroclina de y
% dy = alpha2*(1-y)-beta2*y*x^gamma2/(K2+x^gamma2)=0
% alpha2*(1-y)=beta2*y*x^gamma2/(K2+x^gamma2)
% (1-y)/y=beta2*x^gamma2/(K2+x^gamma2)/alpha2
% 1/y-1=beta2*x^gamma2/(K2+x^gamma2)/alpha2
% 1/y=beta2*x^gamma2/(K2+x^gamma2)/alpha2+1
yss=@(x)1./(beta2*x.^gamma2./(K2+x.^gamma2)./alpha2+1);


%%
ww=0:0.01:1;
   
%%
figure;
subplot(1,3,1)
plot(ww, xss(ww),'LineWidth', 2,  'Color', 'r');
xlabel('y');
ylabel('xss(y): ceroclina de x')
axis square
xlim([0, 1.1]);
ylim([0, 1.1]);

%%
subplot(1,3,2)
plot(ww, yss(ww),'LineWidth', 2,  'Color', 'b');
xlabel('x');
ylabel('yss(x): ceroclina de y')
axis square
xlim([0, 1.1]);
ylim([0, 1.1]);

%%
subplot(1,3,3)
plot(ww, yss(ww),'LineWidth', 2,  'Color', 'b');
hold on
plot(xss(ww), ww, 'LineWidth', 2,  'Color', 'r');
xlabel('x');
ylabel('y');
legend('yss(y): ceroclina de y', 'xss(y): ceroclina de x')
axis square
xlim([0, 1.1]);
ylim([0, 1.1]);
%%

end;


%%

J=Jacobian_matrix_Angeli(0.5, 0.62)
%% evaluar la estabilidad de las raíces

%%%%% ahora vamos a ver la estabiliad de los puntos de equilibrio, que
%%%%% obtuvimos de manera numérica en R

%stable
J=Jacobian_matrix_Angeli(0.9946980, 0.1681565);

% 
J=Jacobian_matrix_Angeli(0.5057743, 0.6195077)