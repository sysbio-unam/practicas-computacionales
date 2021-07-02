close all
clear all
clc

%syms GFP(t) GFP0 k1 k2
%dsolve(diff(GFP) == k1-k2*GFP, GFP(0) == GFP0)

GFP_t=@(k1, k2, t, GFP0)(k1 - exp(-k2.*t).*(k1 - GFP0.*k2))./k2;
%% Producción y degradación de la GFP

% 0 --> GFP [k1]
% GFP --> 0 [k2]


k1=100; 
k2=5;

th=[k1,k2]; %) # parameters (rate constants)

h=@(y, th)[th(1), th(2)*y];

GFP0=5;
n=300;% # number of realizations (tsteps)
iterations=100; 
v=length(th);% # number of reactions (IN THIS CASE: same as number of kinetic parameters)
tvecMax=0;
for iteration=1:1:iterations
tt=0; % # initialize time
x=GFP0; % # initial condition
tvec=zeros(1,n); % # create a vector for the time, n number of realizations
xvec=zeros(1,n); % # vector for state variable 


for i=1:1:n
    totH= sum(h(x, th));%# total hazard
    %tt=tt+exprnd(totH);%sample time to next event from exp probability distribution
    %tt=tt+random('Exponential',totH,1);%log(1/rand).*(1/ totH);
     tt=tt+log(1/rand).*(1/ totH);
    % # choose a reaction to occur
    prob_reac=h(x, th)./totH;
    
    if rand<=prob_reac(1); % choose first reaction
        
    %# update state vector
      x=x+1; % # GFP formed
    else 
      x=x-1; % # decay
    end;
    
    % done, just fill in the vectors:
    tvec(i)=tt;
    xvec(i)=x;
    
end; 

time_vector_regular = 0:0.01:5; 
xt_matrix(iteration,:) = interp1( tvec,xvec, time_vector_regular); 
 
plot(time_vector_regular,xt_matrix(iteration, :));
%plot(tvec,xvec);

xt_final_matrix(iteration)=xvec(end);
hold on

tvecMax=max(tvecMax, tvec(end))
end;

plot(time_vector_regular,mean(xt_matrix),'color','b', 'LineWidth',2)



title(['Gillespie simulations, k1=' num2str(k1) ', k2=' num2str(k2) ', GFP0=' num2str(GFP0) ', realizations=' num2str(n) ', iterations=' num2str(iteration)]);

plot(tvec,GFP_t(k1, k2,tvec, GFP0),'color','r', 'LineWidth',2)

xlabel('time');
ylabel('GFP');
hold on
axis square
line([0, tvecMax], [k1/k2, k1/k2],'color','k', 'LineWidth',2);
xlim([0, tvec(end)]);


figure;
hi=histogram(xt_final_matrix);


%hi.BinWidth = 0.1;
hi.Normalization='Probability';
ylim([0,.5])
hold on 
line([k1/k2, k1/k2], [0,1],'Color','r', 'LineWidth',2)

axis square
xlabel('GFP final value')
ylabel('frequency')
%xlim([5 25])


figure;
plot(time_vector_regular,std(xt_matrix,1),'color','k', 'LineWidth',2)
xlabel('time')
ylabel('standard deviation')
xlim([0, tvec(end)]);
axis square
