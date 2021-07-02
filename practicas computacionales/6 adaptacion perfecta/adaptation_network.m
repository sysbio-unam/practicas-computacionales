function dydt = adaptation_network(t,y, k1, k2, k3, k4,S)

%R=y(1)
%X=y(2)
dydt =[k1*S-k2*y(2)*y(1); % dR/dt=f1(R, X)
       k3*S-k4*y(2)];     % %dX/dt=f2(R,X)
   
   % note that second term is uncoupled