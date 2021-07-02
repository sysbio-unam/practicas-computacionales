function dydt = negative_feedback_network(t,y, k1, k2, k3, k4,S)


dydt =[k1*S-k2*y(2)*y(1); 
       k3*y(1)-k4*y(2)]; 
   
