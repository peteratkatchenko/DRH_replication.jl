function D = NsysEFOPA(X)

global i shocks kappa psi ubarEFO theta inter eps RTFP AA ksi

  N = X;  
    
  c1 = ubarEFO + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(N^ksi);
  c2 = (1-theta)*(((RTFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))));
  
  D(i) = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*exp(shocks(i,4))))^2-N;

   
D = [D(i)];
        

