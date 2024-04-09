function D = NsysTFPOPA(X)

global i shocks kappa psi ubarTFPO theta inter eps AEF AA ksi
  N = X;  
    
  ATFP = exp(shocks(i,2))/((shocks(i,5)*1000)^eps);  
  c1 = ubarTFPO + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(N^ksi);
  c2 = (1-theta)*(((ATFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))));
  
  D(i) = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/(c2))+(2/3)*AEF))^2-N;

   
D = [D(i)];
        

