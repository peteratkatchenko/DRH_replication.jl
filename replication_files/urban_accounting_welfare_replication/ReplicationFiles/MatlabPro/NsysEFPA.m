function D = NsysEFPA(X)

global i shocks kappa psi ubarEF theta inter eps AEF ksi

  N = X;
  
  ATFP = exp(shocks(i,2))/((shocks(i,5)*1000)^eps);  
  c1 = ubarEF + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks(i,3))/(shocks(i,5)*1000)^ksi)*(N^ksi));
  c2 = (1-theta)*(((ATFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))));
  
  D(i) = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*AEF))^2-N;

   
D = [D(i)];
        

