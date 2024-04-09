function D = NsysNSPA(X)

global i kappa psi ubarNS theta inter eps RTFP AA AEF ksi

  N = X;  
    
  c1 = ubarNS + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(N^ksi);
  c2 = (1-theta)*(((RTFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))));
  
  D(i) = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*AEF))^2-N;

   
D = [D(i)];
        

