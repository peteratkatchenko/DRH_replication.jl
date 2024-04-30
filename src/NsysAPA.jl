module NsysAPA_

  export NsysAPA

  function NsysAPA(X)
        
    ATFP = exp(shocks[i,2])/((shocks[i,5]*1000)^eps)  
    c1 = ubarA + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(X^ksi)
    c2 = (1-theta)*(((ATFP*(X^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
        
    f = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*exp(shocks[i,4])))^2-X

    return f
              
  end 

  end #End of module NsysAPA  