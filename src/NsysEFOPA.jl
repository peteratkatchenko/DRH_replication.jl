module NsysEFOPA_

  export NsysEFOPA

  function NsysEFOPA(X) 

    h1 = ubarEFO + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(X^ksi)
    h2 = (1-theta)*(((RTFP*(X^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    f = (pi/(kappa^2))*((log(h2)-h1)/(((1+psi)/h2)+(2/3)*exp(shocks[i,4])))^2-X
    

    return f

  end             
    
  
end #End NsysEFOPA