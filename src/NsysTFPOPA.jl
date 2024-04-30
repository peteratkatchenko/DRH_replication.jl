module NsysTFPOPA_

  export NsysTFOPA

  function NsysTFPOPA(X)

    ATFP = exp(shocks[i,2])/((shocks[i,5]*1000)^eps)
    e1 = ubarTFPO + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(X^ksi)
    e2 = (1-theta)*(((ATFP*(X^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    f = (pi/(kappa^2))*((log(e2)-e1)/(((1+psi)/(e2))+(2/3)*AEF))^2-X

    return f

  end 
            
    
end #End of module NsysTFOPA