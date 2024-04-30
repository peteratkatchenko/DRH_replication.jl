module NsysTFPPA_

  export NsysTFPPA

  function NsysTFPPA(X)

    c1 = ubarTFP + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,4])/(shocks[i,1]*1000)^ksi)*(X^ksi))

    c2 = (1-theta)*(((RTFP*(X^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    f = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*exp(shocks[i,3])))^2-X
  
    return f

  end           
    
end #End of module NsysTFPA