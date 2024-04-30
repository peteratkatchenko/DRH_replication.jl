module NsysEFPA_

  export NsysEFPA

  function NsysEFPA(X)
    
    ATFP = exp(shocks[i,2])/((shocks[i,5]*1000)^eps)

    c1 = ubarEF + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,3])/(shocks[i,5]*1000)^ksi)*(X^ksi))

    c2 = (1-theta)*(((ATFP*(X^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    f = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*AEF))^2-X
    
    return f

  end 

end #End of module NsysEFPA   