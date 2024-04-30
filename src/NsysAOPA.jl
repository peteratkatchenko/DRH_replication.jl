module NsysAOPA_

  export NsysAOPA 

  function NsysAOPA(X)
          
    g1 = ubarAO + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,3])/(shocks[i,5]*1000)^ksi)*(X^ksi))
    g2 = (1-theta)*(((RTFP*(X^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
        
    f = (pi/(kappa^2))*((log(g2)-g1)/(((1+psi)/g2)+(2/3)*AEF))^2-X
 
    return f

  end 
            
    
end #End of module NsysAOPA 