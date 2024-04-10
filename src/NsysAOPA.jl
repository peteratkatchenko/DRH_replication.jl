module NsysAOPA_

  export NsysAOPA 

  function NsysAOPA(X::Vector, dictmain::Dict)
    i = dictmain[:i]
    shocks = dictmain[:shocks]
    kappa = dictmain[:kappa]
    psi = dictmain[:psi]
    ubarAO = dictmain[:ubarAO]
    theta = dictmain[:theta]
    inter = dictmain[:inter]
    eps = dictmain[:eps]
    RTFP = dictmain[:RTFP]
    AEF = dictmain[:AEF]
    ksi = dictmain[:ksi]

    N = X
          
    c1 = ubarAO + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,3])/(shocks[i,5]*1000)^ksi)*(N^ksi))
    c2 = (1-theta)*(((RTFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
        
    D[i] = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*AEF))^2-N
      
        
    D = [D[i]]
    
    return D

  end 
            
    
end #End of module NsysAOPA 