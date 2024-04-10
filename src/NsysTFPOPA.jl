module NsysTFPOPA_

  export NsysTFOPA

  function NsysTFPOPA(X::Vector, dictmain::Dict)

      

    i = dictmain[:i]
    shocks = dictmain[:shocks]
    kappa = dictmain[:kappa]
    psi = dictmain[:psi]
    ubarTFPO = dictmain[:ubarTFPO]
    theta = dictmain[:theta]
    inter = dictmain[:inter]
    eps = dictmain[:eps]
    AEF = dictmain[:AEF]
    AA = dictmain[:AA]
    ksi = dictmain[:ksi]

    N = X
      
    ATFP = exp(shocks[i,2])/((shocks[i,5]*1000)^eps)
    c1 = ubarTFPO + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(N^ksi)
    c2 = (1-theta)*(((ATFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    D[i] = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/(c2))+(2/3)*AEF))^2-N

      
    D = [D(i)]

    return D

  end 
            
    
end #End of module NsysTFOPA