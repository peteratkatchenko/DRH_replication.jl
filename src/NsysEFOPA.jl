module NsysEFOPA_

  export NsysEFOPA

  function d = NsysEFOPA(X::Vector, dictmain::Dict)

    i = dictmain[:i]
    shocks = dictmain[:shocks]
    kappa = dictmain[:kappa]
    psi = dictmain[:psi]
    ubarEFO = dictmain[:ubarEFO]
    theta = dictmain[:theta]
    inter = dictmain[:inter]
    eps = dictmain[:eps]
    RTFP = dictmain[:RTFP]
    AA = dictmain[:AA]
    ksi = dictmain[:ksi]

    N = X
      
    c1 = ubarEFO + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(N^ksi)
    c2 = (1-theta)*(((RTFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    D(i) = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*exp(shocks[i,4])))^2-N
    
      
    D = [D(i)]

    return D

  end             
    
  
end #End NsysEFOPA