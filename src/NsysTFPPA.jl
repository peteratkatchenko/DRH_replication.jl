module NsysTFPPA_

  export NsysTFPPA

  function NsysTFPPA(D::Vector, X::Vector, i::Int64, dictmain::Dict, dictCE::Dict)
    psi = dictmain[:psi]
    ksi = dictmain[:ksi]
    theta = dictmain[:theta]
    eps = dictmain[:eps]
    inter = dictmain[:inter]
    kappa = dictmain[:kappa]
    shocks = dictmain[:shocks]

    ubarTFP = dictCE[:ubarTFP]
    RTFP = dictCE[:RTFP]

    N = X[1]

    c1 = ubarTFP + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,4])/(shocks[i,1]*1000)^ksi)*(N^ksi))

    c2 = (1-theta)*(((RTFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    D[1] = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*exp(shocks[i,3])))^2-N
  
    return D[1]

  end           
    
end #End of module NsysTFPA