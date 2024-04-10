module NsysTFPPA_

  export NsysTFPPA

  function NsysTFPPA(X::Vector, i::Int, dictmain::Dict)

    shocks = dictmain[:shocks]
    kappa = dictmain[:kappa]
    psi = dictmain[:psi]
    ubarTFP = dictmain[:ubarTFP]
    theta = dictmain[:theta]
    inter = dictmain[:inter]
    eps = dictmain[:eps]
    RTFP = dictmain[:RTFP]
    ksi = dictmain[:ksi]

  
    
    N = X  
    
    c1 = ubarTFP + (1 + psi)*log((1 + psi))-psi*log(psi)-((exp(shocks[i,3])/(shocks[i,5]*1000)^ksi)*(N^ksi))

    c2 = (1-theta)*(((RTFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    D[i] = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*exp(shocks[i,4])))^2-N
  
      
    D = [D(i)]

    return D

  end           
    
end #End of module NsysTFPA