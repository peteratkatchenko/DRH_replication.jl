module NsysNSPA_

  export NsysNSPA

  function NsysNSPA(X::Vector, dictmain::Dict)

    i = dictmain[:i]
    kappa = dictmain[:kappa]
    psi = dictmain[:psi]
    ubarNS = dictmain[:ubarNS]
    theta = dictmain[:theta]
    inter = dictmain[:inter]
    eps = dictmain[:eps]
    RTFP = dictmain[:RTFP]
    AA = dictmain[:AA]
    AEF = dictmain[:AEF]
    ksi = dictmain[:ksi]


    N = X
      
    c1 = ubarNS + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(N^ksi);
    c2 = (1-theta)*(((RTFP*(N^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    D[i] = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*AEF))^2-N
  
    
    D = [D(i)]

    return D

  end               

end #End of module NsysNSPA