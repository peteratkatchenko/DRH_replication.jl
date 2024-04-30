module NsysNSPA_

  export NsysNSPA

  function NsysNSPA(X)
      
    j1 = ubarNS + (1 + psi)*log((1 + psi))-psi*log(psi)-AA*(X^ksi);
    j2 = (1-theta)*(((RTFP*(X^eps))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))))
    
    f = (pi/(kappa^2))*((log(c2)-c1)/(((1+psi)/c2)+(2/3)*AEF))^2-X
  

    return f

  end               

end #End of module NsysNSPA