
# This program computes the counter-factual distributions in SA

#Select Country 

China = 0 #Make 0 if US and 1 if China

if China == 0
    
    #US Parameters
    psi = 1.4841
    theta = .3358
    totalhours = 5110
    beta = 0.98
    delta = 0.02
    inter = 0.02
    eps = .00 #(called omega in main text)
    ksi = .00
    ubar = 10
    kappa = 0.002
    USBench= CSV.read("C:\\Users\\peter\\.julia\\dev\\dev_econ_replication\\replication_files\\urban_accounting_welfare_replication\\ReplicationFiles\\MatlabPro\\USBenchmark.txt") 
    shocks = USBench.USBenchmark

else
    
    #China Parameters
    psi = 1.52;
    theta = .52;
    totalhours = 5110;
    beta = 0.90;
    delta = 0.2;
    inter = 0.2;
    eps = .00;
    ksi = .00;
    ubar = 10;
    kappa = 0.001;
    load CSV.read(ChinaBenchmark.txt 
    shocks = ChinaBenchmark;

end

Nbar = sum(1000*shocks(:,5));
Nlar = max(1000*shocks(:,5));

NbarNew = 0;

N = 1000*shocks(:,5)';

    for i = 1:length(shocks(:,1));
        c1(i) = ubar + (1 + psi)*log((1 + psi))-psi*log(psi);
        c2(i) = (1-theta)*((exp(shocks(i,2))^(1/(1-theta)))/((theta/inter)^(theta/(1-theta))));
    
        shocks(i,3) = -log(c2(i))+((N(i)*((kappa^2)/pi))^(1/2))*(((1+psi)/c2(i))+(2/3)*exp(shocks(i,4))) + c1(i);
        shocks(i,3) = log(shocks(i,3));
    end 
    
    
NbarNew = sum(N);


# Counterfactual efficiency
count = 0; 
step = 0.1 + China;
NbarNewTFP = 0;
global ubarTFP RTFP
ubarTFP = ubar;

while abs(Nbar-NbarNewTFP) > 0.02

    for i = 1:length(shocks(:,1));
       
        if China == 0
        
            RTFP = sum(exp(shocks(:,2))./((shocks(:,5)*1000).^eps).*(shocks(:,5)*1000/Nbar));
            
        else
            
            RTFP = prctile(exp(shocks(:,2))./((shocks(:,5)*1000).^eps),49);
            
        end    
        
        if NbarNewTFP == 0;
            X0 = 1000*shocks(i,5);
        else
            X0 = max(NTFP(i),1);
        end
        
        F=@(R)NsysTFPPA(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        NTFP(i) = max(0,real(X));
        
        
    end 
 
    NbarNewTFP = sum(NTFP);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if Nbar-NbarNewTFP > 0.02
        ubarTFP = ubarTFP - step;
    else
        ubarTFP = ubarTFP + step;
    end
    disp(Nbar-NbarNewTFP);
end

NbarNewTFP
max(NTFP)
ubarTFP

# Counterfactual ability
count = 0;
step = 0.1 + China;
NbarNewA = 0;
global ubarA AA
ubarA = ubar;

while abs(Nbar-NbarNewA) > 0.02

    for i = 1:length(shocks(:,1));
        
        if China == 0
        
            AA = sum(exp(shocks(:,3))./((shocks(:,5)*1000).^ksi).*shocks(:,5)*1000/Nbar);
            
        else
            
            AA = prctile(exp(shocks(:,3))./((shocks(:,5)*1000).^ksi),49);
            
        end    
        
        if NbarNewA == 0
            X0 = 1000*shocks(i,5);
        else
            X0 = max(NA(i),1);
        end
        
        F=@(R)NsysAPA(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        NA(i) = max(0,real(X));
            
    end 
 
    NbarNewA = sum(NA);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if Nbar-NbarNewA > 0.02
        ubarA = ubarA - step;
    else
        ubarA = ubarA + step;
    end
    disp(Nbar-NbarNewA);
end

NbarNewA
max(NA)
ubarA

# Counterfactual EF
count = 0;
step = 0.1;
NbarNewEF = 0;
global ubarEF AEF
ubarEF = ubar;

while abs(Nbar-NbarNewEF) > 0.02

    for i = 1:length(shocks(:,1));
        
        if China == 0
        
            AEF = sum(exp(shocks(:,4)).*shocks(:,5)*1000/Nbar);
            
        else
            
            AEF = prctile(exp(shocks(:,4)),49);
            
        end    
        
        if NbarNewEF == 0
            X0 = 1000*shocks(i,5);
        else
            X0 = max(NEF(i),1);
        end
        
        F=@(R)NsysEFPA(R);
        options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
        [X,Eval] = fsolve(F,X0,options);
    
        NEF(i) = max(0,real(X));
        
    end 
 
    NbarNewEF = sum(NEF);
    count = count + 1;
    if count > 10
        step = step/10;
        count = 0;
    end
    
    if Nbar-NbarNewEF > 0.02
        ubarEF = ubarEF - step;
    else
        ubarEF = ubarEF + step;
    end
     disp(Nbar-NbarNewEF);
end

NbarNewEF
max(NEF)
ubarEF

TMA = sum(max((NA - N),0))/Nbar
TMTFP = sum(max((NTFP - N),0))/Nbar
TMEF = sum(max((NEF - N),0))/Nbar

hold on 
SNR = sort(log(1000*shocks(:,5)));
SN = sort(log(N));
SNTFP = sort(log(NTFP));
SNA = sort(log(NA));
SNEF = sort(log(NEF));

for i = 1:length(N)
    lprob(i) = log((length(N) - i + 1)/length(N));
end

subplot(2,2,1)
hold on
grid on
plot(SNR,lprob,'r','LineWidth',2);
plot(SN,lprob,'b','LineWidth',2);
legend('Actual','Modeled','Location','SouthWest');
xlabel('ln(population)')
ylabel('ln(prob > population)')
title(['Model Utility = ',num2str(ubar)])


subplot(2,2,2)
hold on
grid on
plot(SN,lprob,'b','LineWidth',2);
plot(SNTFP,lprob,'g','LineWidth',2);
legend('Actual','Avg. Efficiency','Location','SouthWest');
xlabel('ln(population)')
ylabel('ln(prob > population)')
title(['Counterfactual Utility = ',num2str(ubarTFP) ', Reallocation = ' num2str(TMTFP)])

subplot(2,2,3)
hold on
grid on
plot(SN,lprob,'b','LineWidth',2);
plot(SNA,lprob,'m','LineWidth',2);
legend('Actual','Avg. Amenities','Location','SouthWest');
xlabel('ln(population)')
ylabel('ln(prob > population)')
title(['Counterfactual Utility = ',num2str(ubarA) ', Reallocation = ' num2str(TMA)])

subplot(2,2,4)
hold on
grid on
plot(SN,lprob,'b','LineWidth',2);
plot(SNEF,lprob,'k','LineWidth',2);
legend('Actual','Avg. Exc. Frictions','Location','SouthWest');
xlabel('ln(population)')
ylabel('ln(prob > population)')
title(['Counterfactual Utility = ',num2str(ubarEF) ', Reallocation = ' num2str(TMEF)])
subplot(2,2,1)
text(15.5,1,['Counterfactuals Without One Shock, \kappa = ' num2str(kappa) ' , \omega = ' num2str(eps) ' , \zeta = ' num2str(ksi)])

chA = ((NA - N)./N)';
chTFP = ((NTFP - N)./N)';
chEF = ((NEF - N)./N)';
PercChange = [chA chTFP chEF];
ShocksModel = shocks(:,1:5);

if China == 0 
    save PercChangeUS.txt PercChange -ascii
    save ShocksModelUS.txt ShocksModel -ascii
else
    save PercChangeChina.txt PercChange -ascii
    save ShocksModelChina.txt ShocksModel -ascii
end    
    
#Only US
if China == 0
 
    # Counterfactual TFP only
    count = 0;
    step = 0.1;
    NbarNewTFPO = 0;
    global ubarTFPO
    ubarTFPO = ubar;

    while abs(Nbar-NbarNewTFPO) > 0.02

        for i = 1:length(shocks(:,1));
        
            if NbarNewTFPO == 0
                X0 = 1000*shocks(i,5);
            else
                X0 = max(NTFPO(i),1);
            end
        
            F=@(R)NsysTFPOPA(R);
            options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
            [X,Eval] = fsolve(F,X0,options);
    
            NTFPO(i) = max(0,real(X));
        
        end 
 
        NbarNewTFPO = sum(NTFPO);
        count = count + 1;
        if count > 10
            step = step/10;
            count = 0;
        end
    
        if Nbar-NbarNewTFPO > 0.02
            ubarTFPO = ubarTFPO - step;
        else
            ubarTFPO = ubarTFPO + step;
        end
        disp(Nbar-NbarNewTFPO);
    end

    NbarNewTFPO
    max(NTFPO)
    ubarTFPO

    TMTFPO = sum(max((NTFPO - N),0))/Nbar;

    figure
    subplot(2,2,1)
    hold on
    grid on
    plot(SNR,lprob,'r','LineWidth',2);
    plot(SN,lprob,'b','LineWidth',2);
    legend('Actual','Modeled','Location','SouthWest');
    xlabel('ln(population)')
    ylabel('ln(prob > population)')
    title(['Model Utility = ',num2str(ubar)])

    subplot(2,2,2)
    hold on
    grid on
    SNTFPO = sort(log(NTFPO));
    plot(SN,lprob,'b','LineWidth',2);
    plot(SNTFPO,lprob,'g','LineWidth',2);
    legend('Actual','Efficiency Only','Location','SouthWest');
    xlabel('ln(population)')
    ylabel('ln(prob > population)')
    title(['Counterfactual Utility = ',num2str(ubarTFPO) ', Reallocation = ' num2str(TMTFPO)])

    # Counterfactual Amenities Only
    count = 0;
    step = 0.1;
    NbarNewAO = 0;
    global ubarAO
    ubarAO = ubar;

    while abs(Nbar-NbarNewAO) > 0.02

        for i = 1:length(shocks(:,1));
        
        
            if NbarNewAO == 0
                X0 = 1000*shocks(i,5);
            else
                X0 = max(NAO(i),1);
            end
        
            F=@(R)NsysAOPA(R);
            options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
            [X,Eval] = fsolve(F,X0,options);
    
            NAO(i) = max(0,real(X));
        
        end 
 
        NbarNewAO = sum(NAO);
        count = count + 1;
        if count > 10
            step = step/10;
            count = 0;
        end
    
        if Nbar-NbarNewAO > 0.02
            ubarAO = ubarAO - step;
        else
            ubarAO = ubarAO + step;
        end
        disp(Nbar-NbarNewAO);
    end

    NbarNewAO
    max(NAO)
    ubarAO

    TMAO = sum(max((NAO - N),0))/Nbar;

    subplot(2,2,3)
    hold on
    grid on
    SNAO = sort(log(NAO));
    plot(SN,lprob,'b','LineWidth',2);
    plot(SNAO,lprob,'m','LineWidth',2);
    legend('Actual','Amenities Only','Location','SouthWest');
    xlabel('ln(population)')
    ylabel('ln(prob > population)')
    title(['Counterfactual Utility = ',num2str(ubarAO) ', Reallocation = ' num2str(TMAO)])

    # Counterfactual Excessive Frictions Only
    count = 0;
    step = 0.1;
    NbarNewEFO = 0;
    global ubarEFO
    ubarEFO = ubar;

    while abs(Nbar-NbarNewEFO) > 0.02

        for i = 1:length(shocks(:,1));
        
            if NbarNewEFO == 0
                X0 = 1000*shocks(i,5);
            else
                X0 = max(NEFO(i),1);
            end
        
            F=@(R)NsysEFOPA(R);
            options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
            [X,Eval] = fsolve(F,X0,options);
    
            NEFO(i) = max(0,real(X));
        
        end 
 
        NbarNewEFO = sum(NEFO);
        count = count + 1;
        if count > 10
            step = step/10;
            count = 0;
        end
    
        if Nbar-NbarNewEFO > 0.02
            ubarEFO = ubarEFO - step;
        else
            ubarEFO = ubarEFO + step;
        end
        disp(Nbar-NbarNewEFO);
    end

    NbarNewEFO
    max(NEFO)
    ubarEFO

    TMEFO = sum(max((NEFO - N),0))/Nbar;

    subplot(2,2,4)
    hold on
    grid on
    SNEFO = sort(log(NEFO));
    plot(SN,lprob,'b','LineWidth',2);
    plot(SNEFO,lprob,'k','LineWidth',2);
    legend('Actual','Exc. Frictions Only','Location','SouthWest');
    xlabel('ln(population)')
    ylabel('ln(prob > population)')
    title(['Counterfactual Utility = ',num2str(ubarEFO) ', Reallocation = ' num2str(TMEFO)])

    subplot(2,2,1)
    text(15.5, 1,['Counterfactuals with Only One Shock, \kappa = ', num2str(kappa) ' , \omega = ' num2str(eps) ' , \zeta = ' num2str(ksi)])

    # Counterfactual Without Shocks
    count = 0;
    step = 0.1;
    NbarNewNS = 0;
    global ubarNS
    ubarNS = ubar;

    while abs(Nbar-NbarNewNS) > 0.02

        for i = 1:length(shocks(:,1));
        
            if NbarNewNS == 0
                X0 = 1000*shocks(i,5);
            else
                X0 = max(NNS(i),1);
            end
        
            F=@(R)NsysNSPA(R);
            options=optimset('MaxFunEvals',10000000000,'TolFun',0.0001,'MaxIter',1000000,'Display','off');
            [X,Eval] = fsolve(F,X0,options);
    
            NNS(i) = max(0,real(X));
        
        end 
 
        NbarNewNS = sum(NNS);
        count = count + 1;
        if count > 10
            step = step/10;
            count = 0;
        end
    
        if Nbar-NbarNewNS > 0.02
            ubarNS = ubarNS - step;
        else
            ubarNS = ubarNS + step;
        end
        disp(Nbar-NbarNewNS);
    end

    NbarNewNS
    max(NNS)
    ubarNS

    figure

    hold on
    grid on
    SNNS = sort(log(NNS));
    plot(SN,lprob,'b','LineWidth',2);
    plot(SNNS,lprob,'k','LineWidth',2);
    legend('Actual','No Shocks','Location','SouthWest');
    xlabel('ln(population)')
    ylabel('ln(prob > population)')
    title(['Counterfactual Utility = ',num2str(ubarNS)])

end 
