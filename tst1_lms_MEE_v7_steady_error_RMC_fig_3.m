clear all; close all;
%% MEE-AF-V10
L = 2000;
p = 5;
q = 8;
lamda1 = 0.994;%RMEE
lamda2 = 0.992;%RMEE
lamda3 = 0.980;%RMEE
lamda4 = 0.970;%RMEE
lamda5 = 0.942;%RMEE
% lamda2 = 1;%RMEE for speech signal

sigma1= 0.2;%%RMEE
sigma2 = 0.5;%%RMEE
sigma3 = 1;%%RMEE
sigma4 = 2;%%RMEE
sigma5 = 3;%%RMEE
tic
for mm = 1:100
    vv = randn(1,L) * 0.1;
%     vv = rand(1,L) -0.5;
    v1=randn(1,L)*0.1; v2=randn(1,L)*5;
    rp=rand(1,L);
    %    vv = v1 + (rp>0.95).*v2;
    vv = (rp<=0.95).*v1 + (rp>0.95).*v2;
    vG = exp(-vv.^2/2/sigma1^2).*vv;
    %      vv = exp(-vv.^2/2/sigma1^2).*vv;
    
   wo1 = randn(p,1);wo2 = randn(p,1);wo3 = randn(p,1);
%       wo1 = [0,0,0.9,0,0,0,0.2,0,0,0]';
%     wo = [ kron(wo1, ones(1,L/3)) kron(wo2, ones(1,L/3)) kron(wo3, ones(1,L/3))];
     wo = [ kron(wo1, ones(1,L)) ];
  %   wo = [ kron(randn(p,1), ones(1,L/2)),  kron(randn(p,1), ones(1,L/2)) ];
    uu = randn(p,L);
%     path(path,'e:\work\speech enhancement\dbs');
%     [s1,FS,NBITS]=wavread('sp02.wav');
%     s1 = s1*1;
%     
%     for ii = 1 : p
%         uu(ii, :) = s1(3000+ii : 3000+ii+L-1)*10;
%     end
    for ii = 1 : L
        dd(ii) = wo(:,ii)' * uu(:,ii) + vv(ii);
    end
    
    % dd = wo' * uu + vv;
    
    w_LMS = randn(p,1);   
    w_C_RLS1 =w_LMS;%%MEE
    w_C_RLS2 =w_LMS;%%MEE
    w_C_RLS3 =w_LMS;%%MEE
    w_C_RLS4 =w_LMS;%%MEE
    w_C_RLS5 =w_LMS;%%MEE
    
    
    %% RLS MCC, algorithm 1
   Pn = eye(p)*1;
    for ii = 1 : L
        Err_MCC_RLS1(mm,ii) = (wo(:,ii)  - w_C_RLS1)' * (wo(:,ii)  - w_C_RLS1);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_C_RLS1' * un;
        
        kn = Pn * un / ( exp(en^2/2/sigma1^2)*lamda1 + un' * Pn * un );
        %  kn = Pn * un / ( lamda+ un' * Pn * un );
        Pn = 1/lamda1 * ( Pn - kn * un' * Pn);
        w_C_RLS1 = w_C_RLS1 +kn * en;
    end
     %% RLS MCC, algorithm 2
    Pn = eye(p)*1;
    for ii = 1 : L
        Err_MCC_RLS2(mm,ii) = (wo(:,ii)  - w_C_RLS2)' * (wo(:,ii)  - w_C_RLS2);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_C_RLS2' * un;
        
        kn = Pn * un / ( exp(en^2/2/sigma2^2)*lamda2 + un' * Pn * un );
        %  kn = Pn * un / ( lamda+ un' * Pn * un );
        Pn = 1/lamda2 * ( Pn - kn * un' * Pn);
        w_C_RLS2 = w_C_RLS2 +kn * en;
    end
     %% RLS MCC, algorithm 3
    Pn = eye(p)*1;
    for ii = 1 : L
        Err_MCC_RLS3(mm,ii) = (wo(:,ii)  - w_C_RLS3)' * (wo(:,ii)  - w_C_RLS3);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_C_RLS3' * un;
        
        kn = Pn * un / ( exp(en^2/2/sigma3^2)*lamda3 + un' * Pn * un );
        %  kn = Pn * un / ( lamda3+ un' * Pn * un );
        Pn = 1/lamda3 * ( Pn - kn * un' * Pn);
        w_C_RLS3 = w_C_RLS3 +kn * en;
    end
    
    %% RLS MCC, algorithm 4
    Pn = eye(p)*1;
    for ii = 1 : L
        Err_MCC_RLS4(mm,ii) = (wo(:,ii)  - w_C_RLS4)' * (wo(:,ii)  - w_C_RLS4);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_C_RLS4' * un;
        
        kn = Pn * un / ( exp(en^2/2/sigma4^2)*lamda4 + un' * Pn * un );
        %  kn = Pn * un / ( lamda4+ un' * Pn * un );
        Pn = 1/lamda4 * ( Pn - kn * un' * Pn);
        w_C_RLS4 = w_C_RLS4 +kn * en;
    end
     %% RLS MCC, algorithm 5
    Pn = eye(p)*1;
    for ii = 1 : L
        Err_MCC_RLS5(mm,ii) = (wo(:,ii)  - w_C_RLS5)' * (wo(:,ii)  - w_C_RLS5);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_C_RLS5' * un;
        
        kn = Pn * un / ( exp(en^2/2/sigma5^2)*lamda5 + un' * Pn * un );
        %  kn = Pn * un / ( lamda5+ un' * Pn * un );
        Pn = 1/lamda5 * ( Pn - kn * un' * Pn);
        w_C_RLS5 = w_C_RLS5 +kn * en;
    end
    
     
    %% RLS-MCC, steady-state mean square performance

     vG1 = exp(-vv.^2/2/sigma1^2).*vv;
     vG2 = exp(-vv.^2/2/sigma2^2).*vv;
     vG3 = exp(-vv.^2/2/sigma3^2).*vv;
     vG4 = exp(-vv.^2/2/sigma4^2).*vv;
     vG5 = exp(-vv.^2/2/sigma5^2).*vv;
     Evg1 = mean(exp(- (vv).^2/sigma1^2/2));
     Evg2 = mean(exp(- (vv).^2/sigma2^2/2));
     Evg3 = mean(exp(- (vv).^2/sigma3^2/2));
     Evg4 = mean(exp(- (vv).^2/sigma4^2/2));
     Evg5 = mean(exp(- (vv).^2/sigma5^2/2));
     
    Err_TH_MCC1(mm) = (1-lamda1)/(1+lamda1) * p/ mean(diag(uu*uu'/L)) ...
        *  vG1*vG1'/L / Evg1^2; %1;  
    Err_TH_MCC2(mm) = (1-lamda2)/(1+lamda2) * p/ mean(diag(uu*uu'/L)) ...
        *  vG2*vG2'/L / Evg2^2; %1;
    Err_TH_MCC3(mm) = (1-lamda3)/(1+lamda3) * p/ mean(diag(uu*uu'/L)) ...
        *  vG3*vG3'/L / Evg3^2; %1;
    Err_TH_MCC4(mm) = (1-lamda4)/(1+lamda4) * p/ mean(diag(uu*uu'/L)) ...
        *  vG4*vG4'/L / Evg4^2; %1;
    Err_TH_MCC5(mm) = (1-lamda5)/(1+lamda5) * p/ mean(diag(uu*uu'/L)) ...
        *  vG5*vG5'/L / Evg5^2; %1;
    
end
toc
figure,hold on;
plot(10* log10(mean(Err_MCC_RLS1)),'r'),
plot(10* log10(mean(Err_MCC_RLS2)),'g'),
plot(10* log10(mean(Err_MCC_RLS3)),'b');
plot(10* log10(mean(Err_MCC_RLS4)),'c'),
plot(10* log10(mean(Err_MCC_RLS5)),'k');
plot(10*log10(ones(1,L)*mean(Err_TH_MCC1)),'r--')
plot(10*log10(ones(1,L)*mean(Err_TH_MCC2)),'g--')
plot(10*log10(ones(1,L)*mean(Err_TH_MCC3)),'b--')
plot(10*log10(ones(1,L)*mean(Err_TH_MCC4)),'c--')
plot(10*log10(ones(1,L)*mean(Err_TH_MCC5)),'k--')
legend('RMC1','RMC2','RMC3','RMC4','RMC5', ...
'TH-RMC1', 'TH-RMC2', 'TH-RMC3', 'TH-RMC4', 'TH-RMC5');
xlabel('Iterations');ylabel('MSD');
axis([0 2000 -40 20])
% legend('LMS','RLS','RMC','RMEE');
[10*log10(mean(Err_TH_MCC1))
10*log10(mean(Err_TH_MCC2))
10*log10(mean(Err_TH_MCC3))
10*log10(mean(Err_TH_MCC4))
10*log10(mean(Err_TH_MCC5))]'