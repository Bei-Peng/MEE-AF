clear all; close all;
%% MEE-AF-V10
L = 2000;
p = 5;
q = 8;
lamda2 = 0.996;%RMEE
lamda3 = 0.992;%RLS
% lamda2 = 1;%RMEE for speech signal
lamda1 = 0.942;%RMC
sigma1= 1;%%RMC
sigma2 = sigma1;%%RMEE
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
    % wo = [ kron(randn(p,1), ones(1,L/2)),  kron(randn(p,1), ones(1,L/2)) ];
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
    w_RLS = w_LMS;
    w_C_RLS =w_LMS;  %%MCC
    w_M_RLS =w_LMS;%%MEE
    %% LMS
    mu1 = 0.025;
    for ii = 1 : L
        Err_LMS(mm,ii) = (wo(:,ii) - w_LMS)' *  (wo(:,ii) - w_LMS);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_LMS' * un;
        w_LMS = w_LMS + mu1 * (en) * un;
        %    w_LMS = w_LMS + mu1 * tanh(en) * un;
    end
    
    %% RLS_MSE
    Pn = eye(p)*1;
    for ii = 1 : L
        Err_RLS(mm,ii) = (wo(:,ii)  - w_RLS)' * (wo(:,ii)  - w_RLS);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_RLS' * un;
        
       % kn = Pn * un / ( exp(en^2/2/sigma1^2)*lamda1 + un' * Pn * un );
          kn = Pn * un / ( lamda3+ un' * Pn * un );
        Pn = 1/lamda3 * ( Pn - kn * un' * Pn);
        w_RLS = w_RLS +kn * en;
    end
    
 %% RLS-MCC
    Pn = eye(p)*1;
    for ii = 1 : L
        Err_MCC_RLS(mm,ii) = (wo(:,ii)  - w_C_RLS)' * (wo(:,ii)  - w_C_RLS);
        dn = dd(ii);
        un = uu(:,ii);
        en = dn - w_C_RLS' * un;
        
        kn = Pn * un / ( exp(en^2/2/sigma1^2)*lamda1 + un' * Pn * un );
        %  kn = Pn * un / ( lamda+ un' * Pn * un );
        Pn = 1/lamda1 * ( Pn - kn * un' * Pn);
        w_C_RLS = w_C_RLS +kn * en;
    end
    
    %% RLS MEE, algorithm
    PL = eye(p)*1; 
    eL = zeros(q, 1); 
    PHI = zeros(q, q);
    for ii = q : L
        Err_RLS_MEE(mm,ii) = (wo(:,ii)  - w_M_RLS)' * (wo(:,ii)  - w_M_RLS);
        for jj = 1 : q
            eL(jj) = dd(ii - q + jj) - w_M_RLS' * uu(: , ii - q + jj);
        end
        u0 = uu(:,ii);
        e0 = eL(q);
        phi0 = 0;
        for kk = 2:q
            ek = eL(q-kk+1);
%             ek = eL(kk);
            phi0 = phi0 + lamda2^(kk-1) * exp(- (e0-ek).^2/sigma2^2/2);
        end
%         phi0 = 6.9;
        KL = PL * u0 * inv ( lamda2^2 / phi0 + u0' * PL * u0 );
        PL = 1/lamda2^2 * (PL - KL * u0' * PL);
        w_M_RLS = w_M_RLS + KL * e0;
        
    end
     
    %% RLS-MEE, steady-state mean square performance
    Gsig =0; %phi0
    phi0 = 0;
    v0 = vv(q: end);
    for ii = 2:q
        phi0 = phi0 + lamda2^(ii-1);
    end
    for ii = 2:q
        vi = vv(q-ii+1 : end-ii+1);
         Gsig = Gsig + lamda2^(ii-1) * exp(- (v0-vi).^2/sigma2^2/2);
%        Gsig = Gsig + lamda2^(q - ii +1) * exp(- (v0-vi).^2/sigma2^2/2);
    end
    tmp0 = (v0.*Gsig) ;
    L1 = length(tmp0);
    Evb = [tmp0  ];
    %Evbbv = mean(diag(Evb*Evb'/L1));
    Evbbv =Evb*Evb'/L1;
    Evg = mean(exp(- (vv).^2/sigma2^2/2));
    % a =(1-lamda^2)/(1+lamda^2) * p / mean(diag(uu*uu'/L)) ...
    %        * mean(tmp1) / mean(Gsig.^2)
    a =(1-lamda2^2)/(1+lamda2^2) * p / mean(diag(uu*uu'/L)) ...
        * Evbbv / mean(Gsig)^2; %phi0^2; %
    Err_TH_MEE(mm) = a;
    Err_TH_MCC(mm) = (1-lamda1)/(1+lamda1) * p/ mean(diag(uu*uu'/L)) ...
        *  vG*vG'/L / Evg^2; %1;
     Err_TH_MSE(mm) = (1-lamda3)/(1+lamda3) * p/ mean(diag(uu*uu'/L)) ...
        *  vv*vv'/L / 1;
end
toc
figure,hold on;
plot(10* log10(mean(Err_LMS)),'g'),
plot(10* log10(mean(Err_RLS)),'b'),
plot(10* log10(mean(Err_MCC_RLS)),'r');
plot(10* log10(mean(Err_RLS_MEE)),'k');
plot(10*log10(ones(1,L)*mean(Err_TH_MSE)),'b--')
plot(10*log10(ones(1,L)*mean(Err_TH_MCC)),'r--')
plot(10*log10(ones(1,L)*mean(Err_TH_MEE)),'k--')
legend('LMS','RLS','RMC','RMEE','TH-RLS','TH-RMC','TH-RMEE');
% legend('LMS','RLS','RMC','RMEE');

xlabel('Iterations');ylabel('MSD');
figure,hold on;
plot((mean(Err_LMS)),'g'),
plot((mean(Err_RLS)),'b'),
plot((mean(Err_MCC_RLS)),'r');
plot((mean(Err_RLS_MEE)),'k');
plot((ones(1,L)*mean(Err_TH_MEE)),'k')
plot((ones(1,L)*mean(Err_TH_MCC)),'r')
plot((ones(1,L)*mean(Err_TH_MSE)),'b')
legend('LMS','RLS','RMC','RMEE');
xlabel('Iterations');ylabel('MSD');

figure,plot(v0),
hold on; plot(v0.*Gsig/ phi0^2,'r*')

% figure,
% subplot(211),plot(uu(1,:));title('Speech signal');
% axis([0 3000 -3 3]);
% subplot(212), plot(vv);title('Impulive noise');axis([0 3000 -3 3]);
