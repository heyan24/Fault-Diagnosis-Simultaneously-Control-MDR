tic
warning off
alpha1 = 0.1; alpha2 = 0.1;
m=3000;  %n=1;
pi0=1; pi1=0; pi2=1-pi0-pi1; %初始分布
B1=1/2; B=0.05; B2=1/2; %mean of gamma: A1*B1
df = 5;
COVMAT=eye(m); 
% A=[1/3,1/3,1/3;
%     1/3,1/3,1/3;
%     1/3,1/3,1/3]; %independent
% A=[0.6,0.2,0.2;
%     0.3,0.4,0.3;
%     0.3,0.3,0.4];  %weak correlation
% A=[0.6,0.2,0.2;
%     0.25,0.5,0.25;
%     0.25,0.25,0.5];  %weak correlation
% A=[0.8,0.1,0.1;
%     0.15,0.7,0.15;
%     0.15,0.15,0.7]; %strong correlation
% A=[0.7,0.2,0.1;
%     0.25,0.5,0.25;
%     0.35,0.35,0.3]; %不对称情况 tmMDR = 0.082

A_list = zeros(3,3,7);
%p_list=[0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975]; %null proportions = 0.539, 0.6, 0.667, 0.739, 0.818, 0.905, 0.9512
p_list=[0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975]; %null proportions = 0.7021, 0.7391, 0.7778, 0.8182, 0.8605, 0.905, 0.9512
p_list=[0.975];
for i=1:length(p_list)
    p = p_list(i);
    A_list(:,:,i) = [p,(1-p)/2, (1-p)/2;
                    p/2, 0.5, 0.5-p/2;
                    p/2, 0.5-p/2, 0.5];
end

% for i=1:m
% 	for j=1:m
%         COVMAT(i,j)=0.5^(abs(i-j));
% 	end
% end
% Sum=m/10;
% for ss=1:Sum
%     for i=1:10
%         for j=1:1
%         	COVMAT((ss-1)*10+i,(ss-1)*10+j)=(i==j)+(i~=j)*0.3;
%         end
%     end
% end 

%LAM=[0.16,0.16,0.232,0.255,0.275,0.3,0.35]; %ind, n=2, 0.05
%LAM=[0.15,0.18,0.243,0.255,0.273,0.31,0.36]; %ind, n=2, 0.1
%LAM=[0.15,0.17,0.244,0.256,0.273,0.32,0.36]; %weak, n=2, 0.05
%LAM=[0.16,0.17,0.245,0.26,0.275,0.32,0.36]; %strong, n=2, 0.05
%LAM=[0.19,0.245,0.25,0.265,0.28,0.32,0.36]; %strong, n=2, 0.1
%LAM=[0.19,0.21,0.25,0.265,0.27,0.31,0.35]; %asym, n=2, 0.05
%LAM=[0.15,0.165,0.25,0.267,0.28,0.32,0.35]; %asym, n=2, 1

%LAM=[0.144,0.149,0.159,0.165,0.17,0.172,0.18]; %diff1 A, n=5, 0.05 weak
%LAM=[0.144,0.154,0.163,0.175,0.182,0.191,0.198]; %diff2 A, n=5, 0.05 weak
%LAM=[0.15,0.158,0.165,0.17,0.176,0.183,0.194]; %diff1 A, n=5, 0.05 strong
%LAM=[0.15,0.16,0.166,0.174,0.182,0.204,0.207]; %diff2 A, n=5, 0.05 strong

%LAM=[0.14,0.154,0.163,0.173,0.183,0.192,0.205]; %t, n=5, 0.05 weak
%LAM=[0.155,0.173,0.19,0.225,0.268,0.275,0.36]; %gam, n=5, 0.05 weak
%LAM=[0.145,0.157,0.165,0.175,0.185,0.197,0.21]; %t, n=5, 0.05 strong
%LAM=[0.158,0.178,0.195,0.23,0.27,0.28,0.36]; %gam, n=5, 0.05 strong

%LAM=[0.21,0.18,0.137,0.11,0.095,0.084]; %varies n, A=2, 0.05 weak
%LAM=[0.21,0.155,0.147,0.112,0.095,0.09]; %varies n, A=2, 0.1 weak

%LAM=[0.14,0.155,0.168,0.195,0.222,0.26,0.305]; %corr rho=0.3 0.05
%LAM=[0.137,0.155,0.168,0.193,0.22,0.26,0.307]; %corr rho=0.5 0.05 weak
%LAM=[0.144,0.158,0.172,0.193,0.22,0.26,0.307]; %corr rho=0.5 0.1 weak
%LAM=[0.145,0.165,0.173,0.205,0.23,0.272,0.315]; %corr rho=0.5 0.05 strong
%LAM=[0.155,0.163,0.175,0.195,0.23,0.27,0.31]; %corr rho=0.5 0.1 strong

%LAM=[0.222,0.24,0.257,0.265,0.271,0.272,0.273]; %prop n=5 0.05
LAM=[0.224,0.24,0.257,0.268,0.32,0.272,0.273]; %prop n=5 0.1

LAM=[0.33];
rep1=10; rep2=10; 
%N=[1,2,5,10,15,20];
%N=[1];
n = 5;
sigma=1/sqrt(n); 
%As = 2:0.5:5;
As = [4];
%MDR1=zeros(rep1,length(As)); MDR2=zeros(rep1,length(As)); MDR=zeros(rep1,length(As)); EFP=zeros(rep1,length(As));
MDR1=zeros(rep1,length(p_list)); MDR2=zeros(rep1,length(p_list)); MDR=zeros(rep1,length(p_list)); EFP=zeros(rep1,length(p_list));
%MDR1=zeros(rep1,length(N)); MDR2=zeros(rep1,length(N)); MDR=zeros(rep1,length(N)); EFP=zeros(rep1,length(N));
for r = 1:length(p_list)
lam=LAM(r);
A1 = As(1); A2 = As(1);
A = A_list(:,:,r);
[p0,p1,p2] = stable(A);
%n = N(r);
%sigma=1/sqrt(n); 
Steps = [];
%THETA的产生需要依据HMM模型来产生
for t=1:rep1
    c0=zeros(1,m);
    mu=zeros(1,m);
    THETA=zeros(1,m);
    THETA(1)=1-binornd(1,pi0); %下面根据A产生THETA
    if THETA(1)==1
        ind=binornd(1,pi1/(pi1+pi2));
        THETA(1)=ind*1+(1-ind)*-1;
        mu(1)=ind*(gamrnd(A1,B1)+B)+(1-ind)*(-gamrnd(A2,B2)-B);
    end
    for i=2:m
        if THETA(i-1)==0
            THETA(i)=binornd(1,1-A(1,1));
            if THETA(i)==1
                ind=binornd(1,A(1,2)/(A(1,2)+A(1,3)));
                THETA(i)=ind*1+(1-ind)*-1;
                mu(i)=ind*(gamrnd(A1,B1)+B)+(1-ind)*(-gamrnd(A2,B2)-B);
            end    
        elseif THETA(i-1)==1
            THETA(i)=binornd(1,1-A(2,1));
            if THETA(i)==1
                ind=binornd(1,A(2,2)/(A(2,2)+A(2,3)));
                THETA(i)=ind*1+(1-ind)*-1;
                mu(i)=ind*(gamrnd(A1,B1)+B)+(1-ind)*(-gamrnd(A2,B2)-B);
            end
        elseif THETA(i-1)==-1
            THETA(i)=binornd(1,1-A(3,1));
            if THETA(i)==1
                ind=binornd(1,A(3,2)/(A(3,2)+A(3,3)));
                THETA(i)=ind*1+(1-ind)*-1;
                mu(i)=ind*(gamrnd(A1,B1)+B)+(1-ind)*(-gamrnd(A2,B2)-B);
            end
        end
    end
%     THETA=binornd(1,1-p0,1,m);
%     mu=zeros(1,m);
%     for i=1:m
%         if THETA(i)==1
% 		ind=binornd(1,p1/(p1+p2));
% 		THETA(i)=ind*1+(1-ind)*-1;
% 		mu(i)=ind*(gamrnd(A1,B1)+B)+(1-ind)*(-gamrnd(A2,B2)-B);
%         end
%     end
    %以上产生了THETA与mu的模拟值
	MU=linspace(-5,5,200); MU1=MU(101:199); MU2=MU(1:100);
	a1=0;b1=0;a2=0;b2=0;c=0;a=0; b=0;
	g1=zeros(1,m);g2=zeros(1,m);result=zeros(1,m);

	for k=1:rep2
        %X=gamrnd(3,1,n,m); X=(X-3)/sqrt(3);
        %X=trnd(df,n,m);X=X/sqrt(df/(df-2));
	    X=mvnrnd(mu,COVMAT,n); %产生X的模拟值
        if n>1
            X=mean(X); 
        end
        %X=X+mu;
        [p0_e,p1_e,p2_e]=PEST_dire(X,0,sigma);
        
        Hmu=zeros(1,length(MU)-1);
        for i=1:length(Hmu)
        	Hmu(i)=hmu(1-p0_e,n,lam,X,MU(i));
        end
        Hmu1=Hmu(101:199)*(1-p0_e)/p1_e;
        Hmu2=Hmu(1:100)*(1-p0_e)/p2_e; 
        for i=1:m
            x=X(i);
            %g1_r(i)=sum(normpdf(x,MU1,sigma).*gampdf(MU1,A1,B1)/20);
            %g2_r(i)=sum(normpdf(x,MU2,sigma).*gampdf(-MU2,A2,B2)/20);
            g1(i)=sum(normpdf(x,MU1,sigma).*Hmu1/20);
            g2(i)=sum(normpdf(x,MU2,sigma).*Hmu2/20);
        end
        g0=normpdf(X,0,sigma);%g1=normpdf(X,mu1,sigma);g2=normpdf(X,mu2,sigma);
        
        %下面用EM算法估计pi与A
        Int=[0.8,0.1,0.1;
            0.8,0.1,0.1;
            0.7,0.2,0.1;
            0.7,0.1,0.2];
        [pi0_e,pi1_e,pi2_e,A_e]=EM(m,g0,g1,g2,Int);
         
        %到这里的计算应该都是有用的,g_k是g(X_i|theta_i=k)；
        %后面新的关于H_i^k(X)的计算依据g_k，通过f-b算法完成
        
        %下面利用forward-backward算法计算H
        
        Alpha0=zeros(1,m); Alpha1=zeros(1,m); Alpha2=zeros(1,m);
        Beta0=zeros(1,m); Beta1=zeros(1,m); Beta2=zeros(1,m);
        Alpha0(1)=pi0_e*g0(1); Alpha1(1)=pi1_e*g1(1); Alpha2(1)=pi2_e*g2(1);
        c0(1)=1/(Alpha0(1)+Alpha1(1)+Alpha2(1));
        Alpha0(1)=Alpha0(1)*c0(1); Alpha1(1)=Alpha1(1)*c0(1); Alpha2(1)=Alpha2(1)*c0(1);  
        Beta0(m)=1; Beta1(m)=1; Beta2(m)=1;
        for i=2:m
           Alpha0(i)=g0(i)*(Alpha0(i-1)*A_e(1,1)+Alpha1(i-1)*A_e(2,1)+Alpha2(i-1)*A_e(3,1));
           Alpha1(i)=g1(i)*(Alpha0(i-1)*A_e(1,2)+Alpha1(i-1)*A_e(2,2)+Alpha2(i-1)*A_e(3,2));
           Alpha2(i)=g2(i)*(Alpha0(i-1)*A_e(1,3)+Alpha1(i-1)*A_e(2,3)+Alpha2(i-1)*A_e(3,3));
           c0(i)=1/(Alpha0(i)+Alpha1(i)+Alpha2(i));
           Alpha0(i)=Alpha0(i)*c0(i); Alpha1(i)=Alpha1(i)*c0(i); Alpha2(i)=Alpha2(i)*c0(i); 
        end
        Beta0(m)=Beta0(m)*c0(m); Beta1(m)=Beta1(m)*c0(m); Beta2(m)=Beta2(m)*c0(m); 
        for i=2:m
           Beta0(m-i+1)=(A_e(1,1)*g0(m-i+2)*Beta0(m-i+2)+A_e(1,2)*g1(m-i+2)*Beta1(m-i+2)+A_e(1,3)*g2(m-i+2)*Beta2(m-i+2))*c0(m-i+1);
           Beta1(m-i+1)=(A_e(2,1)*g0(m-i+2)*Beta0(m-i+2)+A_e(2,2)*g1(m-i+2)*Beta1(m-i+2)+A_e(2,3)*g2(m-i+2)*Beta2(m-i+2))*c0(m-i+1);
           Beta2(m-i+1)=(A_e(3,1)*g0(m-i+2)*Beta0(m-i+2)+A_e(3,2)*g1(m-i+2)*Beta1(m-i+2)+A_e(3,3)*g2(m-i+2)*Beta2(m-i+2))*c0(m-i+1);
        end
        H0=Alpha0.*Beta0;
        H1=Alpha1.*Beta1;
        H2=Alpha2.*Beta2;
        H=H0+H1+H2;
        
        H0=H0./H; H1=H1./H; H2=H2./H;
        
        LAMBDA1 = H1./H0; LAMBDA2 = H2./H0;
        [LAM1_sort,Ind1]=sort(LAMBDA1); [LAM2_sort,Ind2]=sort(LAMBDA2);
        H1_sort = H1(Ind1); H2_sort = H2(Ind2);
        result1 = zeros(1,m); result2 = zeros(1,m); 
        for j=1:m
            result1(j)=sum(H1_sort(j:m))/sum(H1);
            result2(j)=sum(H2_sort(j:m))/sum(H2);
        end
        R1 = []; R2 = [];
        r1 = find(result1-1+alpha1<=0,1) -1; R1 = horzcat(R1, r1); 
        r2 = find(result2-1+alpha2<=0,1) -1; R2 = horzcat(R2, r2);
        lam1 = 1/LAM1_sort(r1); lam2 = 1/LAM2_sort(r2);

        N_lam1 = Nhat(1, DRule([lam1, lam2], H0, H1, H2), H1, alpha1);
        N_lam2 = Nhat(-1, DRule([lam1, lam2], H0, H1, H2), H2, alpha2);
        N_tilda1 = Nhat(1, DRule([1/LAM1_sort(r1-1), lam2], H0, H1, H2), H1, alpha1);
        N_tilda2 = Nhat(-1, DRule([lam1, 1/LAM2_sort(r2-1)], H0, H1, H2), H2, alpha2);
        step = 1;
 
        if (N_lam1>=0 && N_lam2>=0)
            mark = 1;
            while ~(N_lam1>=0 && N_lam2>=0 && N_tilda1<0 && N_tilda2<0)
                step = step + 1;
                N_pot1 = zeros(1,r1); N_pot2 = zeros(1,r2);
                for j=1:r1
                   N_pot1(j) = Nhat(1, DRule([1/LAM1_sort(j),lam2],H0,H1,H2), H1, alpha1);
                end
                for j=1:r2
                   N_pot2(j) = Nhat(-1, DRule([lam1,1/LAM2_sort(j)], H0,H1,H2), H2, alpha2);
                end
                r1 = find(N_pot1>0,1) ; r2 = find(N_pot2>0,1) ;
                R1 = horzcat(R1, r1); R2 = horzcat(R2, r2);
                lam1 = 1/LAM1_sort(r1); lam2 = 1/LAM2_sort(r2);
                
                if (r1==1 || r2==1)
                   break 
                end

                N_lam1 = Nhat(1, DRule([lam1, lam2], H0, H1, H2), H1, alpha1);
                N_lam2 = Nhat(-1, DRule([lam1, lam2], H0, H1, H2), H2, alpha2);
                N_tilda1 = Nhat(1, DRule([1/LAM1_sort(r1-1), lam2], H0, H1, H2), H1, alpha1);
                N_tilda2 = Nhat(-1, DRule([lam1, 1/LAM2_sort(r2-1)], H0, H1, H2), H2, alpha2);

            end
        elseif (N_lam1*N_lam2 < 0)
            mark = 0;
            if N_lam1>0 
                r1 = r1 - 1;
            elseif N_lam2>0 
                r2 = r2 - 1;
            end
        else 
            mark = 2;
        end

        Steps = horzcat(Steps, step);
        
        lam_star = [lam1, lam2];
        [delta, R11, R22] = DRule(lam_star, H0, H1, H2);
        
        a=a+sum(abs(THETA).*(THETA~=delta));
        b=b+sum(abs(THETA));
        a1=a1+sum((THETA == 1).*(THETA~=delta));
        a2=a2+sum((THETA == -1).*(THETA~=delta));
        b1=b1+sum(THETA == 1); b2=b2+sum(THETA == -1); 
        c=c+sum(abs(delta).*(1-abs(THETA)));
    end
    MDR(t,r)=a/b;
    MDR1(t,r)=a1/b1;
    MDR2(t,r)=a2/b2;
    EFP(t,r)=c/rep2;
end
fprintf('A=%.2f, Step=%.2f, MDR1=%.4f, MDR2=%.4f, MDR=%.4f,EFP=%.2f \n',As(1), mean(Steps), mean(MDR1(:,r)), mean(MDR2(:,r)),mean(MDR(:,r)),mean(EFP(:,r)))
end
%fprintf('n=%d, MDR1=%.4f, MDR2=%.4f, EFP=%.2f',n, mean(MDR1), mean(MDR2), mean(EFP))
Datadriven_result=[mean(MDR1)',mean(MDR2)',mean(MDR)',mean(EFP)'];
%csvwrite("../results/Datadriven2_prop2.csv", Datadriven_result);
toc
%end
load train
sound(y,Fs)

