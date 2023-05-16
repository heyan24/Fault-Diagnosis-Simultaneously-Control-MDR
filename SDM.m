tic
warning off
alpha=0.075; 
m=3000; 
pi0=1; pi1=0; pi2=1-pi0-pi1; %初始分布
B1=1/2; B=0.05; B2=1/2;

df=5;
COVMAT=eye(m);
% A=[0.75,0.125,0.125;
%     0.6,0.3,0.1;
%     0.6,0.1,0.3];

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


%for i=1:m
%	for j=1:m
%        COVMAT(i,j)=0.5^(abs(i-j));
%	end
%end
% Sum=m/10;
% for ss=1:Sum
%     for i=1:10
%         for j=1:10
%         	COVMAT((ss-1)*10+i,(ss-1)*10+j)=(i==j)+(i~=j)*0.3;
%         end
%     end
% end

%Lam=[0.28 0.3 0.32];n=1; %0.28
%Lam=[0.18 0.19 0.2 0.22];n=2; %0.22
%Lam=[0.15 0.16 0.17 0.18 0.19 0.20 0.22];n=5; %0.16
%Lam=[0.125 0.135 0.12 0.13 0.14];n=10; %0.14
%Lam=[0.105 0.10 0.11 0.12 0.13];n=15; %0.12
%Lam=[0.095 0.10 0.105];n=20; %0.105
%Lam=0.28;

%LAM=[0.24,0.22,0.145,0.11,0.1,0.09]; %矩阵1，A1=2
%LAM=[0.24,0.22,0.16,0.125,0.11,0.10]; %矩阵1，A1=3
%LAM=[0.24,0.22,0.135,0.106,0.093,0.083]; %矩阵2，A1=2
%LAM=[0.24,0.2,0.15,0.115,0.1,0.09]; %矩阵2，A1=3
%LAM=[0.24,0.22,0.143,0.111,0.097,0.087]; %矩阵3，A1=2
%LAM=[0.24,0.212,0.155,0.122,0.109,0.099]; %矩阵3，A1=3
%LAM=[0.24,0.2,0.14,0.135,0.12,0.11]; %Gamma
%LAM=[0.24,0.2,0.205,0.17,0.155,0.145]; %Gamma A1=3
%LAM=[0.22,0.2,0.15,0.115,0.098,0.087]; %t
%LAM=[0.22,0.18,0.17,0.123,0.105,0.095]; %t A1=3
%LAM=[0.24,0.22,0.157,0.123,0.11,0.103]; %相关0.3，A1=3，矩阵1
%LAM=[0.24,0.22,0.157,0.125,0.112,0.103]; %相关0.5，A1=3，矩阵1
%LAM=[0.24,0.22,0.15,0.12,0.105,0.095]; %m=1000

rep1=100; rep2=10;

LAM=[0.217,0.22,0.227,0.241,0.258,0.278,0.3]; % prop n=5 0.05
%LAM=[0.214,0.219,0.227,0.24,0.26,0.28,0.3]; % prop n=5 0.1

LAM=[0.286];

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
for t=1:rep1
    mu=zeros(1,m);
    THETA=zeros(1,m);
    THETA(1)=binornd(1,1-pi0); %下面根据A产生THETA
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
    
    MU=linspace(-5,5,200); MU1=MU(101:199); MU2=MU(1:100);
    a1=0;b1=0;a2=0;b2=0;c=0;a=0;b=0;
	g1=zeros(1,m); g2=zeros(1,m); result=zeros(1,m);
    
    for k=1:rep2
        %X=trnd(df,n,m);X=X/sqrt(df/(df-2));
        %X=gamrnd(3,1,n,m); X=(X-3)/sqrt(3);
        X=mvnrnd(mu,COVMAT,n);% observed data
        if n>1
        	X=mean(X);
        end
        %X=X+mu;
        [p0_e,p1_e,p2_e]=PEST_dire(X,0,sigma);
        %[g,~]=ksdensity(X,X);
        
        Hmu=zeros(1,length(MU)-1);
        for i=1:length(Hmu)
        	Hmu(i)=hmu(1-p0_e,n,lam,X,MU(i));
        end
        Hmu1=Hmu(101:199)*(1-p0_e)/p1_e;
        Hmu2=Hmu(1:100)*(1-p0_e)/p2_e; 
        for i=1:m
            x=X(i);
            g1(i)=sum(normpdf(x,MU1,sigma).*Hmu1/20);
            g2(i)=sum(normpdf(x,MU2,sigma).*Hmu2/20);
        end
        g0=normpdf(X,0,sigma);
        G0=p0_e*g0; G1=p1_e*g1; G2=p2_e*g2;
        G=G0+G1+G2;
        MAX=max(G1,G2);
        MAX_IND=MAX==G1;
        POS1=find(MAX_IND==1);POS2=find(MAX_IND==0);
        LAMBDA=MAX./(G0);
        [LAM_sort,Ind]=sort(LAMBDA);
        MAX=MAX./G;
        MAX_sort=MAX(Ind); SUM=sum((G1+G2)./G);
        for j=1:m
        	result(j)=sum(MAX_sort(j:m))/SUM;
        end
        cut=find(result-1+alpha<=0,1);
        OC=Ind(cut:m);
        delta=zeros(1,m);
        delta(intersect(OC,POS1))=1;delta(intersect(OC,POS2))=-1;
        
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
fprintf('A=%.2f, MDR1=%.4f, MDR2=%.4f, MDR=%.4f,EFP=%.2f \n',As(1), mean(MDR1(:,r)), mean(MDR2(:,r)),mean(MDR(:,r)),mean(EFP(:,r)))
end

SDM_result=[mean(MDR1)',mean(MDR2)',mean(MDR)',mean(EFP)'];
%csvwrite("../results/SDM2_prop2.csv",SDM_result);
toc
load train
sound(y,Fs)
