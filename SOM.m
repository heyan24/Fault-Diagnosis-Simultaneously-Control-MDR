tic
warning off
alpha=0.1; 
m=3000;
pi0=1; pi1=0; pi2=1-pi0-pi1; %初始分布
B1=1/2; B=0.05; B2=1/2; %mean of gamma: A1*B1
df=5;
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
%     0.35,0.35,0.3]; %不对称情况 19/30 vs. 11/30 tmMDR = 0.082

A_list = zeros(3,3,7);
%p_list=[0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975]; %null proportions = 0.539, 0.6, 0.667, 0.739, 0.818, 0.905, 0.9512
p_list=[0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975]; %null proportions = 0.6667, 0.7021, 0.7391, 0.7778, 0.8182, 0.8605, 0.905, 0.9512
%p_list=[0.825];
for i=1:length(p_list)
    p = p_list(i);
    A_list(:,:,i) = [p,(1-p)/2, (1-p)/2;
                    p/2, 0.5, 0.5-p/2;
                    p/2, 0.5-p/2, 0.5];
end

%[p0,p1,p2]=stable(A);
rep1=10; rep2=1; 
% for i=1:m
% 	for j=1:m
%        COVMAT(i,j)=0.5^(abs(i-j));
% 	end
% end

%N=[1,2,5,10,15,20];
n = 10;
sigma=1/sqrt(n);
%As = 2:0.5:5;
As = [4];
%MDR1=zeros(rep1,length(As)); MDR2=zeros(rep1,length(As)); MDR=zeros(rep1,length(As)); EFP=zeros(rep1,length(As));
MDR1=zeros(rep1,length(p_list)); MDR2=zeros(rep1,length(p_list)); MDR=zeros(rep1,length(p_list)); EFP=zeros(rep1,length(p_list));
%MDR1=zeros(rep1,length(N)); MDR2=zeros(rep1,length(N)); MDR=zeros(rep1,length(N)); EFP=zeros(rep1,length(N));
for r = 1:length(p_list)
A1 = As(1); A2 = As(1);
A = A_list(:,:,r);
[p0,p1,p2] = stable(A);
%n = N(r);
%sigma=1/sqrt(n); 
for t=1:rep1
    c0=zeros(1,m);
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
    
    MU1=linspace(B,6+B,120);MU2=linspace(-6-B,-B,120);
    a1=0;b1=0;a2=0;b2=0;c=0;a=0;b=0;
	g1=zeros(1,m);g2=zeros(1,m);result=zeros(1,m);
    for k=1:rep2
        %X=gamrnd(3,1,n,m); X=(X-3)/sqrt(3);
        %X=trnd(df,n,m);X=X/sqrt(df/(df-2));
		X=mvnrnd(mu,COVMAT,n);% observed data
        if n>1
        	X=mean(X);
        end
        %X=X+mu;
        for i=1:m
			x=X(i);
            g1(i)=sum(normpdf(x,MU1,sigma).*gampdf(MU1-B,A1,B1)/20);
            g2(i)=sum(normpdf(x,MU2,sigma).*gampdf(-MU2-B,A2,B2)/20);
        end
        g0=normpdf(X,0,sigma);%g1=normpdf(X,mu1,sigma);g2=normpdf(X,mu2,sigma);
        G0=p0*g0; G1=p1*g1; G2=p2*g2;
        G=G0+G1+G2;
        MAX=max(G1,G2);
        MAX_IND=MAX==G1;
        POS1=find(MAX_IND==1);POS2=find(MAX_IND==0);
        LAMBDA=MAX./(G0);
        [LAM_sort,Ind]=sort(LAMBDA);
        MAX=MAX./G;
        MAX_sort=MAX(Ind);
        SUM=sum((G1+G2)./G);
        
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
SOM_result=[mean(MDR1)',mean(MDR2)',mean(MDR)',mean(EFP)'];
%csvwrite("../results/SOM2_prop2.csv",SOM_result);
%end
toc
load train
sound(y,Fs)