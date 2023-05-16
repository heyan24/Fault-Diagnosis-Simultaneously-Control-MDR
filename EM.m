%Int为pi与A的初始值，g为已经得到的X的条件分布值
function [pi0_e,pi1_e,pi2_e,A_e]=EM(m,g0,g1,g2,Int)
pi0_e=Int(1,1);pi1_e=Int(1,2);pi2_e=Int(1,3);
A_e=[Int(2,1),Int(2,2),Int(2,3);
    Int(3,1),Int(3,2),Int(3,3);
    Int(4,1),Int(4,2),Int(4,3)];

%C=1;
l=0;
limit=100;
%err=1e-30;
while l<limit
%while n<=10000
    %E-step
    c0=zeros(1,m);
    d0=zeros(1,m);
    Alpha0=zeros(1,m); Alpha1=zeros(1,m); Alpha2=zeros(1,m);
    Beta0=zeros(1,m); Beta1=zeros(1,m); Beta2=zeros(1,m);
    Alpha0(1)=pi0_e*g0(1); Alpha1(1)=pi1_e*g1(1); Alpha2(1)=pi2_e*g2(1);
    c0(1)=1/(Alpha0(1)+Alpha1(1)+Alpha2(1));
    Alpha0(1)=Alpha0(1)*c0(1); Alpha1(1)=Alpha1(1)*c0(1); Alpha2(1)=Alpha2(1)*c0(1);  
    Beta0(m)=1; Beta1(m)=1; Beta2(m)=1;
    d0(m)=1/(Beta0(m)+Beta1(m)+Beta2(m));
    for i=2:m
       Alpha0(i)=g0(i)*(Alpha0(i-1)*A_e(1,1)+Alpha1(i-1)*A_e(2,1)+Alpha2(i-1)*A_e(3,1));
       Alpha1(i)=g1(i)*(Alpha0(i-1)*A_e(1,2)+Alpha1(i-1)*A_e(2,2)+Alpha2(i-1)*A_e(3,2));
       Alpha2(i)=g2(i)*(Alpha0(i-1)*A_e(1,3)+Alpha1(i-1)*A_e(2,3)+Alpha2(i-1)*A_e(3,3));
       c0(i)=1/(Alpha0(i)+Alpha1(i)+Alpha2(i));
       Alpha0(i)=Alpha0(i)*c0(i); Alpha1(i)=Alpha1(i)*c0(i); Alpha2(i)=Alpha2(i)*c0(i); 
    end
    Beta0(m)=Beta0(m)*d0(m); Beta1(m)=Beta1(m)*d0(m); Beta2(m)=Beta2(m)*d0(m); 
    %Beta0(m)=Beta0(m); Beta1(m)=Beta1(m); Beta2(m)=Beta2(m); 
    for i=2:m
       Beta0(m-i+1)=(A_e(1,1)*g0(m-i+2)*Beta0(m-i+2)+A_e(1,2)*g1(m-i+2)*Beta1(m-i+2)+A_e(1,3)*g2(m-i+2)*Beta2(m-i+2));%*c0(m-i+1);
       Beta1(m-i+1)=(A_e(2,1)*g0(m-i+2)*Beta0(m-i+2)+A_e(2,2)*g1(m-i+2)*Beta1(m-i+2)+A_e(2,3)*g2(m-i+2)*Beta2(m-i+2));%*c0(m-i+1);
       Beta2(m-i+1)=(A_e(3,1)*g0(m-i+2)*Beta0(m-i+2)+A_e(3,2)*g1(m-i+2)*Beta1(m-i+2)+A_e(3,3)*g2(m-i+2)*Beta2(m-i+2));%*c0(m-i+1);
       d0(m-i+1)=1/(Beta0(m-i+1)+Beta1(m-i+1)+Beta2(m-i+1));
       Beta0(m-i+1)=Beta0(m-i+1)*d0(m-i+1);Beta1(m-i+1)=Beta1(m-i+1)*d0(m-i+1);Beta2(m-i+1)=Beta2(m-i+1)*d0(m-i+1);
    end
    
    %M-step
    r0=Alpha0.*Beta0; r1=Alpha1.*Beta1; r2=Alpha2.*Beta2;
    s=r0+r1+r2;
    r0=r0./s; r1=r1./s; r2=r2./s;
    k00=zeros(0,m-1); k01=zeros(0,m-1); k02=zeros(0,m-1);
    k10=zeros(0,m-1); k11=zeros(0,m-1); k12=zeros(0,m-1);
    k20=zeros(0,m-1); k21=zeros(0,m-1); k22=zeros(0,m-1);
    for i=1:(m-1)
        k00(i)=A_e(1,1)*r0(i)*g0(i+1)*Beta0(i+1)*d0(i)/Beta0(i);
        k01(i)=A_e(1,2)*r0(i)*g1(i+1)*Beta1(i+1)*d0(i)/Beta0(i);
        k02(i)=A_e(1,3)*r0(i)*g2(i+1)*Beta2(i+1)*d0(i)/Beta0(i);
        k10(i)=A_e(2,1)*r1(i)*g0(i+1)*Beta0(i+1)*d0(i)/Beta1(i);
        k11(i)=A_e(2,2)*r1(i)*g1(i+1)*Beta1(i+1)*d0(i)/Beta1(i);
        k12(i)=A_e(2,3)*r1(i)*g2(i+1)*Beta2(i+1)*d0(i)/Beta1(i);
        k20(i)=A_e(3,1)*r2(i)*g0(i+1)*Beta0(i+1)*d0(i)/Beta2(i);
        k21(i)=A_e(3,2)*r2(i)*g1(i+1)*Beta1(i+1)*d0(i)/Beta2(i);
        k22(i)=A_e(3,3)*r2(i)*g2(i+1)*Beta2(i+1)*d0(i)/Beta2(i);
    end
    pi0_new=r0(1); pi1_new=r1(1); pi2_new=r2(1);
    A_new=[sum(k00)/sum(r0(1:m-1)),sum(k01)/sum(r0(1:m-1)),sum(k02)/sum(r0(1:m-1));
        sum(k10)/sum(r1(1:m-1)),sum(k11)/sum(r1(1:m-1)),sum(k12)/sum(r1(1:m-1));
        sum(k20)/sum(r2(1:m-1)),sum(k21)/sum(r2(1:m-1)),sum(k22)/sum(r2(1:m-1))];

    %judge=[abs(pi0_new-pi0_e)<err,abs(pi1_new-pi1_e)<err,abs(pi2_new-pi2_e)<err;
    %    abs(A_new-A_e)<repmat(err,3,3)];
%     if sum(sum(judge))==12
%        pi0_e=pi0_new; pi1_e=pi1_new; pi2_e=pi2_new;
%        reg1=sum(A_new(1,:)); reg2=sum(A_new(2,:)); reg3=sum(A_new(3,:));
%        A_new(1,:)=A_new(1,:)/reg1; 
%        A_new(2,:)=A_new(2,:)/reg2;
%        A_new(3,:)=A_new(3,:)/reg3;
%        A_e=A_new;
%        break
%     end
    pi0_e=pi0_new; pi1_e=pi1_new; pi2_e=pi2_new;
%        reg1=sum(A_new(1,:)); reg2=sum(A_new(2,:)); reg3=sum(A_new(3,:));
%        A_new(1,:)=A_new(1,:)/reg1; 
%        A_new(2,:)=A_new(2,:)/reg2;
%        A_new(3,:)=A_new(3,:)/reg3;
    A_e=A_new;
    l=l+1;
    %c=c+1;
%     if c>=limit
%        reg1=sum(A_new(1,:)); reg2=sum(A_new(2,:)); reg3=sum(A_new(3,:));
%        A_new(1,:)=A_new(1,:)/reg1; 
%        A_new(2,:)=A_new(2,:)/reg2;
%        A_new(3,:)=A_new(3,:)/reg3;
%        A_e=A_new;
%        break
%     end
end

  

end