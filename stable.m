function [p0,p1,p2]=stable(A)
Q=A; 
Q(1,3)=1; Q(2,3)=1; Q(3,3)=1;
Q(1,1)=A(1,1)-1; Q(2,2)=A(2,2)-1;
C=[0;0;1];
P=linsolve(Q',C);
p0=P(1); p1=P(2); p2=P(3);
%disp([p0,p1,p2]);
end

