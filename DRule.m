% This function calculates the theoretical decision rule given lambda and H
% that minimizes the loss function.

function [delta, R1, R2] = DRule(lambda, H0, H1, H2)
m = length(H0);
delta = zeros(1,m);
R1 = H0 - lambda(1)*H1;
R2 = H0 - lambda(2)*H2;
delta(R1<=0 & R1<=R2) = 1;
delta(R2<=0 & R2<=R1) = -1;
%delta(R1<=0) = 1;
end

