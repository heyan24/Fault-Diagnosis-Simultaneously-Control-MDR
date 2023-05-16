% This function implements the algorithm that get lambda_star, 
% the estiamtes of the optimal lambda.

function lam_star = get_lamstar(H1, H2, H0, alpha1, alpha2)
m = length(H0);
LAMBDA1 = H1./H0; LAMBDA2 = H2./H0;
[LAM1_sort,Ind1]=sort(LAMBDA1); [LAM2_sort,Ind2]=sort(LAMBDA2);
H1_sort = H1(Ind1); H2_sort = H2(Ind2);
result1 = zeros(1,m); result2 = zeros(1,m); 
for j=1:m
    result1(j)=sum(H1_sort(j:m))/sum(H1);
    result2(j)=sum(H2_sort(j:m))/sum(H2);
end
r1 = find(result1-1+alpha1<=0,1) - 1;
r2 = find(result2-1+alpha2<=0,1) - 1;
lam1 = 1/LAM1_sort(r1); lam2 = 1/LAM2_sort(r2);

N_lam1 = Nhat(1, DRule([lam1, lam2], H0, H1, H2), H1, alpha1);
N_lam2 = Nhat(-1, DRule([lam1, lam2], H0, H1, H2), H2, alpha2);
N_tilda1 = Nhat(1, DRule([1/LAM1_sort(r1+1), lam2], H0, H1, H2), H1, alpha1);
N_tilda2 = Nhat(-1, DRule([lam1, 1/LAM2_sort(r2+1)], H0, H1, H2), H2, alpha2);

while ~(N_lam1<=0 && N_lam2<=0 && N_tilda1>0 && N_tilda2>0)
    N_pot1 = zeros(1,m-r1+1); N_pot2 = zeros(1,m-r2+1);
    for j=r1:m
       N_pot1(j-r1+1) = Nhat(1, DRule([1/LAM1_sort(j),lam2], H0,H1,H2), H1, alpha1);
    end
    for j=r2:m
       N_pot2(j-r2+1) = Nhat(-1, DRule([lam1,1/LAM2_sort(j)], H0,H1,H2), H2, alpha2);
    end
    r1 = find(N_pot1>=0,1) - 1; r2 = find(N_pot2>=0,1) - 1;
    lam1 = 1/LAM1_sort(r1); lam2 = 1/LAM2_sort(r2);
    
    N_lam1 = Nhat(1, DRule([lam1, lam2], H0, H1, H2), H1, alpha1);
    N_lam2 = Nhat(-1, DRule([lam1, lam2], H0, H1, H2), H2, alpha2);
    N_tilda1 = Nhat(1, DRule([1/LAM1_sort(r1+1), lam2], H0, H1, H2), H1, alpha1);
    N_tilda2 = Nhat(-1, DRule([lam1, 1/LAM2_sort(r2+1)], H0, H1, H2), H2, alpha2);
end

lam_star = [lam1, lam2];
end