% This function returns N_hat in the teesting procedure.

function N = Nhat(k, delta, H, alpha)
N = sum((1-1*(delta==k)-alpha).*H);%/length(H);
end