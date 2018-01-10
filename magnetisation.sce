fig_mag = scf(); clf(fig_mag);

N = 5;
n0 = 800;
T = 50;
J = ones(N,N,2)/T;
h = zeros(N,N)/T;

//Test pour gibbs séquentiel
n = n0;
n_burn = n/10;
n_simu = n + n_burn;

X1 = ising_gibbs_seq_chain(J,h,n_simu);
X1 = X1(:,:,(n_burn+1):(n_simu+1));
m = zeros(1,n+1);
for k = 1:n+1
    m(k) = mean(double(X1(:,:,k)));
end
subplot(2,2,1);
plot(1:n+1,cumsum(m) ./ (1:(n+1)),[1 n+1],[0 0]);
title("Ergodicité pour Gibbs séquentiel");

//Test pour gibbs randomisé
n = n0*N^2;
n_burn = n/10;
n_simu = n + n_burn;

X2 = ising_gibbs_rand_chain(J,h,n_simu);
X2 = X2(:,:,(n_burn+1):(n_simu+1));
m = zeros(1,n+1);
for k = 1:n+1
    m(k) = mean(double(X2(:,:,k)));
end
subplot(2,2,2);
plot(1:n+1,cumsum(m) ./ (1:(n+1)),[1 n+1],[0 0]);
title("Ergodicité pour Gibbs randomisé");

//Test pour gibbs randomisé
n = n0*N^2;
n_burn = n/10;
n_simu = n + n_burn;

X3 = ising_MH_chain(J,h,n_simu);
X3 = X3(:,:,(n_burn+1):(n_simu+1));
m = zeros(1,n+1);
for k = 1:n+1
    m(k) = mean(double(X3(:,:,k)));
end
subplot(2,2,3);
plot(1:n+1,cumsum(m) ./ (1:(n+1)),[1 n+1],[0 0]);
title("Ergodicité pour Metropolis-Hastings");

//Test pour coupling from the past
n = n0/10;
m = zeros(1,n+1);
for k = 1:n+1
    m(k) = mean(double(ising_coupling_MH(J,h,%f)));
end
subplot(2,2,4);
plot(1:n+1,cumsum(m) ./ (1:(n+1)),[1 n+1],[0 0]);
title("Ergodicité pour Coupling From The Past par Metropolis-Hastings");
