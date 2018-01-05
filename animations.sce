//rand("seed",75); //seed fixée

N = 10;
T = 1.1;
J = ones(N,N,2)/T;
h = zeros(N,N)/T;
n1 = 100;
n2 = 1000;
n3 = 1000;

printf("Temps d''exécution :\n");

tic();
X1 = ising_gibbs_seq_chain(J,h,n1);
t1 = toc();
printf("Gibbs séquentiel : "+string(t1)+"s\n");

tic();
X2 = ising_gibbs_rand_chain(J,h,n2);
t2 = toc();
printf("Gibbs aléatoire : "+string(t2)+"s\n");

tic();
X3 = ising_MH_chain(J,h,n3);
t3 = toc();
printf("Metropolis-Hastings : "+string(t3)+"s\n");

tic();
X4 = ising_coupling_MH(J,h);
t4 = toc();
printf("Coupling Metropolis-Hastings : "+string(t4)+"s\n");

tic();
X5 = ising_coupling_gibbs(J,h);
t5 = toc();
printf("Coupling Metropolis-Hastings : "+string(t5)+"s\n");





/*********
Animations
*********/

scf(1);
for k = 1:size(X1,3)
    drawlater;
    clf(1);
    Matplot(X1(:,:,k)+1);
    title("Gibbs séquentiel, k = "+string(k-1));
    drawnow;
    sleep(100);
end

scf(2);
for k = 1:size(X2,3)
    drawlater;
    clf(2);
    Matplot(X2(:,:,k)+1);
    title("Gibbs aléatoire, k = "+string(k-1));
    drawnow;
end

scf(3);
for k = 1:size(X3,3)
    drawlater;
    clf(3);
    Matplot(X3(:,:,k)+1);
    title("Metropolis-Hastings, k = "+string(k-1));
    drawnow;
end

scf(4);
clf(4);
Matplot(X4+1);
title("Coupling Metropolis-Hastings");

scf(5);
clf(5);
Matplot(X5+1);
title("Coupling Gibbs");
