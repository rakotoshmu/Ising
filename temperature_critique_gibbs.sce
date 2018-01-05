//rand("seed",75); //seed fixée

/*
Température critique : T = 1
Sur mon ordinateur, pour N = 100, n fixé, ising_gibbs_seq met un peu plus de n minutes à terminer
*/
N = 100;
J0 = ones(N,N,2);
h0 = zeros(N,N);
n = 5;

scf(1); clf(1);
printf("Temps d''exécution :\n");

tic();
T1 = 1.1;
//X1 = ising_gibbs_rand(J0/T1,h0/T1,n);
X1 = ising_gibbs_seq(J0/T1,h0/T1,n);
t1 = toc();
printf("T = "+string(T1)+" : "+string(t1)+"s\n");
subplot(2,2,1);
Matplot(X1+1);
title("T = "+string(T1)+" (n = "+string(n)+")");

tic();
T2 = 1.01;
//X2 = ising_gibbs_rand(J0/T2,h0/T2,n);
X2 = ising_gibbs_seq(J0/T2,h0/T2,n);
t2 = toc();
printf("T = "+string(T2)+" : "+string(t2)+"s\n");
subplot(2,2,2);
Matplot(X2+1);
title("T = "+string(T2)+" (n = "+string(n)+")");

tic();
T3 = 1;
//X3 = ising_gibbs_rand(J0/T3,h0/T3,n);
X3 = ising_gibbs_seq(J0/T3,h0/T3,n);
t3 = toc();
printf("T = "+string(T3)+" : "+string(t3)+"s\n");
subplot(2,2,3);
Matplot(X3+1);
title("T = "+string(T3)+" (n = "+string(n)+")");

tic();
T4 = 0.99;
//X4 = ising_gibbs_rand(J0/T4,h0/T4,n);
X4 = ising_gibbs_seq(J0/T4,h0/T4,n);
t4 = toc();
printf("T = "+string(T4)+" : "+string(t4)+"s\n");
subplot(2,2,4);
Matplot(X4+1);
title("T = "+string(T4)+" (n = "+string(n)+")");
