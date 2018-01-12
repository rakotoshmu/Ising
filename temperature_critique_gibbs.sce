/*
Température critique :
T_c théorique ~ 2.269
*/
N = 100;
J = ones(N,N,2);
h = zeros(N,N);
n = 50;
//T1 = 2.28; T2 = 2.27; T3 = 2.26; T4 = 2.25;
T1 = 2.4; T2 = 2.3; T3 = 2.2; T4 = 2.1;
//T1 = 3; T2 = 2.5; T3 = 2; T4 = 1.5;

fig_tcg = scf(); clf(fig_tcg);
printf("Temps d''exécution :\n");

tic();
//X1 = ising_gibbs_rand(J/T1,h/T1,n);
X1 = ising_gibbs_seq(J/T1,h/T1,n);
t1 = toc();
printf("T = "+string(T1)+" : "+string(t1)+"s\n");
subplot(2,2,1);
Matplot(X1+1);
title("T = "+string(T1)+" (n = "+string(n)+")");

tic();
//X2 = ising_gibbs_rand(J/T2,h/T2,n);
X2 = ising_gibbs_seq(J/T2,h/T2,n);
t2 = toc();
printf("T = "+string(T2)+" : "+string(t2)+"s\n");
subplot(2,2,2);
Matplot(X2+1);
title("T = "+string(T2)+" (n = "+string(n)+")");

tic();
//X3 = ising_gibbs_rand(J/T3,h/T3,n);
X3 = ising_gibbs_seq(J/T3,h/T3,n);
t3 = toc();
printf("T = "+string(T3)+" : "+string(t3)+"s\n");
subplot(2,2,3);
Matplot(X3+1);
title("T = "+string(T3)+" (n = "+string(n)+")");

tic();
//X4 = ising_gibbs_rand(J/T4,h/T4,n);
X4 = ising_gibbs_seq(J/T4,h/T4,n);
t4 = toc();
printf("T = "+string(T4)+" : "+string(t4)+"s\n");
subplot(2,2,4);
Matplot(X4+1);
title("T = "+string(T4)+" (n = "+string(n)+")");
