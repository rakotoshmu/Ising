//rand("seed",75); //seed fixée

/*Température critique : T_c = 1
Pour T<1, les algorithmes de couplage ne terminent pas (en temps court)*/
N = 10;
J = ones(N,N,2);
h = zeros(N,N);
T1 = 1.1;
T2 = 1.01;
T3 = 1;
T4 = 0.99;

/*Température critique ?
Pour T<T_c, les algorithmes de couplage ne terminent pas (en temps court)
N = 10;
J = ones(N,N,2);
h = 5*ones(N,N);
T1 = 10;
T2 = 7.5;
T3 = 5;
T4 = 2.5;*/


fig_tcc = scf(); clf(fig_tcc);
printf("Temps d''exécution :\n");

tic();
X1 = ising_coupling_MH(J/T1,h/T1,%t);
//X1 = ising_coupling_gibbs(J/T1,h/T1,%t);
t1 = toc();
printf("T = "+string(T1)+" : "+string(t1)+"s\n");
subplot(2,2,1);
Matplot(X1+1);
title("T = "+string(T1));

tic();
X2 = ising_coupling_MH(J/T2,h/T2,%t);
//X2 = ising_coupling_gibbs(J/T2,h/T2,%t);
t2 = toc();
printf("T = "+string(T2)+" : "+string(t2)+"s\n");
subplot(2,2,2);
Matplot(X2+1);
title("T = "+string(T2));

tic();
X3 = ising_coupling_MH(J/T3,h/T3,%t);
//X3 = ising_coupling_gibbs(J/T3,h/T3,%t);
t3 = toc();
printf("T = "+string(T3)+" : "+string(t3)+"s\n");
subplot(2,2,3);
Matplot(X3+1);
title("T = "+string(T3));

tic();
X4 = ising_coupling_MH(J/T4,h/T4,%t);
//X4 = ising_coupling_gibbs(J/T4,h/T4,%t);
t4 = toc();
printf("T = "+string(T4)+" : "+string(t4)+"s\n");
subplot(2,2,4);
Matplot(X4+1);
title("T = "+string(T4));
