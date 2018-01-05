//rand("seed",75); //seed fixée

/*Température critique : T = 1*/
N = 10;
J = ones(N,N,2);
h = zeros(N,N);

scf(1); clf(1);
printf("Temps d''exécution :\n");

tic();
T1 = 1.1;
X1 = ising_coupling_MH(J/T1,h/T1);
//X1 = ising_coupling_gibbs(J/T1,h/T1);
t1 = toc();
printf("T = "+string(T1)+" : "+string(t1)+"s\n");
subplot(2,2,1);
Matplot(X1+1);
title("T = "+string(T1));

tic();
T2 = 1.01;
X2 = ising_coupling_MH(J/T2,h/T2);
//X2 = ising_coupling_gibbs(J/T2,h/T2);
t2 = toc();
printf("T = "+string(T2)+" : "+string(t2)+"s\n");
subplot(2,2,2);
Matplot(X2+1);
title("T = "+string(T2));

tic();
T3 = 1;
X3 = ising_coupling_MH(J/T3,h/T3);
//X3 = ising_coupling_gibbs(J/T3,h/T3);
t3 = toc();
printf("T = "+string(T3)+" : "+string(t3)+"s\n");
subplot(2,2,3);
Matplot(X3+1);
title("T = "+string(T3));

tic();
T4 = 0.99;
X4 = ising_coupling_MH(J/T4,h/T4);
//X4 = ising_coupling_gibbs(J/T4,h/T4);
t4 = toc();
printf("T = "+string(T4)+" : "+string(t4)+"s\n");
subplot(2,2,4);
Matplot(X4+1);
title("T = "+string(T4));
