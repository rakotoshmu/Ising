/********************************
Comparaison entre les algorithmes
********************************/

N = 50; //taille du réseau
T = 3; //température sur-critique pour que les algorithmes de couplage terminent
J = ones(N,N,2)/T;
h = zeros(N,N)/T;

//nombres de balayages
n1 = 100;
n2 = n1 * N^2;
n3 = n2;

fig_comp = scf(); clf(fig_comp);

printf("Temps d''exécution :\n");

//Échantillonneur de Gibbs par balayage séquentiel
printf("Gibbs séquentiel :\n");
tic();
X1 = ising_gibbs_seq(J,h,n1);
t1 = toc();
printf("\t"+string(t1)+"s\n");
//Résultat
subplot(2,3,1);
Matplot(X1+1);
title("Gibbs, "+string(n1)+" balayages séquentiels");

//Échantillonneur de Gibbs par balayage randomisé
printf("Gibbs aléatoire :\n");
tic();
X2 = ising_gibbs_rand(J,h,n2);
t2 = toc();
printf("\t"+string(t2)+"s\n");
//Résultat
subplot(2,3,2);
Matplot(X2+1);
title("Gibbs, balayage aléatoire, "+string(n2)+" transitions");

//Metropolis-Hastings
printf("Metropolis-Hastings :\n");
tic();
X3 = ising_MH(J,h,n3);
t3 = toc();
printf("\t"+string(t3)+"s\n");
//Résultat
subplot(2,3,3);
Matplot(X3+1);
title("Metropolis-Hastings, "+string(n3)+" transitions");

//Coupling From The Past par Metropolis-Hastings
printf("Coupling Metropolis-Hastings :\n");
tic();
X4 = ising_coupling_MH(J,h,%f);
t4 = toc();
printf("\t"+string(t4)+"s\n");
//Résultat
subplot(2,3,4);
Matplot(X4+1);
title("Coupling Metropolis-Hastings");

//Coupling From The Past par Gibbs
printf("Coupling Gibbs :\n");
tic();
X5 = ising_coupling_gibbs(J,h,%f); //lancer l'algo sans feedback
t5 = toc();
printf("\t"+string(t5)+"s\n");
//Résultat
subplot(2,3,5);
Matplot(X5+1);
title("Coupling Gibbs");

//Donner la température
subplot(2,3,6);
title("T = "+string(T));
