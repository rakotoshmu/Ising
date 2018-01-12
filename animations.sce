/*************************
Animations des algorithmes
*************************/

N = 20; //taille du réseau
T = 2; //température sur-critique pour que les algorithmes de couplage terminent
J = ones(N,N,2)/T;
h = zeros(N,N)/T;
n1 = 10;
n2 = n1 * N^2;
n3 = n2;

printf("Temps d''exécution :\n");

//Échantillonneur de Gibbs par balayage séquentiel
tic();
X1 = ising_gibbs_seq_chain(J,h,n1);
t1 = toc();
printf("Gibbs séquentiel : "+string(t1)+"s\n");
//Animation
scf(1);
for k = 1:size(X1,3)
    drawlater;
    clf(1);
    Matplot(X1(:,:,k)+1);
    title("Gibbs séquentiel, k = "+string(k-1));
    drawnow;
    sleep(100);
end

//Échantillonneur de Gibbs par balayage randomisé
tic();
X2 = ising_gibbs_rand_chain(J,h,n2);
t2 = toc();
printf("Gibbs aléatoire : "+string(t2)+"s\n");
//Animation
scf(2);
for k = 1:10:size(X2,3)
    drawlater;
    clf(2);
    Matplot(X2(:,:,k)+1);
    title("Gibbs aléatoire, k = "+string(k-1));
    drawnow;
end

//Metropolis-Hastings
tic();
X3 = ising_MH_chain(J,h,n3);
t3 = toc();
printf("Metropolis-Hastings : "+string(t3)+"s\n");
//Animation
scf(3);
for k = 1:10:size(X3,3)
    drawlater;
    clf(3);
    Matplot(X3(:,:,k)+1);
    title("Metropolis-Hastings, k = "+string(k-1));
    drawnow;
end

//Coupling From The Past par Metropolis-Hastings
tic();
X4 = ising_coupling_MH(J,h,%t); //afficher l'animation "attente de coalition"
t4 = toc();
printf("Coupling Metropolis-Hastings : "+string(t4)+"s\n");
//Résultat
scf(4);
clf(4);
Matplot(X4+1);
title("Coupling Metropolis-Hastings");

//Coupling From The Past par Gibbs
tic();
X5 = ising_coupling_gibbs(J,h,%t); //afficher l'animation "attente de coalition"
t5 = toc();
printf("Coupling Gibbs : "+string(t5)+"s\n");
//Résultat
scf(5);
clf(5);
Matplot(X5+1);
title("Coupling Gibbs");
