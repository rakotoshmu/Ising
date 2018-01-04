//rand("seed",75); //seed fixée

N = 20;
J = ones(N,N,2);
h = zeros(N,N);
n1 = 100;
n2 = 2000;
n3 = 2000;

printf("Temps d''exécution :\n");
tic();
X1 = ising_gibbs_seq(J,h,n1);
t1 = toc();
printf("Gibbs séquentiel : "+string(t1)+"\n");
tic();
X2 = ising_gibbs_rand(J,h,n2);
t2 = toc();
printf("Gibbs aléatoire : "+string(t2)+"\n");
tic();
X3 = ising_MH(J,h,n3,0.25);
t3 = toc();
printf("Métropolis-Hastings : "+string(t3)+"\n");





/*********
Animations
*********/

//scf(1);
//for k = 1:size(X1,3)
//    drawlater;
//    clf(1);
//    Matplot(X1(:,:,k)+1);
//    title("Gibbs séquentiel, k = "+string(k-1));
//    drawnow;
//    sleep(100);
//end
//scf(2);
//for k = 1:size(X2,3)
//    drawlater;
//    clf(2);
//    Matplot(X2(:,:,k)+1);
//    title("Gibbs aléatoire, k = "+string(k-1));
//    drawnow;
//end
//f = scf(3);
//f.color_map = cmap;
//for k = 1:size(X3,3)
//    drawlater;
//    clf(3);
//    Matplot(X3(:,:,k)+1);
//    title("Metropolis-Hastings, k = "+string(k-1));
//    drawnow;
//end
scf(4); clf(4);
subplot(2,2,1);
k = 1;
Matplot(X2(:,:,k)+1);
title("Gibbs séquentiel, k = "+string(k-1));
subplot(2,2,2);
k = 501;
Matplot(X2(:,:,k)+1);
title("Gibbs séquentiel, k = "+string(k-1));
subplot(2,2,3);
k = 1001;
Matplot(X2(:,:,k)+1);
title("Gibbs séquentiel, k = "+string(k-1));
subplot(2,2,4);
k = 2001;
Matplot(X2(:,:,k)+1);
title("Gibbs séquentiel, k = "+string(k-1));
