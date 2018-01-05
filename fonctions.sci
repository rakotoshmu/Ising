//clear;

/*
Dans tout le fichier, les paramètres du modèle d'Ising seront :
- N entier le côté du réseau carré de dimension 2 sur lequel on simule le modèle d'Ising
- J de taille N x N x 2 matrice des forces d'intéractions
- h de taille N x 1 matrice du champ magnétique extérieur
*/





/******************
Metropolis-Hastings
******************/

function V = V_u(J,h,x,i,j)
    /*
    Calcule et renvoie V_u(x) où u = (i,j) et x est un état (de taille N x N)
    N,J,h les paramètres du modèle d'Ising
    */
    N = size(h,1);
    V = h(i,j);
    if i>1 then
        V = V + J(i-1,j,1) * x(i-1,j);
    end
    if i<N then
        V = V + J(i,j,1) * x(i+1,j);
    end
    if j>1 then
        V = V + J(i,j-1,2) * x(i,j-1);
    end
    if j<N then
        V = V + J(i,j,2) * x(i,j+1);
    end
endfunction

function X = ising_MH_chain(J,h,n)
    /*
    Simule le modèle d'Ising selon un algorithme de Metropolis-Hasting, en partant d'un état initial aléatoire uniforme
    N,J,h les paramètres du modèle d'Ising
    n entier le temps de simulation de la chaine
    p dans [0,1] probabilité de changer chaque coordonnée pour la proposition de M-H
    Renvoie X de taille N x N x n+1 (X(:,:,k) état de la chaîne de Markov à l'instant k+1)
    */
    N = size(h,1);

    X = int8(ones(N,N,n+1));
    X(:,:,1) = 2*int8((grand(N,N,"def")<0.5))-1;
    I = ceil(N*grand(2,n,"def")); //indices aléatoires
    S = 2*int8(grand(1,n,"def")<1/2)-1; //spin aléatoires
    U = grand(1,n,"def"); //uniformes sur [0,1]
    for k = 1:n
        X(:,:,k+1) = X(:,:,k);
        if X(I(1,k),I(2,k),k)~=S(k) then
            v = V_u(J,h,X(:,:,k),I(1,k),I(2,k));
            if U(k) < exp(double(2*S(k)*v)) then
                X(I(1,k),I(2,k),k+1) = S(k);
            end
        end
    end
endfunction

function X = ising_MH(J,h,n)
    /*
    Comme ising_MH_chain, mais ne renvoie que l'état final
    */
    N = size(h,1);

    X = 2*int8((grand(N,N,"def")<0.5))-1;
    I = ceil(N*grand(2,n,"def")); //indices aléatoires
    S = 2*int8(grand(1,n,"def")<1/2)-1; //spin aléatoires
    U = grand(1,n,"def"); //uniformes sur [0,1]
    for k = 1:n
        if X(I(1,k),I(2,k))~=S(k) then
            v = V_u(J,h,X,I(1,k),I(2,k));
            if U(k) < exp(double(2*S(k)*v)) then
                X(I(1,k),I(2,k)) = S(k);
            end
        end
    end
endfunction





/**********************
Échantilloneur de Gibbs
**********************/

function Y = ising_gibbs_step(J,h,X,i,j,u)
    /*
    Effectue un pas de l'échantilloneur de Gibbs selon la coordonnée (i,j)
    Renvoie Y le résultat de ce pas
    N,J,h les paramètres du modèle d'Ising
    X de taille N x N l'état de départ
    i,j coordonnées
    u dans [0,1] tiré uniformément
    */
    N = size(h,1);

    p = 1/(1+exp(double(-2*V_u(J,h,X,i,j))));

    Y = X;
    Y(i,j) = 2*int8(u<p)-1;
endfunction

function X = ising_gibbs_seq_chain(J,h,n)
    /*
    Simule le modèle d'Ising par l'échantilloneur de Gibbs avec balayage séquentiel, en partant d'un état initial aléatoire uniforme
    N,J,h les paramètres du modèle d'Ising
    n entier le temps de simulation de la chaine
    Renvoie X de taille N x N x n+1 (X(:,:,k) état de la chaîne de Markov à l'instant k+1)
    */
    N = size(h,1);

    X = int8(ones(N,N,n+1));
    X(:,:,1) = 2*int8((grand(N,N,"def")<0.5))-1;
    for k = 1:n
        Xtemp = X(:,:,k);
        U = grand(N,N,"def");
        //balayage séquentiel
        for i = 1:N
            for j = 1:N
                Xtemp = ising_gibbs_step(J,h,Xtemp,i,j,U(i,j));
            end
        end
        X(:,:,k+1) = Xtemp;
    end
endfunction

function X = ising_gibbs_seq(J,h,n)
    /*
    Comme ising_gibbs_seq_chain, mais ne renvoie que l'état final
    */
    N = size(h,1);

    X = 2*int8((grand(N,N,"def")<0.5))-1;
    for k = 1:n
        U = grand(N,N,"def");
        //balayage séquentiel
        for i = 1:N
            for j = 1:N
                X = ising_gibbs_step(J,h,X,i,j,U(i,j));
            end
        end
    end
endfunction

function X = ising_gibbs_rand_chain(J,h,n)
    /*
    Simule le modèle d'Ising par l'échantilloneur de Gibbs avec balayage séquentiel, en partant d'un état initial aléatoire uniforme
    N,J,h les paramètres du modèle d'Ising
    n entier le temps de simulation de la chaine
    Renvoie X de taille N x N x n+1 (X(:,:,k) état de la chaîne de Markov à l'instant k+1)
    */
    N = size(h,1);

    X = int8(ones(N,N,n+1));
    X(:,:,1) = 2*int8((grand(N,N,"def")<0.5))-1;
    I = ceil(N*grand(2,n,"def")); //indices aléatoires
    U = grand(1,n,"def"); //uniformes
    for k = 1:n //balayage aléatoire
        X(:,:,k+1) = ising_gibbs_step(J,h,X(:,:,k),I(1,k),I(2,k),U(k));
    end
endfunction

function X = ising_gibbs_rand(J,h,n)
    /*
    Comme ising_gibbs_rand_chain, mais ne renvoie que l'état final
    */
    N = size(h,1);

    X = 2*int8((grand(N,N,"def")<0.5))-1;
    I = ceil(N*grand(2,n,"def")); //indices aléatoires
    U = grand(1,n,"def"); //uniformes
    for k = 1:n //balayage aléatoire
        X = ising_gibbs_step(J,h,X,I(1,k),I(2,k),U(k));
    end
endfunction





/********************
Couplage par le passé
********************/

function X = ising_coupling_MH(J,h)
    /*
    Simule le modèle d'Ising par couplage par le passé sur Metropolis-Hastings
    N,J,h les paramètres du modèle d'Ising
    Renvoie X de taille N x N
    */
    N = size(h,1);

    X = int8(ones(N,N));
    Y = -int8(ones(N,N));
    counter = 0;
    fig = scf();

    while max(abs(X-Y))>0 do
        i = ceil(N*grand(1,1,"def")); j = ceil(N*grand(1,1,"def")); //indices aléatoires
        s = 2*int8(grand(1,1,"def")<1/2)-1; //spin aléatoire
        u = grand(1,1,"def"); //uniforme sur [0,1]

        if X(i,j)~=s then
            v = V_u(J,h,X,i,j);
            if u < exp(double(2*s*v)) then
                X(i,j) = s;
            end
        end

        if Y(i,j)~=s then
            v = V_u(J,h,Y,i,j);
            if u < exp(double(2*s*v)) then
                Y(i,j) = s;
            end
        end

        counter = counter + 1;
        /* Permet d'afficher l'avancement du couplage
        */
        drawlater;
        clf(fig);
        subplot(2,2,1);
        title("Nombre d''itérations : "+string(counter));
        subplot(2,2,2);
        title("Nombre de spins différents : "+string(sum(double(abs(X-Y))/2)));
        subplot(2,2,3);
        Matplot(X+1); title("1");
        subplot(2,2,4);
        Matplot(Y+1); title("-1");
        drawnow;
    end

    close(fig);
    printf("\t"+string(counter)+" itérations pour coupling from the past via MH\n");
endfunction

function X = ising_coupling_gibbs(J,h)
    /*
    Simule le modèle d'Ising par couplage par le passé sur l'échantillonneur de Gibbs
    N,J,h les paramètres du modèle d'Ising
    Renvoie X de taille N x N
    */
    N = size(h,1);

    X = int8(ones(N,N));
    Y = -int8(ones(N,N));
    counter = 0;
    fig = scf();

    while max(abs(X-Y))>0 do
        U = grand(N,N,"def");
        //balayage séquentiel
        for i = 1:N
            for j = 1:N
                X = ising_gibbs_step(J,h,X,i,j,U(i,j));
                Y = ising_gibbs_step(J,h,Y,i,j,U(i,j));
            end
        end
        counter = counter + 1;
        /* Permet d'afficher l'avancement du couplage
        */
        drawlater;
        clf(fig);
        subplot(2,2,1);
        title("Nombre d''itérations : "+string(counter));
        subplot(2,2,2);
        title("Nombre de spins différents : "+string(sum(double(abs(X-Y))/2)));
        subplot(2,2,3);
        Matplot(X+1); title("1");
        subplot(2,2,4);
        Matplot(Y+1); title("-1");
        drawnow;
        
    end

    close(fig);
    printf("\t"+string(counter)+" itérations pour coupling from the past via Gibbs\n");
endfunction





/**********************
Simulation exacte naïve
**********************/

function X = num2etat(N,m)
    /*
    Calcule le m-ième état dans l'énumération lexicographique des états
    N entier la taille du côté du réseau
    m entier entre 1 et 2^(N^2) le numéro de l'état
    Renvoie X le m-ième état
    */
    //a = matrix(1:N^2,N,N) est une matrice de taille N x N les entiers de 1 à N^2
    //b = bitget(m,a) est une matrice de taille NxN listant les bit de m
    //c = 2*b-1 transforme ces {0,1} et {-1,1}
    X = 2*int8(bitget(m,matrix(1:N^2,N,N)))-1;
endfunction

function m = etat2num(X)
    /*
    Calcule le numéro de l'état dans l'énumération lexicographique des états
    N entier la taille du côté du réseau
    X de taille N x N à valeur {1,-1}
    Renvoie m le numéro de l'état
    */
    //a = matrix(X,1,-1) est le vecteur ligne composé des éléments de X
    //b = (a+1)/2 transforme les {-1,1} et {0,1} donc b(k) est le k-ième bit de m
    m = sum(2.^(0:N^2-1) .* (matrix(X,1,-1)+1)/2);
endfunction

function y = pi(J,h,x)
    /*
    Calcule la loi non normalisée
    N,J,h les paramètres du modèle d'Ising
    x de taille N x N
    Renvoie y = Z_T * pi(x)
    */
    N = size(h,1);

    //intercation entre les voisins de même ordonnée
    s1 = sum(J(1:N-1,:,1).*x(1:N-1,:).*x(2:N,:)); 
    //intercation entre les voisins de même abscisse
    s2 = sum(J(:,1:N-1,2).*x(:,1:N-1).*x(:,2:N)); 
    //champ magnétique extérieur
    s3 = sum(h.*x);

    y = exp(double(s1+s2+s3));
endfunction

function p = ising_law(J,h)
    /*
    Renvoie la loi de probabilité du modèle d'Ising
    p(m)/sum(p) est la probabilité du m-ième état
    N,J,h les paramètres du modèle d'Ising
    */
    N = size(h,1);
    d = 2^(N^2); //nombre d'états

    p = zeros(1,d);
    for m = 1:d
        x = num2etat(N,m);
        p(m) = pi(J,h,x); //probabilité non normalisée
    end
endfunction

function X = ising_exact(J,h,n)
    /*
    Simule un n-échantillon du modèle d'Ising de manière naïve (probabilités discrètes), à ne lancer qu'avec N petit
    N,J,h les paramètres du modèle d'Ising
    n le nombre de tirage
    Renvoie X de taille N x N x n
    */
    N = size(h,1);

    p = ising_law(J,h);
    c = cumsum(p);
    for k = 1:n
        //m suit la loi p
        [t,m] = max(1*(grand(1,1,"def")*c(n)<c));
        X(:,:,k) = num2etat(N,m);
    end
endfunction
