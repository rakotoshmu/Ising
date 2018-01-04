clear;

/*
Dans tout le fichier, les paramètres du modèle d'Ising seront :
- N entier le côté du réseau carré de dimension 2 sur lequel on simule le modèle d'Ising
- J de taille N x N x 2 matrice des forces d'intéractions
- h de taille N x 1 matrice du champ magnétique extérieur
*/





/******************
Metropolis-Hastings
******************/

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

function X = ising_MH(J,h,n,p)
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
    M = 2*int8((grand(N,N,n,"def")<1-p))-1; //matrices par lesquelles multiplier X pour avoir la transition
    for k = 1:n
        Y = X(:,:,k) .* M(:,:,k);
        b = int8(grand(1,1,"def") < min(1,pi(J,h,Y)/pi(J,h,X(:,:,k))));
        X(:,:,k+1) = (1-b)*X(:,:,k) + b*Y;
    end
endfunction





/**********************
Échantilloneur de Gibbs
**********************/

function Y = ising_gibbs_step(J,h,X,i,j)
    /*
    Effectue un pas de l'échantilloneur de Gibbs selon la coordonnée (i,j)
    Renvoie Y le résultat de ce pas
    N,J,h les paramètres du modèle d'Ising
    X de taille N x N l'état de départ
    i,j coordonnées
    */
    N = size(h,1);

    //l = \lambda_u où u = (i,j)
    l = h(i,j);
    if i>1 then
        l = l + J(i-1,j,1) * X(i-1,j);
    end
    if i<N then
        l = l + J(i,j,1) * X(i+1,j);
    end
    if j>1 then
        l = l + J(i,j-1,2) * X(i,j-1);
    end
    if j<N then
        l = l + J(i,j,2) * X(i,j+1);
    end
    p = 1/(1+exp(-2*double(l)));

    Y = X;
    Y(i,j) = 2*int8((grand(1,1,"def")<p))-1;
endfunction

function X = ising_gibbs_seq(J,h,n)
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
        //balayage séquentiel
        for i = 1:N
            for j = 1:N
                Xtemp = ising_gibbs_step(J,h,Xtemp,i,j);
            end
        end
        X(:,:,k+1) = Xtemp;
    end
endfunction

function X = ising_gibbs_rand(J,h,n)
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
    for k = 1:n //balayage aléatoire
        X(:,:,k+1) = ising_gibbs_step(J,h,X(:,:,k),I(1,k),I(2,k));
    end
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

function X = ising_exact(J,h,n)
    /*
    Simule un n-échantillon du modèle d'Ising de manière naïve (probabilités discrètes), à ne lancer qu'avec N petit
    N,J,h les paramètres du modèle d'Ising
    n le nombre de tirage
    Renvoie X de taille N x N x n
    */
    N = size(h,1);
    d = 2^(N^2); //nombre d'états
    
    p = zeros(1,d);
    for m = 1:d
        x = num2etat(m);
        p(m) = pi(J,h,x); //probabilité non normalisée
    end
    c = cumsum(p);
    for k = 1:n
        //m suit la loi p
        [t,m] = max(1*(grand(1,1,"def")*c(n)<c));
        X(:,:,k) = num2etat(m);
    end
endfunction
