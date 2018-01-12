/*Température critique :
T_c théorique ~ 2.269*/
N = 4;

T = linspace(0.01,5);
n = length(T);
tic();
U = ising_energy_std(N);
t = toc();
printf("Énergie calculée en "+string(t)+"s\n");
p = zeros(T);

for k = 1:n
    temp = exp(-U/T(k));
    p(k) = temp(length(U))/sum(temp);
end

scf(1); clf(1);
plot(T,p);
title("Probabilité des états (1) et (-1) en fonction de la température");
