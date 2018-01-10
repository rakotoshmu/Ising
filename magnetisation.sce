//La magnétisation est d'espérance nulle
N = 5;
J = ones(N,N,2);
h = zeros(N,N);
n = 1000;
m = 1;
magnet = 0;

for k = 1:50
    magnet = magnet + mean(double(ising_gibbs_rand(J,h,n)))/m;
end
