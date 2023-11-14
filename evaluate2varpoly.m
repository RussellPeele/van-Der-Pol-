function P_sigma = evaluate2varpoly(P, sigma1, sigma2, N)

P_sigma = [0;0;0];
for n = 0:N
    for m = 0:N
    P_sigma = P_sigma + P(:, indexp(n,m,N))*sigma1^n*sigma2^m;
    end
end