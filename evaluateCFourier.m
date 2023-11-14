function fx = evaluateCFourier(A, x, w, N)

sum = 0;
for n = -N:N
    sum = sum + A(fourierindex(n, N))*exp(w*1i*n*x); 
end
fx = sum;