function fx = evaluateFourier(x, a0, A, B, N)

sum = a0;
for n = 1:N
    sum = sum + A(n)*sin(n*x) + B(n)*cos(n*x); 
end
fx = sum;