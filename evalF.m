function Fwab = evalF(w,a, b,azero,mu,N)
Daf1 = zeros(1,N);
for j=0:N-1
    Daf1(j+1)=1;
end



    Fwab=[0; zeros(N,1);zeros(N,1)];
    L=w*makeML(N);
    Fwab1 = sum(a)-azero;
    Fwab2=L*a'-b';
    trunc = mu*(conv(conv(a,a,'same'),b,'same'));
    Fwab3= L*b'-mu*b' +a' + trunc';
   
    
   
  

Fwab = [Fwab1,Fwab2',Fwab3'];

end