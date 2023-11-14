function Moperator = makeM_ab(a,b)
%for j=-2:2
 %   a(3+j)=j
%end

N=length(a);
ordera = (N-1)/2;
ab = conv(a,b,'same');

M =length(ab);
orderab=(M-1)/2;
M_aconv = zeros(M); 
for j=0:M-1
   for k=0:M-1
       
            M_aconv(j+1, k+1)=accessFourierCoeff(j-k,ab,orderab);
      
    end
end

Moperator=M_aconv;
end