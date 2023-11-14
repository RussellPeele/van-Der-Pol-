function Moperator = makeM_aa(a)
%for j=-2:2
 %   a(3+j)=j
%end

N=length(a);
ordera = (N-1)/2;
aa = conv(a,a,'same');

M =length(aa);
orderaa=(M-1)/2;
M_aconv = zeros(M); 
for j=0:M-1
   for k=0:M-1
       
            M_aconv(j+1, k+1)=accessFourierCoeff(j-k,aa,orderaa);
      
    end
end


Moperator=M_aconv;
end