%function to build DF matrix
function DF = buildDF(w,a,b,mu)

%Daf1

n=length(a); 
Daf1 = zeros(1,n);
for j=0:n-1
    Daf1(j+1)=1;
end

%Dwf2, Dwf3
Dwf2 = zeros(n,1);
Dwf3 = zeros(n,1);
k=(n-1)/2;
order = (n-1)/2;
for j=-k:k
    aj=j*1i*accessFourierCoeff(j,a,order);
    Dwf2(fourierindex(j,order))=aj;
    bj=j*1i*accessFourierCoeff(j,b,order); 
    Dwf3(fourierindex(j,order))=bj;

end

%Daf2 = wML, matrix of L(x)_n = iwnx





Daf2=w*makeML(n);
L=Daf2;
Dbf2=-eye(n);

%Daf3 and M_aconv 

M_aconv = makeM_aa(a);
M_abconv = makeM_ab(a,b);
Daf3 = eye(n)+2*mu*M_abconv;

%Dbf3 = ML-muI +muM_aconv

Dbf3 = L-mu*eye(n)+mu*M_aconv;


%DF(w,a,b) = Dwf1, Daf1, Dbf1; Dwf2... 

Dcell = {0, Daf1, zeros(1,n); Dwf2, Daf2, Dbf2; Dwf3, Daf3, Dbf3};

DF = cell2mat(Dcell);




end