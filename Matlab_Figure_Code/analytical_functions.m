%exact solution
pA = @(a,b,k,l) max(0, min(1, (a./k - b./l + a.*(k - 2)./k) ./ (b.*(l-2)./l + a.*(k-2)./k)));
rBA = @(a,b,k,l) (k-2)./(k-1) .* (1-pA(a,b,k,l));
qAB = @(a,b,k,l) (l-2)./(l-1) .* pA(a,b,k,l);

pAWM = @(a,b) a ./ (a + b);
Wwm = @(a, b) pAWM(a, b) .* a .* (1 - pAWM(a, b)) + (1 - pAWM(a, b)) .* b .* pAWM(a, b);
Wpa = @(a,b,k,l) pA(a,b,l,k) .* a .* rBA(a,b,k,l) + (1 - pA(a,b,l,k)) .* b .* qAB(a,b,k,l) ;
Wrel = @(a,b,k,l) Wpa(a,b,k,l) ./ Wwm(a,b);

muRelExtict = @(k, l) l .* (k-1) ./ k; %use if muB > muA
%muRel2 = @(k, l) l ./ (k .* (l-1)); %use if muA > muB

