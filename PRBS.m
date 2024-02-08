function z=PRBS(init,g)
% z=prbs(init,g)
% 2^n-1-bit PRBS based en condiciones iniciales ejm: [1 0 1 1 0 1 1] para
% n=7
% representa al polinomio generador g --> x^7+x^6+1 ejm:  g=[7 6]

z=init;
n=length(init);
for i=(n+1):(2^n-1)
q=z(i-g(1));
    for j=2:length(g)
        q=xor(q,z(i-g(j)));
    end
z=[z q];
end