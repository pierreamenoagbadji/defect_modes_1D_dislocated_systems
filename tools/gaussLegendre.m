function [ points, poids ] = gaussLegendre( Nquad, a, b )

% Genere une formule de quadrature de Gauss-Legendre a Nquad points pour
% calculer une integrale definie sur un intervalle [a,b]

% Matrice de Jacobi equivalente
alpha = zeros(1,Nquad);
j = 1:Nquad-1;
beta = j ./ sqrt( (2*j-1) .* (2*j+1) );

J = diag( alpha ) + diag( beta, 1 ) + diag( beta, -1 );

% Valeurs et Vecteurs propres
[V, D] = eig(J);

% Noeuds et Poids
points =  0.5*(b + a) + 0.5*(b - a) * diag(D).';
poids = (b - a) * abs(V(1,:)).^2;

end

