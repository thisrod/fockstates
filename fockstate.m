function y = fockstate(n,x,flag)
%FOCKSTATE harmonic oscillator wave function
%
%   FOCKSTATE(N,X) evaluates the harmonic oscillator eigenstate F_n(x)
%   at the elements of X.
%   
%   FOCKSTATE(N,X,'matrix'), where X is a vector of length m, returns
%   an m-by-(N+1) matrix whose ith column is F_{i-1}(x).
%   
%   The eigenstates are
%   
%    F_n(x) = 1/sqrt(2^n*factorial(n))*sqrt(pi))*exp(-x^2/2)*H(n,x),
%   
%   where H is a Hermite polynomial.  These satisfy the Schrodinger
%   equation
%   
%     -(1/2) F''_n(x) + x^2/2 F_n(x) = (n+1/2) F_n(x).
%   
%   In terms of the annihilation operator a, the quadrature x is 
%   x = (a+a')/sqrt(2).
%   
%   The Hermite polynomial H_n can be evaluated by the function
%   
%     H = @(n,x) sqrt(2^n*factorial(n)*sqrt(pi))*exp(-x^2/2) ...
%       .*FOCKSTATE(n,x)
%   
%   The function FOCKSTATE is stable well past the point where the
%   Hermite polynomials overflow the double precision floating point
%   numbers.
%   
%   See also: FOCKSUM

if nargin == 3 && strcmp(flag, 'matrix')
	y = focksum(kron(speye(n+1), ones(1,length(x))), repmat(x(:), n+1, 1));
	y = reshape(y, length(x), n+1);
else
	xs = x(:);
	y = focksum([zeros(1,n) 1], xs);
	y = reshape(y, size(x));
end

end
