function y = focksum(cs,x)
%FOCKSUM linear combination of harmonic oscillator wave functions
%
%   FOCKSUM(CS,X) where CS and X are vectors, sums a series with
%   coefficients CS and evaluates the wave functions at the points X.
%   
%   FOCKSUM(CS,X) where X is a vector and CS is a matrix, sums the
%   series with coefficients CS(:,i) at the point X(i).
%   
%   The terms of the series are CS(n+1)*F_n(x), where F is the harmonic
%   oscillator wave function evaluated by FOCKSTATE.
%   
%   F_n is a wave function for the harmonic oscillator position quadrature
%   (a+a')/sqrt(2), where a is the lowering operator.  A wave function
%   for the general quadrature
%   
%     Y = (a*e^(-it) + a'*e^(it))/sqrt(2)
%   
%   can be evaluated by propagating CS through t/2*pi of a cycle, i.e.,
%   Y = FOCKSUM(CS.*exp(-i*(0:n)*t), XS).  For the momentum quadrature
%   Y = i*(a'-a)/sqrt(2), use FOCKSUM(CS.*(-i.^(0:n)), XS).
% 
%   See also: FOCKSTATE

% For the algorithm, see Clenshaw, "A note on the summation of
% Chebyshev series", Mathematics of Computation 9, p118 (1955).

assert(isvector(x))
if isvector(cs), cs = cs(:); end

N = size(cs,1)-1;
xs = x(:);

% this evaluates a polynomial F_n = H_n/sqrt(2^n*n!), whose
% recurrence relation can be easily derived from that of H_n.

B = zeros(length(xs), 2);	% [ b_{N+1} b_{N+2} ]
logB = zeros(length(xs), 1);	% scale by 2^-logB to prevent overflow
for n = N:-1:0
	a = -xs*sqrt(2/(n+1));	% a_n from F_n recurrence relation
	b = sqrt((n+1)/(n+2));	% and b_{n+1}
	B = [cs(n+1,:).' - a.*B(:,1) - b.*B(:,2),  B(:,1)];
	[F, E] = log2(B(:,1));
	logB = logB + E;
	B = [F B(:,2)./2.^E];
end

y = pi^(-1/4)*exp(-xs.^2/2+log(2)*logB).*B(:,1);	% reverse the scaling
y = reshape(y, size(x));

end