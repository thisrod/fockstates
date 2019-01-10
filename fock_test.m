% tests of fockstate and focksum

xx = -4:0.1:4;
figure, hold on
plot(xx,xx.^2/2,'--k')
for n = 0:5
	cs = [zeros(1,n) 1];
	plot(xx, 0.5*focksum(cs, xx)+1/2+n, '-k')
end

figure, hold on
xx = -30:0.01:30;
plot(xx,xx.^2/2,'--k')
for n = 300:305
	cs = [zeros(1,n) 1];
	plot(xx, 0.5*focksum(cs, xx)+1/2+n, '-k')
end
ylim([300 306])

figure
xx = 1413.5:0.01:1414.5;
plot(xx, 0.5*fockstate(1e6, xx), '-k')
title 'Millionth eigenstate at turn around point'

xr = 40*rand(1,20) - 20;
xx = -20:0.1:20;

ns = [0:19 20:20:300];
A = fockstate(300,xr,'matrix');
A = A(:,ns+1);

for i = 1:length(ns)
	n = ns(i);
	tic, hh = hermiteF(n,xr); thh(i) = toc;
	cs = [zeros(1,n) 1];
	tic, fs = focksum(cs,xr); tfs(i) = toc;
	assert(all(fs == focksum(cs.',xr)), ...
		'Row and column vectors of coefficients give different sums')
	err(i) = norm(fs-hh, Inf) / norm(hh, Inf);
	sumerr(i) = norm(fs(:)-A(:,i), inf);
end

figure
spy(sumerr)
xlabel n
title 'Discrepencies between fockstate matrix and focksum'

figure
semilogy(ns, err, '.k')
ylabel 'relative error in H_n'
title 'Comparision of focksum with hermiteH'

figure
semilogy(ns, thh, ':k', ns, tfs, '-k')
xlabel n, ylabel t/s, legend hermiteH focksum
title 'Run times'

% test matrix cs
xr = 8*rand(1,8) - 4;
cs = kron(eye(2),ones(1,length(xr)));
qs = focksum(cs, [xr xr]);
assert(all(qs(1:length(xr)) == focksum([1], xr)))
assert(all(qs(length(xr)+1:end) == focksum([0 1], xr)))

function y = hermiteF(n,x)
	y = pi^(-1/4)*exp(-x.^2/2).*hermiteH(n,x) /sqrt(2^n*factorial(n));
end
