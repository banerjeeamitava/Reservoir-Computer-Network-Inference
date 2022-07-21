function A = generate_reservoir(size, radius, degree)
% rng(1,'twister');
rng(3,'twister');

sparsity = degree/size;

A = sprand(size, size, sparsity);

e = max(abs(eigs(A)));

A = (A./e).*radius;
