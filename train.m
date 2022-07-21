function w_out = train(params, states, data)

beta = params.beta;
% rng(2,'twister');
 rng(5,'twister');

idenmat = beta*speye(params.N);

% states(2:2:params.N,:) = states(2:2:params.N,:).^2;
w_out = data*transpose(states)*pinv(states*transpose(states)+idenmat);
