function [x, wout, A, win,states] = train_reservoir(resparams, data)

A = generate_reservoir(resparams.N, resparams.radius, resparams.degree);
% A=zeros(resparams.N);%Setting A to zero
q = resparams.N/resparams.num_inputs;
win = zeros(resparams.N, resparams.num_inputs);
for i=1:resparams.num_inputs
%     rng(i,'twister')
     rng(i+51,'twister')

    ip = resparams.sigma*(-1 + 2*rand(q,1));
    win((i-1)*q+1:i*q,i) = ip;
end
%win = resparams.sigma*(-1 + 2*rand(resparams.N, resparams.num_inputs));

states = reservoir_layer(A, win, data, resparams);

wout = train(resparams, states, data(:,1:resparams.train_length));

x = states(:,end);

%[output,x] = predict(A,win,resparams,states(:,end), wout);
% %%
% 
% figure()
% subplot(2,1,1)
% imagesc(output)
% subplot(2,1,2)
% imagesc(data(:,resparams.train_length+1:resparams.train_length + resparams.predict_length) - output)
