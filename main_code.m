clear
rng('shuffle')

%% Basic variables for the setting  
ModelParams=struct();
ModelParams.tau = 0.02; %% time step
ModelParams.nstep = 100000; % number of time steps to generate
ModelParams.N =100;  %number of Oscillators


%% Lorenz Model parameters and their heterogeneity

ModelParams.a=10.0+0.0*rand(1,ModelParams.N);
ModelParams.b=28+0.0*rand(1,ModelParams.N);
ModelParams.c=8/3+0.0*rand(1,ModelParams.N);
ModelParams.NoiseD=10^(-3);%Noise strength


%% Reading the image file and creating the Adjacency matrix from that
Aim=double(imread('false.png'));
Aim=Aim(:,:,1);
ModelParams.Adj=abs((Aim-mean(Aim,'all'))/(3*max(Aim,[],'all')));
hvsd = @(xx) [0.0*(xx == 0) + (xx > 0)];%step function
ModelParams.Adj(:,:)=ModelParams.Adj(:,:).*hvsd(abs(ModelParams.Adj(:,:))-0.2*ones(size(ModelParams.Adj(:,:))));%filtering elements smaller than a threshold

%% Initial conditions for Lorenz model
init_cond=struct();
init_cond.x = 0.1*(-1+2*rand(1,ModelParams.N)); %random initial condition for x(i)
init_cond.y = 0.1*(-1+2*rand(1,ModelParams.N)); %random initial condition for y(i)
init_cond.z = 0.1*(-1+2*rand(1,ModelParams.N)); %random initial condition for z(i)

rng('shuffle')

data = transpose(rk_lor_solve(init_cond,ModelParams));%Data from the simulation of the Lorenz Model
data_new=data(:,20000:end);%Discard initial transients


ModelParams.N=3*ModelParams.N;
measured_vars = 1:1:ModelParams.N;
num_measured = length(measured_vars);
measurements = data(measured_vars, :);% + z;
%% Simulation of the Lorenz model is done

%% Construction of the reservoir, specification of parameters
resparams=struct();
[num_inputs,~] = size(measurements);
resparams.radius = 0.9; % spectral radius
resparams.degree = 3; % connection degree
approx_res_size = 3000; % reservoir size
resparams.N = floor(approx_res_size/num_inputs)*num_inputs; % actual reservoir size divisible by number of inputs
resparams.sigma = 0.1; % input weight scaling
resparams.train_length = 30000; % number of points used to train
resparams.num_inputs = num_inputs; 
resparams.predict_length = 50; % number of predictions after training
resparams.predict_length_max=resparams.predict_length;
resparams.beta = 0.0001; %regularization parameter

%% Training the reservoir 
[~, w_out, A, win,r] = train_reservoir(resparams, data_new(:, 1:resparams.train_length));%Training
% w_out = reservoir-to-prediction coupling matrix
% win = input-to-reservoir coupling matrix
% A = Reservoir connectivity matrix
% r= Reservoir state vector at different times
%% Inferring the Jacobian, please refer to our paper for the equations used

conn1=zeros(size(win,2));%Matrix to store the Jacobian elements

B=A+win*w_out;

av_length=1000;%Length of time for averaging the estimating Jacobian

for it=resparams.train_length-av_length:resparams.train_length
    
%     disp(it)
    
    x=r(:,it);
    A2=B*x;
    mat1=zeros(size(w_out));
    for i1=1:resparams.N
        mat1(:,i1)=w_out(:,i1)*(sech(A2(i1)))^2;

    end
     conn1=conn1+abs(mat1*(win+A*pinv(w_out)));%Estimated Jacobian matrix at each time step (we take mod before time averaging)

end

conn1=conn1/(ModelParams.tau*av_length); %time averaging

%% Extract the X-Y sector the full 3N*3N  Jacobian
xend=ModelParams.N-2;
yend=ModelParams.N-1;
diffxy1(:,:)=conn1(1:3:xend,2:3:yend); %% The X-Y connectivity Jacobian


%% plot actual and inferred Jacobian side-by-side
figure
s1 = (1:1:ModelParams.N/3);
s2 = (1:1:ModelParams.N/3);

subplot(1,2,1)
imagesc(s1,s2,10*ModelParams.Adj)%Actual Jacobian, has a factor of 10 coming from the parameter a in the Lorenz model

caxis([0 max(max(ModelParams.Adj*10,diffxy1),[],'all')])
colorbar

colormap('gray')
pbaspect([1 1 1])
title('$L=100$','Interpreter', 'Latex','fontsize',40)

subplot(1,2,2)
imagesc(s1,s2,diffxy1)%Inferred Jacobian

caxis([0 max(max(ModelParams.Adj*10,diffxy1),[],'all')])
colorbar

colormap('gray')
pbaspect([1 1 1])
