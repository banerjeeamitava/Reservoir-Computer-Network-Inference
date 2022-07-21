function soln = rk_lor_solve(init_cond,ModelParams)
rng(13,'twister');

%initialize structures
t=zeros(ModelParams.nstep+1,1);
t(1) = 0;
x = zeros(ModelParams.N,ModelParams.nstep+1);
y = zeros(ModelParams.N,ModelParams.nstep+1);
z = zeros(ModelParams.N,ModelParams.nstep+1);
x(:,1) = init_cond.x;
y(:,1) = init_cond.y;
z(:,1) = init_cond.z;

for i=1:ModelParams.nstep
    %update time
    t(i+1)=t(i)+ModelParams.tau;
    %update vector
    %y is handle of solution vector
    [x(:,i+1),y(:,i+1),z(:,i+1)]= rk4_lor(@f_lor_x,@f_lor_y,@f_lor_z, x(:,i)+sqrt(2*ModelParams.NoiseD*ModelParams.tau)*randn(size(x(:,i))),y(:,i)+sqrt(2*ModelParams.NoiseD*ModelParams.tau)*randn(size(y(:,i))),z(:,i)+sqrt(2*ModelParams.NoiseD*ModelParams.tau)*randn(size(z(:,i))), t(i), ModelParams);
end
yy=zeros(3*ModelParams.N,ModelParams.nstep+1);

yy(1:3:3*ModelParams.N-2,:)=x;
yy(2:3:3*ModelParams.N-1,:)=y;
yy(3:3:3*ModelParams.N,:)=z;

soln = yy';

