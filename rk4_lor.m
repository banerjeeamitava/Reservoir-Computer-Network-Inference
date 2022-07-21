function [x,y,z] = rk4_lor(f_lor_x,f_lor_y,f_lor_z, x,y,z, t, ModelParams)
%RK4 Runge-Kutta 4th order integration
%   RK4(F,Y,T,DT,...) integrates the differential equation y' = f(t,y) from
%   time T to time T+DT, where Y is the value of the solution vector at
%   time T. F is a function handle. For a scalar T and a vector Y, F(T,Y)
%   must return a column vector corresponding to f(t,y). Additional
%   arguments will be passed to the function F(T,Y,...).

k1x = feval(f_lor_x, t,        x,y,z          , ModelParams);
k1y=feval(f_lor_y, t,        x,y,z          , ModelParams);
k1z=feval(f_lor_z, t,        x,y,z          , ModelParams);

k2x = feval(f_lor_x, t + ModelParams.tau/2, x + k1x*ModelParams.tau/2,y + k1y*ModelParams.tau/2,z + k1z*ModelParams.tau/2, ModelParams);
k2y = feval(f_lor_y, t + ModelParams.tau/2, x + k1x*ModelParams.tau/2,y + k1y*ModelParams.tau/2,z + k1z*ModelParams.tau/2, ModelParams);
k2z = feval(f_lor_z, t + ModelParams.tau/2, x + k1x*ModelParams.tau/2,y + k1y*ModelParams.tau/2,z + k1z*ModelParams.tau/2, ModelParams);

k3x = feval(f_lor_x, t + ModelParams.tau/2, x + k2x*ModelParams.tau/2, y + k2y*ModelParams.tau/2,z + k2z*ModelParams.tau/2,ModelParams);
k3y = feval(f_lor_y, t + ModelParams.tau/2, x + k2x*ModelParams.tau/2, y + k2y*ModelParams.tau/2,z + k2z*ModelParams.tau/2,ModelParams);
k3z = feval(f_lor_z, t + ModelParams.tau/2, x + k2x*ModelParams.tau/2, y + k2y*ModelParams.tau/2,z + k2z*ModelParams.tau/2,ModelParams);

k4x = feval(f_lor_x, t + ModelParams.tau, x + k3x*ModelParams.tau, y + k3y*ModelParams.tau, z + k3z*ModelParams.tau, ModelParams);
k4y = feval(f_lor_y, t + ModelParams.tau, x + k3x*ModelParams.tau, y + k3y*ModelParams.tau, z + k3z*ModelParams.tau, ModelParams);
k4z = feval(f_lor_z, t + ModelParams.tau, x + k3x*ModelParams.tau, y + k3y*ModelParams.tau, z + k3z*ModelParams.tau, ModelParams);

x = x + ModelParams.tau/6 * (k1x + 2*k2x + 2*k3x + k4x);
y = y + ModelParams.tau/6 * (k1y + 2*k2y + 2*k3y + k4y);
z = z + ModelParams.tau/6 * (k1z + 2*k2z + 2*k3z + k4z);


