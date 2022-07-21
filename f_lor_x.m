function dfx = f_lor_x(t,x,y,z,ModelParams)
% coup=ModelParams.Adj+[0 ModelParams.amp*sin(ModelParams.fr*t);ModelParams.amp*cos(ModelParams.fr*t) 0];
% dfx= -1*ModelParams.a'.*x + ModelParams.a'.*(coup*y);
dfx= -1*ModelParams.a'.*x + ModelParams.a'.*(ModelParams.Adj*y);%+ModelParams.Adj*x;


