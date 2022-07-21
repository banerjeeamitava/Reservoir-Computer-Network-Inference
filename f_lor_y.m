function dfy = f_lor_y(t,x,y,z,ModelParams)

dfy = ModelParams.b'.*x - y -x.*z;%+ModelParams.Adj*y;


