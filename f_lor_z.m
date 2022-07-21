function dfz = f_lor_z(t,x,y,z,ModelParams)

dfz = -ModelParams.c'.*z + x.*y;%+ModelParams.Adj*z;

