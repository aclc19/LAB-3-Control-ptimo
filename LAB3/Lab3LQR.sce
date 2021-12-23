for i=1:length(scs_m.objs)
    if typeof(scs_m.objs(i))=="Block" & scs_m.objs(i).gui=="SUPER_f" then
        scs_m = scs_m.objs(i).model.rpar;
        break;
    end
end
L=[15;15]; // punto de trabajo
U=[0];
sys = lincos(scs_m ,L,U)

A=sys.A
B=sys.B
C=sys.C
D=sys.D

P=syslin('c',A,B,C,D)    //The plant (continuous-time)
Q_xx=diag([1 1]); //Weights on states
R_uu   = 0.1; //Weight on input peso de la entrada
Kc=lqr(P,Q_xx,R_uu);

M=syslin('c',(A+B*Kc),B,C,D);
lam=spec(A+B*Kc)
G=ss2tf(M)// funcion de transferencia 
scf(1)
bode (G,0.01 ,100)
title('comportamiento en frecuencia - Lazo abierto')
GC=G/(1+G)
scf(2)
bode (GC,0.01 ,100)
title('comportamiento en frecuencia - Lazo cerrado')
scf(3)
dt=0.1;
t=0:dt:200;
y=csim('step',t,M,L);
clf;plot(t',y');xlabel(_("time (s)"))
L=legend(["$dy_1$","$dy_2$"]);L.font_size=4;
