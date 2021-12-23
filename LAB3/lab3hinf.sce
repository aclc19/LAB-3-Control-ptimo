// HROV linearized model
// Author: juan C. Cutipa-Luque
// Source:
// x=[u v r x y apsi] 
// tau= [tau1 tau2]
// y= [u v r x y apsi] 
clf();         // close current figure
clear          // clear all pasta variables
xdel(winsid()) // close all windows
A=[-0.0656    0.0000    0.0000         0         0         0
    0.0656   -0.0656        0         0         0         0
    0.0000    0.0000    0.0000        0         0          0
    0.0000         0         0         0         0         0
    0.0000         0         0         0         0         0
         0         0         0        0         0          0];
 
B=[ 0.2964    0.0000
    0.0000    0.0000
   -0.0000    0.0000
         0         0
         0         0
         0         0];

C =[ 1     0     0     0     0     0
     0     1     0     0     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0];

D=[  0     0
     0     0
     0     0
     0     0
     0     0
     0     0];

// Reducing the model
// X=[u, v, r]', U=[tau_u, tau_r]', Y=[u r]'
ap=A(1:3,1:3);  
bp=B(1:3,:);
cp=[C(1,1:3);
    C(3,1:3)];
dp=[D(1,1:2);
    D(3,1:2)];
// Controllability and Observability
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(ap,bp)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap, cp)
rankO=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a

/*             Plot singular values of LTI the model                      */
G = syslin('c', ap, bp, cp, dp);//Planta modelo LTI
poles=spec(ap)
tr  = trzeros(G)

w = logspace(-3,3);
sv = svplot(G,w);

scf(0);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");
ZG=ss2tf(G,w)
/*                                 Scaling                                 */
d2r=%pi/180;
r2d=180/%pi;

su = diag( [1/200,1/(200*0.8)] );   // scaling input
                                    // tau_u_max=300 [N] tau_r_max=100 [Nm]
sx = diag([1/1,1/0.5,1/(50*d2r)]);  // scaling state 

sy = diag([1/1,1/(50*d2r)]);        // scaling output

ap_s = sx*ap*inv(sx)
bp_s = sx*bp*inv(su)
cp_s = sy*cp*inv(sx)
dp_s = sy*dp*inv(su)

//g=minreal(g)
G = syslin('c', ap_s, bp_s, cp_s, dp_s);
w = logspace(-3,3);
sv = svplot(G,w);

scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");


ms=1.7;// 0.3;%1.5;    % guarantee overshot Mp < 6dB = 20*log10(2) 
wbs=0.23;//0.05;%0.23;
ee=1e-3;//1e-4
ki=1; // used to give more accurate adjustment to the cut-off frequency wbs
      // by default set it to 1
//           --------     WT Data    ------------
mt=1.3;//1.00;    % guarantee overshot Mp < 2dB = 20*log10(1.26)
wbt=4.1;//9.1;%4.1;
ee=1e-3;//1e-4

//           --------     WS     ------------

s=poly(0,'s');
ws1=(s/ms+wbs)/(s+wbs*ee),
ws2=ws1;
ws=[ws1,0;0,ws2]
//Ws=syslin('c',ws)
Ws=blockdiag(ws1,ws2)

//           --------     WT     ------------

s=poly(0,'s');
wt1=(s+wbt/mt)/(ee*s+wbt),
wt2=wt1;
wt=[wt1,0;0,wt2]
//Wt=syslin('c',wt)
Wt=blockdiag(wt1,wt2)


//           --------     WR     ------------
s=poly(0,'s');
wr1=s/s,
wr2=wr1;
wr=[wr1,0;0,wr2]



// ------------------ Plot weighting functions

svs = svplot(Ws,w);
svt = svplot(Wt,w);
scf(2);
plot2d("ln", w, [-20*log(svs')/log(10) -20*log(svt')/log(10)])
xgrid(12)
xtitle("Singular values plot inv(Ws) and inv(Wt)","Frequency (rad/s)", "Amplitude (dB)");


[P,r]=augment(G,'ST');
//[P,r]=augment(g,'SRT');
P = blockdiag(Ws,Wt,eye(G))*P;
//P=minreal(P);

//trick to tackle when "D12 is not full rank"
P.D(1,3)=0.0001;
P.D(2,4)=0.0001;

r=[2,2]
romin=0.0001
romax=2000
nmax=100;
//[K,ro]=h_inf(P,r,romin,romax,nmax)
// alternatives
//[AK,BK,CK,DK,(RCOND)] = hinf(P.A,P.B,P.C,P.D,2,2,4)
K = ccontrg(P, [2,2], 1.3) // this is good for me and for this system

// -------------- Analysis of the Feeedback Control System

[Se,Re,Te]=sensi(G,K) // S=(I+GK)^-1, T=I-S=GK(I+GK)^-1

// ------------------  weighting functions

svS = svplot(Se,w);
svT = svplot(Te,w);

// --------------- Open loop and Closed loop analysis --------------

sysOL=G*K
sysCL=G*K*inv(eye(2,2)+G*K)

//LA=ss2tf(sysOL,w) Lazo abierto
OL=wt(1,1)*ZG(1,1);
//CL=ss2tf(sysCL,w) Lazo cerrado
CL=ws(1,1)*ZG(1,1);
//Gráficas XCOS

//SISTEMA CONTROLADO
scf(3);
t=0:0.1:100;
SC=csim('step',t,ZG(1,1))
plot(t',SC')
xtitle('RPTA EN EL TIEMPO');xgrid;
//legend('abierto','cerrado')
// Time responses in XCOS
