(* ::Package:: *)

(*saturation function: second argument is to be considered as a (fixed) parameter*)
d0\[Sigma]c[xx_,\[Sigma]\[Sigma]_]:=(\[Sigma]\[Sigma] xx)/Sqrt[\[Sigma]\[Sigma]^2+xx xx];
d1\[Sigma]c[xx_,\[Sigma]\[Sigma]_]:=Simplify[Evaluate[D[d0\[Sigma]c[xxx,\[Sigma]\[Sigma]],xxx]]]/.{xxx-> xx};
d2\[Sigma]c[xx_,\[Sigma]\[Sigma]_]:=Simplify[Evaluate[D[d1\[Sigma]c[xxx,\[Sigma]\[Sigma]],xxx]]]/.{xxx-> xx};

\[Sigma][xx_,\[Sigma]\[Sigma]_]:={d0\[Sigma]c[xx[[1]],\[Sigma]\[Sigma]],d0\[Sigma]c[xx[[2]],\[Sigma]\[Sigma]],d0\[Sigma]c[xx[[3]],\[Sigma]\[Sigma]]}

(*1st derivative saturation function*)
D\[Sigma][xx_,\[Sigma]\[Sigma]_]:=DiagonalMatrix[{d1\[Sigma]c[xx[[1]],\[Sigma]\[Sigma]],d1\[Sigma]c[xx[[2]],\[Sigma]\[Sigma]],d1\[Sigma]c[xx[[3]],\[Sigma]\[Sigma]]}]

(*2nd derivative saturation function*)
DD\[Sigma][xx_,\[Sigma]\[Sigma]_]:={
	DiagonalMatrix[{d2\[Sigma]c[xx[[1]],\[Sigma]\[Sigma]],0,0}],
	DiagonalMatrix[{0,d2\[Sigma]c[xx[[2]],\[Sigma]\[Sigma]],0}],
	DiagonalMatrix[{0,0,d2\[Sigma]c[xx[[3]],\[Sigma]\[Sigma]]}]
}

V\[Sigma][p_,v_,\[Sigma]p_,\[Sigma]v_,kp_,kv_,\[Beta]_] := kp (\[Sigma]p Sqrt[p^2+\[Sigma]p^2]-\[Sigma]p^2)  +\[Beta] d0\[Sigma]c[p,\[Sigma]p]d0\[Sigma]c[v,\[Sigma]v]+1/2 v^2;
DV\[Sigma][p_,v_,\[Sigma]p_,\[Sigma]v_,kp_,kv_,\[Beta]_]:={
	kp d0\[Sigma]c[p,\[Sigma]p] +\[Beta] d1\[Sigma]c[p,\[Sigma]p]d0\[Sigma]c[v,\[Sigma]v],
	\[Beta] d1\[Sigma]c[v,\[Sigma]v]d0\[Sigma]c[p,\[Sigma]p]+v
};

DDV\[Sigma][p_,v_,\[Sigma]p_,\[Sigma]v_,kp_,kv_,\[Beta]_]:={
	{kp d1\[Sigma]c[p,\[Sigma]p] +\[Beta] d2\[Sigma]c[p,\[Sigma]p] d0\[Sigma]c[v,\[Sigma]v],\[Beta] d1\[Sigma]c[p,\[Sigma]p]d1\[Sigma]c[v,\[Sigma]v]},
	{\[Beta] d1\[Sigma]c[v,\[Sigma]v] d1\[Sigma]c[p,\[Sigma]p] ,\[Beta] d2\[Sigma]c[v,\[Sigma]v] d0\[Sigma]c[p,\[Sigma]p]+1}
};

Wmatrix\[Sigma][p_,v_,\[Sigma]p_,\[Sigma]v_,kp_,kv_,\[Beta]_]:=\[Sigma]v/Sqrt[v^2+\[Sigma]v^2] {
	{kp \[Beta] Sqrt[v^2+\[Sigma]v^2]/\[Sigma]v d1\[Sigma]c[v,\[Sigma]v],1/2 \[Beta] kv d1\[Sigma]c[v,\[Sigma]v]},
	{1/2 \[Beta] kv d1\[Sigma]c[v,\[Sigma]v],kv- \[Beta] d1\[Sigma]c[p,\[Sigma]p]}
};

W\[Sigma][p_,v_,\[Sigma]p_,\[Sigma]v_,kp_,kv_,\[Beta]_]:=-{d0\[Sigma]c[p,\[Sigma]p],v}.Wmatrix\[Sigma][p,v,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]].{d0\[Sigma]c[p,\[Sigma]p],v};

BoundedDoubleIntegrator[
{px_,py_,pz_,vx_,vy_,vz_},(*Position and velocity*)
ControllerGains_
]:=
Module[{

(*gains*)
kp =kp/.ControllerGains,
kv=kv/.ControllerGains,
\[Sigma]p =\[Sigma]p/.ControllerGains,
\[Sigma]v=\[Sigma]v/.ControllerGains,
\[Beta]=\[Beta]/.ControllerGains,

u,Du,DDu,
V,DV,DDV,
W,Wmatrix,
P,Z

},
(*controller for double integrator*)
u=-kp \[Sigma][{px,py,pz},\[Sigma]p]-kv \[Sigma][{vx,vy,vz},\[Sigma]v];

(*array (3,6)*)
Du=Join[-kp D\[Sigma][{px,py,pz},\[Sigma]p],-kv D\[Sigma][{vx,vy,vz},\[Sigma]v],2];
(*array (3,6,6)*)
DDu=Join[
Join[-kp DD\[Sigma][{px,py,pz},\[Sigma]p],ConstantArray[0,{3,3,3}],2],(*first, derivative w.r.t. position*)
Join[ConstantArray[0,{3,3,3}],-kv DD\[Sigma][{vx,vy,vz},\[Sigma]v],2],(*second, derivative w.r.t. velocity*)
3
];

V = V\[Sigma][px,vx,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]]+
V\[Sigma][py,vy,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]]+
V\[Sigma][pz,vz,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]];

(**)
DV={
DV\[Sigma][px,vx,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[1]],
DV\[Sigma][py,vy,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[1]],
DV\[Sigma][pz,vz,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[1]],
DV\[Sigma][px,vx,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[2]],
DV\[Sigma][py,vy,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[2]],
DV\[Sigma][pz,vz,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[2]]
};

P={
{1,0,0,0,0,0},{0,0,0,1,0,0},
{0,1,0,0,0,0},{0,0,0,0,1,0},
{0,0,1,0,0,0},{0,0,0,0,0,1}
};
Z=ConstantArray[0,{2,2}];
(*DDV=ConstantArray[0,{6,6}];*)
(*DDV=P.Join[
Join[DDV\[Sigma][px,vx,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]],Z,Z,2],
Join[Z,DDV\[Sigma][py,vy,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]],Z,2],
Join[Z,Z,DDV\[Sigma][pz,vz,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]],2]
].P;
*)
DDV={
{DDV\[Sigma][px,vx,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[1,1]],0,0,DDV\[Sigma][px,vx,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[1,2]],0,0},
{0,DDV\[Sigma][py,vy,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[1,1]],0,0,DDV\[Sigma][py,vy,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[1,2]],0},
{0,0,DDV\[Sigma][pz,vz,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[1,1]],0,0,DDV\[Sigma][pz,vz,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[1,2]]},
{DDV\[Sigma][px,vx,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[2,1]],0,0,DDV\[Sigma][px,vx,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[2,2]],0,0},
{0,DDV\[Sigma][py,vy,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[2,1]],0,0,DDV\[Sigma][py,vy,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[2,2]],0},
{0,0,DDV\[Sigma][pz,vz,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[2,1]],0,0,DDV\[Sigma][pz,vz,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]][[2,2]]}
};

W = W\[Sigma][px,vx,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]]+
W\[Sigma][py,vy,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]]+
W\[Sigma][pz,vz,\[Sigma]p,\[Sigma]v,kp,kv,\[Beta]];

{{u,Du,DDu},{V,DV,DDV},W}
]
