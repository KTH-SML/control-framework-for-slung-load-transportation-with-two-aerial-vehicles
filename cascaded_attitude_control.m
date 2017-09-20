(* ::Package:: *)

Function1st2ndDerivatives[
t_, (*time instant*)
x_, (*state before cascade*)
{n1_,n2_,n3_,\[Omega]1_,\[Omega]2_,\[Omega]3_}, (*attitude unit vector and angular velocity*)
{u_,D1u_,D2u_,D1D1u_,D2D1u_,D1D2u_,D2D2u_}, (*vector whose direction attitude unit vector needs to track*)
{X_,D1X_,D2X_},(*dynamics of state x*)
method_
]:=
Module[{

	(*unit vector and angular velocity*)
	n={n1,n2,n3},
	\[Omega]={\[Omega]1,\[Omega]2,\[Omega]3},
	
	(*useful variables*)
	xd0,xd1,xd2, (*state and derivatives*)
	ud0,ud1,ud2, (*vector u and derivatives*)
	yd0, yd1, (*yd0 = [x,n] = [xd0,n]*)
	
	(*Skew symmetric matrix in R3*)
	Skew={{0,-#[[3]],#[[2]]},{#[[3]],0,-#[[1]]},{-#[[2]],#[[1]],0}}&
},

	xd0=x;
	yd0 = Join[xd0,n];
	ud0 =u[t,xd0];

	xd1=X[t,yd0];
	yd1=Join[xd1,Skew[\[Omega]].n];
	
	If[method[[1]]=="IgnoreInput2ndDerivative", 
		ud1={0,0,0};
		ud2={0,0,0};
	];
	If[method[[1]]=="Numerical", 
		Module[{
			dt = method[[2]],
			du = D1u[#1,#2[[1;;12]]]+D2u[#1,#2[[1;;12]]].X[t,#2]&
			},
			ud1=(u[t+dt, xd0 + dt xd1] - u[t,xd0])/dt;
			ud2 = (du[t+dt, yd0 + dt yd1] - du[t, yd0])/dt;
		];
	];
	If[method[[1]]=="NumericalSloppyAlternative", 
		Module[{dt = method[[2]]},
			ud1=(u[t+dt, xd0 + dt xd1] - u[t,xd0])/dt;
			(*central second derivative*)		
			ud2=(u[t+dt, xd0 + dt xd1] - 2 u[t, xd0] + u[t - dt, xd0 - dt xd1])/dt^2;
		];
	];
	If[method[[1]]=="Precise", 
		ud1 = D1u[t,xd0]+D2u[t,xd0].xd1;
		xd2 = D1X[t,yd0]+D2X[t,yd0].yd1;
		ud2 = D1D1u[t,xd0] + 2D2D1u[t,xd0].xd1 + (D2D2u[t,xd0].xd1).xd1 + D2u[t,xd0].xd2;
	];
	
	{ud0,ud1,ud2}
]


TorqueAttitudeControlPD[
t_, (*time instant*)
x_, (*state before cascade*)
{n1_,n2_,n3_,\[Omega]1_,\[Omega]2_,\[Omega]3_}, (*attitude unit vector and angular velocity*)
{u_,D1u_,D2u_,D1D1u_,D2D1u_,D1D2u_,D2D2u_}, (*vector whose direction attitude unit vector needs to track*)
{X_,D1X_,D2X_},(*dynamics of state x*)
{V_,W_,DVEx_},(*Lyapunov and error gradient of Lyapunov*)
AttitudeControllerGains_,
methodDerivative_,
simpleTracking_
]:=
Module[{

	(*unit vector and angular velocity*)
	n={n1,n2,n3},
	\[Omega]={\[Omega]1,\[Omega]2,\[Omega]3},
	
	(*useful variables*)
	ud0,ud1,ud2, (*vector u and derivatives*)
	
	nd,\[Omega]d,\[Tau]d, (*desired unit vector, angular velocity and torque*)
	en,e\[Omega], (*error unit vector and angular velocity*)
	\[Tau]ff,\[Tau]cl, (*feedforward and final control law for torque*)

	V\[Theta],W\[Theta], (*lyapunov and derivative associated to unit vector tracking problem*)
	VV,WW,XX, (*final lyapunov, its derivative and vector field*)

	(*AttitudeControllerGains={kp\[Rule] 0.7,kd\[Rule] 1.1, \[Beta]\[Rule] 0.81};*)
	kp=kp/.AttitudeControllerGains,
	kd=kd/.AttitudeControllerGains,
	\[Beta]=\[Beta]/.AttitudeControllerGains,

	(*Skew symmetric matrix in R3*)
	Skew={{0,-#[[3]],#[[2]]},{#[[3]],0,-#[[1]]},{-#[[2]],#[[1]],0}}&,
	(*orthogonal projection operator*)
	OP=IdentityMatrix[3]-KroneckerProduct[#,#]&

},
	
	If[simpleTracking=="on",
		(*if there is no state x that the vector u depends on*)
		(*recall objective: n \[Rule] u/Norm[u]*)
		ud0 = u[t];
		ud1 = u'[t];
		ud2 = u''[t];
		,
		Module[{uds},	
			uds = Function1st2ndDerivatives[t,x,Join[n,\[Omega]],
				{u,D1u,D2u,D1D1u,D2D1u,D1D2u,D2D2u},
				{X,D1X,D2X},
				methodDerivative];	
			ud0 = uds[[1]];
			ud1 = uds[[2]];
			ud2 = uds[[3]];
		];
	];

	nd=#/Sqrt[#.#]&[ud0];
	\[Omega]d=Skew[nd].(#2/Sqrt[#1.#1])&[ud0,ud1];
	\[Tau]d=Skew[nd].(#3/Sqrt[#1.#1]-2 #2.#1/#1.#1 #2/Sqrt[#1.#1])&[ud0,ud1,ud2];

	en=Skew[nd].n;
	e\[Omega]=Skew[n].(\[Omega]-\[Omega]d);

	\[Tau]ff=OP[n].\[Tau]d+Skew[n].\[Omega](n.\[Omega]d);
	\[Tau]cl = \[Tau]ff -kp en-kd OP[n].(\[Omega]-\[Omega]d);

	V\[Theta] = kp(1-n.nd)-\[Beta]  n.Skew[nd].(\[Omega]-\[Omega]d)+1/2 (#.#)&[e\[Omega]];
	W\[Theta] = -kp \[Beta] (#.#)&[en]-\[Beta] e\[Omega] .(n.\[Omega]d IdentityMatrix[3] + kd Skew[n]).en -(kd-\[Beta] n.nd) (#.#)&[e\[Omega]];


	If[simpleTracking=="on",
	
		VV = V\[Theta];
		WW = W\[Theta];

		XX=
		Join[
			Skew[\[Omega]].n,
			\[Tau]cl (*=OP[n].\[Tau]cl*)
		];
		,

		VV = V[t,x] + V\[Theta];
		WW = W[t,x]+ DVEx[t,Join[x,n]].Skew[n].ud0 + W\[Theta];

		XX=
		Join[
			X[t,Join[x,n]],
			Skew[\[Omega]].n,
			\[Tau]cl (*=OP[n].\[Tau]cl*)
		];

	];
	
	{\[Tau]cl,VV,WW,XX}
]


TorqueBacksteppingController[
t_, (*time instant*)
x_, (*state before cascade*)
{n1_,n2_,n3_,\[Omega]1_,\[Omega]2_,\[Omega]3_}, (*attitude unit vector and angular velocity*)
{u_,D1u_,D2u_,D1D1u_,D2D1u_,D1D2u_,D2D2u_},
{X_,D1X_,D2X_},
{V_,W_,DVEX_,D1DVEX_,D2DVEX_},
ControllerGains_,
methodDerivative_,
cascadeBackstepping_,
simpleTracking_]:=
Module[{

	(*state*)
	n={n1,n2,n3},
	\[Omega]={\[Omega]1,\[Omega]2,\[Omega]3},
	
	x1 = x, 
	x2 = Join[x,{n1,n2,n3}], 
	(*x3 = Join[x,{n1,n2,n3,\[Omega]1,\[Omega]2,\[Omega]3}]*)
	
	(*part 1*)
	d0u,ncl,X1,
	
	(*part 2*)
	d1u,d0bs,\[Omega]ncl,\[Omega]cl,X2,V2,W2,
	
	(*part 3*)
	d2u,d1bs,\[Tau]\[Omega]ncl,\[Tau]\[Omega]cl,\[Tau]cl,X3,V3,W3,
	
	(*gains*)
	kp = kp/.ControllerGains,
	kd = kd/.ControllerGains,
	kVp = kVp/.ControllerGains,
	kVd = kVd/.ControllerGains,

	(*defining some useful functions*)
	(*Orthogonal Projection Operator, in R3*)
	OP=IdentityMatrix[3]-KroneckerProduct[#, #]&,
	(*Skew symmetric matrix in R3*)
	Skew={{0,-#[[3]],#[[2]]},{#[[3]],0,-#[[1]]},{-#[[2]],#[[1]],0}}&

},

	If[simpleTracking=="on",
		(*if there is no state x that the vector u depends on*)
		(*recall objective: n \[Rule] u/Norm[u]*)
		d0u = u[t];
		d1u = u'[t];
		d2u = u''[t];
		,
		Module[{uds},	
			uds = Function1st2ndDerivatives[t,x,Join[n,\[Omega]],
				{u,D1u,D2u,D1D1u,D2D1u,D1D2u,D2D2u},
				{X,D1X,D2X},
				methodDerivative];	
			d0u = uds[[1]];
			d1u = uds[[2]];
			d2u = uds[[3]];
		];
	];
	
	(*0th step*)
	ncl=d0u/Sqrt[d0u.d0u];
	
	
	If[cascadeBackstepping=="off" || simpleTracking=="on",	
		(*1st step*)
		\[Omega]ncl=Skew[ncl].(d1u/Sqrt[d0u.d0u]);
		\[Omega]cl = \[Omega]ncl - kp Skew[ncl].n;

		If[simpleTracking=="on",
			X2 = Skew[\[Omega]].n;
			V2 = kVp(1-n.ncl);
			W2 = - kVp kp (Skew[n].ncl).(Skew[n].ncl);					
			,
			X1 = X[t,x2];
			X2 = Join[X1,Skew[\[Omega]].n];
			V2 = V[t,x1] + kVp(1-n.ncl);
			W2 = W[t,x1] + DVEX[t,x2].Skew[n].d0u - kVp kp (Skew[n].ncl).(Skew[n].ncl);		
			(*V2 = V[t,x1] + kVp(1-n.ncl);
			W2 = W[t,x1] + DVEX[t,x2].Skew[n].d0u - kVp kp (Skew[n].ncl).(Skew[n].ncl) - kVp (Skew[n].ncl).(\[Omega] - \[Omega]cl);*)
		];
			
		(*2nd step*)
		\[Tau]\[Omega]ncl=Skew[ncl].(d2u/Sqrt[d0u.d0u]-2 d1u.d0u/d0u.d0u d1u/Sqrt[d0u.d0u]);
		\[Tau]\[Omega]cl = \[Tau]\[Omega]ncl-kp Skew[Skew[\[Omega]ncl].ncl].n-kp Skew[ncl].Skew[\[Omega]].n;
		\[Tau]cl = OP[n].(\[Tau]\[Omega]cl+(\[Omega]cl.n)Skew[n].\[Omega]cl-kd(\[Omega]-\[Omega]cl)+kVp/kVd Skew[n].ncl);

		X3 = Join[X2,OP[n].\[Tau]cl];
		V3 = V2 + kVd/2 (Skew[n].(\[Omega]-\[Omega]cl)).(Skew[n].(\[Omega]-\[Omega]cl));
		W3 = W2 - kVd kd (Skew[n].(\[Omega]-\[Omega]cl)).(Skew[n].(\[Omega]-\[Omega]cl));
		
	];
	
	(*the method below will never be implemented if X=Null*)
	If[cascadeBackstepping=="on",	
		(*1st step*)
		\[Omega]ncl=Skew[ncl].(d1u/Sqrt[d0u.d0u]);
		(*backstepping term*)
		d0bs = DVEX[t,x2];
		\[Omega]cl = \[Omega]ncl - kp Skew[ncl].n + Sqrt[d0u.d0u]/kVp d0bs;
		
		X1 = X[t,x2];
		X2 = Join[X1,Skew[\[Omega]].n];
		V2 = V[t,x1] + kVp(1-n.ncl);
		W2 = W[t,x1] - kVp kp (Skew[n].ncl).(Skew[n].ncl);
		
		(*2nd step*)
		\[Tau]\[Omega]ncl=Skew[ncl].(d2u/Sqrt[d0u.d0u]-2 d1u.d0u/d0u.d0u d1u/Sqrt[d0u.d0u]);
		
		d1bs = D1DVEX[t,x2] + D2DVEX[t,x2].X2;
		\[Tau]\[Omega]cl = \[Tau]\[Omega]ncl-kp Skew[Skew[\[Omega]ncl].ncl].n-kp Skew[ncl].Skew[\[Omega]].n + d0u.d1u/Sqrt[d0u.d0u] 1/kVp d0bs + Sqrt[d0u.d0u]/kVp d1bs;
		\[Tau]cl = OP[n].(\[Tau]\[Omega]cl+(\[Omega]cl.n)Skew[n].\[Omega]cl-kd(\[Omega]-\[Omega]cl)+kVp/kVd Skew[n].ncl);
	
		X3 = Join[X2,OP[n].\[Tau]cl];
		V3 = V2 + kVd/2 (Skew[n].(\[Omega]-\[Omega]cl)).(Skew[n].(\[Omega]-\[Omega]cl));
		W3 = W2 - kVd kd (Skew[n].(\[Omega]-\[Omega]cl)).(Skew[n].(\[Omega]-\[Omega]cl));
	];
	
	{
	{\[Omega]cl,V2,W2,X2},(*if you wish control at first level -- angular velocity control*)
	{\[Tau]cl,V3,W3,X3}(*if you wish control at second level -- torque control*)
	}
]
