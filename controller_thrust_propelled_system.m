(* ::Package:: *)

Get[FileNameJoin[{NotebookDirectory[],"bounded_double_integrator_cw.m"}]];
Get[FileNameJoin[{NotebookDirectory[],"cascaded_attitude_control.m"}]];

ThrustPropelledController[
	t_,(*time instant*)
	{p1_,p2_,p3_,v1_,v2_,v3_,n1_,n2_,n3_,\[Omega]1_,\[Omega]2_,\[Omega]3_},(*state of thrust propelled system*)
	gravity_,(*time-varying gravity*)
	BoundedDoubleIntegratorGains_,(*gains of double integrator*)
	{AttitudeCascadeStrategy_,AttitudeControllerGains_,methodDerivative_,strategyOptions_}(*either "PD" or "BackStepping"*)
]:=
Module[{
	p={p1,p2,p3},(*position*)
	v={v1,v2,v3},(*velocity*)
	pv={p1,p2,p3,v1,v2,v3},(*Position and velocity*)

	n={n1,n2,n3},
	\[Omega]={\[Omega]1,\[Omega]2,\[Omega]3},

	(*variables related to the double integrator control*)
	BDI,uDI,DuDI,DDuDI, VDI,DVDI,DDVDI,WDI,

	(**)
	U3d,D1U3d,D2U3d,D1D1U3d,D2D1U3d,D1D2U3d,D2D2U3d,
	XDIcomplete,D1XDIcomplete,D2XDIcomplete,
	
	DVDIEX,
	
	(*control law for thrust and torque*)
	Tcl,\[Tau]cl,ucl,

	V,W,X,Ex,

	(*BoundedDoubleIntegratorGains={kp\[Rule] 0.2,kv\[Rule] 0.3,\[Sigma]p\[Rule] 0.4,\[Sigma]v\[Rule] 0.5,\[Beta]\[Rule] 0.6}*)
	(*AttitudeControllerGains={kp\[Rule] 0.7,kd\[Rule] 1.1,\[Beta]\[Rule] 0.1};*)
	(*or*)
	(*AttitudeControllerGains={kp\[Rule] 0.7,kd\[Rule] 1.1, kVp\[Rule] 100, kVd\[Rule] 200};*)

	(*Skew symmetric matrix in R3*)
	Skew={{0,-#[[3]],#[[2]]},{#[[3]],0,-#[[1]]},{-#[[2]],#[[1]],0}}&,
	OP=IdentityMatrix[3]-KroneckerProduct[#,#]&

},

	BDI = BoundedDoubleIntegrator[#,BoundedDoubleIntegratorGains]&;
	uDI=BDI[#][[1]][[1]]&;
	DuDI=BDI[#][[1]][[2]]&;
	DDuDI=BDI[#][[1]][[3]]&;

	(*make these functions of Real cross Reals6: even though it does not depend on time*)
	VDI=BDI[#2][[2]][[1]]&;
	DVDI=BDI[#2][[2]][[2]]&;
	DDVDI=BDI[#2][[2]][[3]]&;
	WDI=BDI[#2][[3]]&;

	DVDIEX = DVDI[#t, #pv].Join[ConstantArray[0,{3,3}],Skew[#n]]&[<|"t"-> #1,"pv"->#2[[1;;6]], "n"->#2[[7;;9]] |>]&;

	U3d = uDI[#2] + gravity[#1]&;
	D1U3d = gravity'[#1]&;
	D2U3d = DuDI[#2]&;
	D1D1U3d = gravity''[#1]&;
	D2D1U3d = ConstantArray[0,{3,6}]&;
	D1D2U3d = ConstantArray[0,{3,6}]&;
	D2D2U3d = DDuDI[#2]&;
	
	Tcl = U3d[t,pv].n;

	XDIcomplete = Join[
		#v, 
		(uDI[Join[#p,#v]] + gravity[#t]).#n #n -gravity[#t]
	]&[<|"t"->#1, "p"->#2[[1;;3]], "v"->#2[[4;;6]], "n"->#2[[7;;9]] |>]&;

	D1XDIcomplete = Join[
		{0,0,0}, 
		(gravity'[#t]).#n #n -gravity'[#t]
	]&[<|"t"->#1, "p"->#2[[1;;3]], "v"->#2[[4;;6]], "n"->#2[[7;;9]] |>]&;
	
	D2XDIcomplete = Join[
		Join[ConstantArray[0,{3,3}],IdentityMatrix[3], ConstantArray[0,{3,3}],2],
		Join[
			KroneckerProduct[#n,#n].DuDI[Join[#p,#v]], (*derivative w.r.t. p and v*)
			(uDI[Join[#p,#v]] + gravity[#t]).#n IdentityMatrix[3] + KroneckerProduct[#n,(uDI[Join[#p,#v]] + gravity[#t])], (*derivative w.r.t. p and n*)
			2
		]
	]&[<|"t"->#1, "p"->#2[[1;;3]], "v"->#2[[4;;6]], "n"->#2[[7;;9]] |>]&;
	
	If[AttitudeCascadeStrategy=="PD",
		Module[{aux},
			aux = TorqueAttitudeControlPD[
					t,{p1,p2,p3,v1,v2,v3},Join[n,\[Omega]],
					{U3d,D1U3d,D2U3d,D1D1U3d,D2D1U3d,D1D2U3d,D2D2U3d},
					{XDIcomplete,D1XDIcomplete,D2XDIcomplete},
					{VDI,WDI,DVDIEX},
					AttitudeControllerGains,
					methodDerivative,
					"off"
				];
				
			\[Tau]cl = aux[[1]];
			V = aux[[2]];
			W = aux[[3]];
			X = aux[[4]];
			
		];
	];
	
	If[AttitudeCascadeStrategy=="BackStepping",
	
		Module[{D1DVDIEX,D2DVDIEX},
			
			D1DVDIEX = {0,0,0}&[<|"t"-> #1,"pv"->#2[[1;;6]], "n"->#2[[7;;9]] |>]&;
			D2DVDIEX = Join[
						-Skew[#n].(DDVDI[#t, #pv][[4;;6,;;]]),  (*derivative w.r.t. p and v*)
						Skew[DVDI[#t, #pv][[4;;6]]],  (*derivative w.r.t. n*)
						2
					] &[<|"t"-> #1,"pv"->#2[[1;;6]], "n"->#2[[7;;9]] |>]&;
			
			Module[{aux},
				aux = TorqueBacksteppingController[
						t,{p1,p2,p3,v1,v2,v3},Join[n,\[Omega]],
						{U3d,D1U3d,D2U3d,D1D1U3d,D2D1U3d,D1D2U3d,D2D2U3d},
						{XDIcomplete,D1XDIcomplete,D2XDIcomplete},
						{VDI,WDI,DVDIEX,D1DVDIEX,D2DVDIEX},
						AttitudeControllerGains,
						methodDerivative,
						strategyOptions,
						"off"
					];
				
				\[Tau]cl = aux[[2]][[1]];
				V = aux[[2]][[2]];
				W = aux[[2]][[3]];
				X = aux[[2]][[4]];	
			];
			
		];
	];

	ucl={Tcl,\[Tau]cl};
	{ucl,V,W,X}
]
