(* ::Package:: *)

BeginPackage["LG`"];


LGMag::usage = "LGMag[l, p] Returns image of magnitude of LG mode l, p. Default options Width \[Rule] 500, Height \[Rule] 500, WaveVector \[Rule] 1.0, BeamWaist \[Rule] 10.0, ScalingFactor \[Rule] 1.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False.";
LGPhase::usage = "LGPhase[l, p] Returns image of phase of LG mode l, p. Default options Width \[Rule] 500, Height \[Rule] 500, WaveVector \[Rule] 1.0, BeamWaist \[Rule] 10.0, ScalingFactor \[Rule] 1.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False.";
TwoLGMag::uage = "TwoLGMag[l0, p0, l1, p1] Returns image of magnitude of sum of two LG modes. Default options  Width \[Rule] 500, Height \[Rule] 500, WaveVector \[Rule] 1.0, BeamWaist \[Rule] 10.0, ScalingFactor \[Rule] 1.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False."  ;
TwoLGPhase::usage = "TwoLGPhase[l0, p0, l1, p1] Returns image of phase of sum of two LG modes. Default options  Width \[Rule] 500, Height \[Rule] 500, WaveVector \[Rule] 1.0, BeamWaist \[Rule] 10.0, ScalingFactor \[Rule] 1.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False.";
StartLG::usage = "StartLG[path, binpath] Install WSTP executable. Must be run before anything else. 'path' is the absolute path of repository. path/binpath should resolve to path of the WSTP binary.";
BlazeImage::usage = "BlazeImage[image, angle] Returns 'image' with blazing added on top. 
Default options Frequency \[Rule] 0.5, PhaseOffset \[Rule] 0, Size \[Rule] {None, None}, Start \[Rule] {0,0}, Clamped \[Rule] False, and Raw \[Rule] False.";
Spherical::usage = "Spherical[width, height] Returns image of quadratic phase profile. Default options ScalingFactor \[Rule] 1, PhaseOffset \[Rule] 0, and Raw \[Rule] False.";
GetLGLink::usage = "GetLGLink[] returns the Link object.";
CloseLG::usage = "CloseLG[] close the Link if it exists."
StartLGDebug::usage = "StartLGDebug[name] connect to WSTP link with given link name (for debugging)";


Begin["`Private`"];


LGpath = None;
LGlink = None;
DEFAULTS = {
	Width -> 500,
	Height -> 500,
	WaveVector -> 1.0,
	BeamWaist -> 10.0,
	ScalingFactor -> 1.0,
	Plane -> "XY",
	Z -> 0.01,
	X -> 0,
	Y -> 0,
	Clamp-> True
};


StartLG[path_, binpath_:"bin_win32/Release/lg.exe"] := 
	Module[{},
		If[LGlink != None, Uninstall[LGlink]; LGlink = None;];
		LGlink = Install[FileNameJoin[{path, binpath}]];
		LGpath = FileNameJoin[{path, binpath}];
	];


StartLGDebug[name_] := 
	Module[{},
		If[LGlink =!= None, Uninstall[LGlink]; LGlink = None;];
		LGlink = Install[name, LinkMode -> Connect];
		LGpath = None;
	];


CloseLG[] :=
	Module[{},
		If[LGlink =!= None, 
			If[LGpath === None, LinkClose[LGlink], Uninstall[LGlink]];
			LGlink = None;
		];
	];


GetLGLink[] := LGlink;


Spherical[width_,height_,
	OptionsPattern[{
		ScalingFactor-> 1,
		PhaseOffset-> 0,
		Raw-> False
	}]] := 
	Module[
		{scale,window, data, offset,w, h, bytes},
		scale = N[OptionValue[ScalingFactor]];
		offset = N[OptionValue[PhaseOffset]];
		
		w = IntegerPart@width;
		h = IntegerPart@height;

		bytes = Check[Global`SphericalPhase[scale,offset,w,h], Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


BlazeImage[image_, angle_,
	OptionsPattern[{
		Frequency-> 0.5,
		PhaseOffset-> 0,
		Size->{None, None}, 
		Start-> {0,0},
		Clamped-> False,
		Raw-> False
	}]] := 
	Module[ 
		{scale,w,h,x,y,data,window,ang,offset, fc,cl},
		scale = N[OptionValue[Frequency]];
		data=If[ImageQ@image,ImageData[image], image];
		window=OptionValue[Size];
		w=IntegerPart@If[NumberQ@window[[1]],window[[1]], Dimensions[data][[2]]];
		h=IntegerPart@If[NumberQ@window[[2]],window[[2]], Dimensions[data][[1]]];
		x =IntegerPart[OptionValue[Start][[1]]];
		y=IntegerPart[OptionValue[Start][[2]]];
		offset=N[OptionValue[PhaseOffset]];
		ang =N[angle];
		cl= If[OptionValue[Clamped], 1, 0];

		bytes=Check[Global`Blaze[ang,scale,offset,x,y,w,h,cl,data], Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


LGMag[l_, p_, 
	OptionsPattern[Join[
		{Raw-> False},
		DEFAULTS
	]]] := 
	Module[
		{scale,w0,w,h,cl,z,bytes,k,li,pi,plane,transverse,usex,v},
		scale = N[OptionValue[ScalingFactor]];
		w=IntegerPart[OptionValue[Width]];
		h=IntegerPart[OptionValue[Height]];
		cl=If[OptionValue[Clamp],1, 0];
		li =IntegerPart[l];
		pi=IntegerPart[p];

		k=N[OptionValue[WaveVector]];
		w0=N[OptionValue[BeamWaist]];
		z=N[OptionValue[Z]];
		
		plane=OptionValue[Plane];
		usex = If[plane == "ZX", 1, 0];
		v = N@If[usex == 1, OptionValue[Y], OptionValue[X]];
		transverse= (usex == 1 || plane == "ZY");

		bytes = Check[
					If[!transverse,
						Global`LaguerreGuassianMag[w, h, scale, li, pi, k,w0,z,cl],
						Global`LaguerreGuassianMagZX[w, h, scale, li, pi, k, w0, v, usex,cl]
					],
					Abort[]
				];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


LGPhase[l_, p_, 
	OptionsPattern[Join[
		{Raw-> False},
		DEFAULTS
	]]]:= 
	Module[
		{scale,w0,w,h,cl,z,bytes,k,li,pi,plane, usex, transverse,v},
		scale = N[OptionValue[ScalingFactor]];
		w = IntegerPart[OptionValue[Width]];
		h = IntegerPart[OptionValue[Height]];
		cl = If[OptionValue[Clamp],1, 0];
		li = IntegerPart[l];
		pi = IntegerPart[p];

		k = N[OptionValue[WaveVector]];
		w0 = N[OptionValue[BeamWaist]];
		z = N[OptionValue[Z]];
		
		plane=OptionValue[Plane];
		usex = If[plane == "ZX", 1, 0];
		v = N@If[usex == 1, OptionValue[Y], OptionValue[X]];
		transverse= (usex == 1 || plane == "ZY");

		bytes = Check[
					If[!transverse,
						Global`LaguerreGuassianPhase[w, h, scale, li, pi, k,w0,z,cl],
						Global`LaguerreGuassianPhaseZX[w, h, scale, li, pi, k, w0, v,usex, cl]
					],
					Abort[]
				];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


TwoLGPhase[l0_, p0_,l1_,p1_, 
	OptionsPattern[Join[
		{Raw-> False},
		DEFAULTS
	]]]:= 
	Module[
		{scale,w0,w,h,cl,z,bytes,k,l1i,l0i,p1i,p0i,usex,transverse,v,plane},
		scale = N[OptionValue[ScalingFactor]];
		w=IntegerPart[OptionValue[Width]];
		h=IntegerPart[OptionValue[Height]];
		cl=If[OptionValue[Clamp],1, 0];
		l0i =IntegerPart[l0];
		p0i =IntegerPart[p0];
		l1i =IntegerPart[l1];
		p1i =IntegerPart[p1];

		k=N[OptionValue[WaveVector]];
		w0=N[OptionValue[BeamWaist]];
		z=N[OptionValue[Z]];
		
		plane=OptionValue[Plane];
		usex = If[plane == "ZX", 1, 0];
		v = N@If[usex == 1, OptionValue[Y], OptionValue[X]];
		transverse= (usex == 1 || plane == "ZY");

		bytes=Check[
			If[!transverse,
				Global`TwoLaguerreGuassianPhase[w, h, scale, l0i, p0i,l1i,p1i, k,w0,z,cl],
				Global`TwoLaguerreGuassianPhaseZX[w, h, scale, l0i, p0i,l1i,p1i, k,w0,v,usex,cl]
			],
			Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


TwoLGMag[l0_, p0_,l1_,p1_, 
	OptionsPattern[Join[
		{Raw-> False},
		DEFAULTS
	]]]:= 
	Module[
		{scale,w0,w,h,cl,z,bytes,k,l1i,l0i,p1i,p0i,usex,transverse,v,plane},
		
		scale = N[OptionValue[ScalingFactor]];
		w = IntegerPart[OptionValue[Width]];
		h = IntegerPart[OptionValue[Height]];
		cl = If[OptionValue[Clamp],1, 0];
		l0i = IntegerPart[l0];
		p0i = IntegerPart[p0];
		l1i = IntegerPart[l1];
		p1i = IntegerPart[p1];

		k = N[OptionValue[WaveVector]];
		w0 = N[OptionValue[BeamWaist]];
		z = N[OptionValue[Z]];
	
		plane=OptionValue[Plane];
		usex = If[plane == "ZX", 1, 0];
		v = N@If[usex == 1, OptionValue[Y], OptionValue[X]];
		transverse= (usex == 1 || plane == "ZY");

		bytes=Check[
			If[!transverse,
				Global`TwoLaguerreGuassianMag[w, h, scale, l0i, p0i,l1i,p1i, k,w0,z,cl],
				Global`TwoLaguerreGuassianMagZX[w, h, scale, l0i, p0i,l1i,p1i, k,w0,v,usex,cl]
			],
			Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


End[];


EndPackage[];
