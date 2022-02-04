(* ::Package:: *)

BeginPackage["LG`"];


LGMag::usage = "LGMag[l, p] Returns image of magnitude of LG mode l, p. Default options Size \[Rule] {256, 256}, WaveVector \[Rule] 2\[Pi], BeamWaist \[Rule] \!\(\*SuperscriptBox[\(10\), \(3\)]\), Delta \[Rule] 10.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False.";
LGPhase::usage = "LGPhase[l, p] Returns image of phase of LG mode l, p. Default options Size \[Rule] {256, 256}, WaveVector \[Rule] 2\[Pi], BeamWaist \[Rule] \!\(\*SuperscriptBox[\(10\), \(3\)]\), Delta \[Rule] 10.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False.";
TwoLGMag::uage = "TwoLGMag[l0, p0, l1, p1] Returns image of magnitude of sum of two LG modes. Default options  Size \[Rule] {256, 256}, WaveVector \[Rule] 2\[Pi], BeamWaist \[Rule] \!\(\*SuperscriptBox[\(10\), \(3\)]\), Delta \[Rule] 10.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False."  ;
TwoLGPhase::usage = "TwoLGPhase[l0, p0, l1, p1] Returns image of phase of sum of two LG modes. Default options  Size \[Rule] {256, 256}, WaveVector \[Rule] 2\[Pi], BeamWaist \[Rule] \!\(\*SuperscriptBox[\(10\), \(3\)]\), Delta \[Rule] 10.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False.";
StartLG::usage = "StartLG[path, binpath] Install WSTP executable. Must be run before anything else. 'path' is the absolute path of repository. path/binpath should resolve to path of the WSTP binary.";
BlazeImage::usage = "BlazeImage[image, angle] Returns 'image' with blazing added on top. 
Default options Frequency \[Rule] 0.5, PhaseOffset \[Rule] 0, Size \[Rule] {None, None}, Start \[Rule] {0,0}, Clamped \[Rule] False, and Raw \[Rule] False.";
Spherical::usage = "Spherical[width, height] Returns image of quadratic phase profile. Default options ScalingFactor \[Rule] 1, PhaseOffset \[Rule] 0, and Raw \[Rule] False.";
GetLGLink::usage = "GetLGLink[] returns the Link object.";
CloseLG::usage = "CloseLG[] close the Link if it exists."
StartLGDebug::usage = "StartLGDebug[name] connect to WSTP link with given link name (for debugging)";
PropagateField::usage = "PropagateField[field, dz] returns field propagated downstream by dz. Default options WaveVector \[Rule] 2\[Pi], Delta \[Rule] 10.0, FourierSpace \[Rule] False, AmpRange \[Rule] Null, AmpColorMap \[Rule] ThermometerColors, Raw \[Rule] False";


Begin["`Private`"];


LGpath = None;
LGlink = None;
DEFAULTS = {
	Size -> {256, 256},
	WaveVector -> 2\[Pi],
	BeamWaist -> 10^3,
	Delta -> 10.0,
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
		Delta-> 1,
		PhaseOffset-> 0,
		Raw-> False
	}]] := 
	Module[
		{scale,window, data, offset,w, h, bytes},
		scale = N[OptionValue[Delta]];
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
		{deltaw, deltah,w0,w,h,cl,z,bytes,k,li,pi,plane,transverse,usex,v},
		If[ListQ@OptionValue[Delta],
			{deltaw, deltah} = N/@OptionValue[Delta],
			deltaw = N[OptionValue[Delta]];
			deltah = deltaw;
		];
		If[ListQ@OptionValue[Size], 
			{w, h} = IntegerPart/@OptionValue[Size],
			w = IntegerPart[OptionValue[Size]]; 
			h = w;
		];
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
						Global`LaguerreGuassianMag[w, h, deltaw, deltah, li, pi, k,w0,z,cl],
						Global`LaguerreGuassianMagZX[w, h, deltaw, deltah, li, pi, k, w0, v, usex,cl]
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
		{deltaw,deltah,w0,w,h,cl,z,bytes,k,li,pi,plane, usex, transverse,v},
		If[ListQ@OptionValue[Delta],
			{deltaw, deltah} = N/@OptionValue[Delta],
			deltaw = N[OptionValue[Delta]];
			deltah = deltaw;
		];
		If[ListQ@OptionValue[Size], 
			{w, h} = IntegerPart/@OptionValue[Size],
			w = IntegerPart[OptionValue[Size]]; 
			h = w;
		];
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
						Global`LaguerreGuassianPhase[w, h, deltaw, deltah, li, pi, k,w0,z,cl],
						Global`LaguerreGuassianPhaseZX[w, h, deltaw, deltah, li, pi, k, w0, v,usex, cl]
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
		{deltaw, deltah,w0,w,h,cl,z,bytes,k,l1i,l0i,p1i,p0i,usex,transverse,v,plane},
		If[ListQ@OptionValue[Delta],
			{deltaw, deltah} = N/@OptionValue[Delta],
			deltaw = N[OptionValue[Delta]];
			deltah = deltaw;
		];
		If[ListQ@OptionValue[Size], 
			{w, h} = IntegerPart/@OptionValue[Size],
			w = IntegerPart[OptionValue[Size]]; 
			h = w;
		];
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
				Global`TwoLaguerreGuassianPhase[w, h, deltaw, deltah, l0i, p0i,l1i,p1i, k,w0,z,cl],
				Global`TwoLaguerreGuassianPhaseZX[w, h, deltaw, deltah, l0i, p0i,l1i,p1i, k,w0,v,usex,cl]
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
		{deltaw, deltah, w0,w,h,cl,z,bytes,k,l1i,l0i,p1i,p0i,usex,transverse,v,plane},
		If[ListQ@OptionValue[Delta],
			{deltaw, deltah} = N/@OptionValue[Delta],
			deltaw = N[OptionValue[Delta]];
			deltah = deltaw;
		];
		If[ListQ@OptionValue[Size], 
			{w, h} = IntegerPart/@OptionValue[Size],
			w = IntegerPart[OptionValue[Size]]; 
			h = w;
		];
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
				Global`TwoLaguerreGuassianMag[w, h, deltaw, deltah, l0i, p0i,l1i,p1i, k,w0,z,cl],
				Global`TwoLaguerreGuassianMagZX[w, h, deltaw, deltah, l0i, p0i,l1i,p1i, k,w0,v,usex,cl]
			],
			Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


PropagateField[field_, dz_, 
	OptionsPattern[Join[
		{Raw-> False, FourierSpace -> False, AmpRange -> Null, 
		AmpColorMap -> "ThermometerColors", TitleSuffix -> "", AmpColorFunctionScaling -> False, AmpImageSize -> Automatic, PhaseImageSize -> Automatic},
		DEFAULTS
	]]]:= 
	Module[
		{delta,dim,ft,z,k,data,p1,p2,amp,phase,arange,amax},
		
		delta = N[OptionValue[Delta]];
		dim = Dimensions[field];
		If[Length@dim != 3 || dim[[1]] != 2 || dim[[2]] != dim[[3]], 
			Message[Propagate::error, "field must be 2 x N x N"]; Abort[]];
		
		ft = If[OptionValue[FourierSpace],1, 0];
		z = N[dz];

		k = N[OptionValue[WaveVector]];
		{amp, phase} = Global`Propagate[k, delta, z, ft, field];
		
		If[OptionValue[Raw], 
			{amp,phase},
			arange = If[ListQ@OptionValue[AmpRange], OptionValue[AmpRange],
				amax = Max[amp];
				If[amax > 1, {0, amax}, {0, 1}]
			];
			p1 = ArrayPlot[amp, PlotLegends -> Automatic, PlotLabel -> "Amplitude"<>OptionValue[TitleSuffix], 
				ColorFunction -> OptionValue[AmpColorMap], PlotRange -> arange, ColorFunctionScaling -> OptionValue[AmpColorFunctionScaling],
				ImageSize -> OptionValue[AmpImageSize]];
			p2 = ArrayPlot[phase, PlotLegends -> BarLegend[{Hue, {0,1}}], PlotLabel -> "Phase"<>OptionValue[TitleSuffix], 
				ColorFunction -> Hue, PlotRange -> {0, 1}, ColorFunctionScaling -> False,
				ImageSize -> OptionValue[PhaseImageSize]];
			{p1, p2}
		]
	];


End[];


EndPackage[];
