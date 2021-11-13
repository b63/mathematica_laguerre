(* ::Package:: *)

BeginPackage["LG`"];


LGMag::usage = "LGMag[l, p] Returns image of magnitude of LG mode l, p. Default options Width \[Rule] 500, Height \[Rule] 500, WaveVector \[Rule] 1.0, BeamWaist \[Rule] 10.0, ScalingFactor \[Rule] 1.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False.";
LGPhase::usage = "LGPhase[l, p] Returns image of phase of LG mode l, p. Default options Width \[Rule] 500, Height \[Rule] 500, WaveVector \[Rule] 1.0, BeamWaist \[Rule] 10.0, ScalingFactor \[Rule] 1.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False.";
TwoLGMag::uage = "TwoLGMag[l0, p0, l1, p1] Returns image of magnitude of sum of two LG modes. Default options  Width \[Rule] 500, Height \[Rule] 500, WaveVector \[Rule] 1.0, BeamWaist \[Rule] 10.0, ScalingFactor \[Rule] 1.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False."  ;
TwoLGPhase::usage = "TwoLGPhase[l0, p0, l1, p1] Returns image of phase of sum of two LG modes. Default options  Width \[Rule] 500, Height \[Rule] 500, WaveVector \[Rule] 1.0, BeamWaist \[Rule] 10.0, ScalingFactor \[Rule] 1.0, Z \[Rule] 0.01, Clamp \[Rule] True, and Raw \[Rule] False.";
StartLG::usage = "StartLG[path, binpath] Install WSTP executable. Must be run before anything else. 'path' is the absolute path of repository. path/binpath should resolve to path of the WSTP binary.";
BlazeImage::usage = "BlazeImage[image, angle] Returns 'image' with blazing added on top. Default options AmplitudeOffset \[Rule] 0, Frequency \[Rule] 0, PhaseOffset \[Rule] 0, Size \[Rule] {None, None}, Start \[Rule] {0,0}, Clamped \[Rule] False, and Raw \[Rule] False.";
Spherical::usage = "Spherical[width, height] Returns image of quadratic phase profile. Default options ScalingFactor \[Rule] 1, PhaseOffset \[Rule] 0, and Raw \[Rule] False.";
GetLGLink::usage = "GetLGLink[] returns the Link object.";


Begin["`Private`"];


LGpath = None;
LGlink = None;
DEFAULTS = {
	Width-> 500,
	Height-> 500,
	WaveVector-> 1.0,
	BeamWaist-> 10.0,
	ScalingFactor-> 1.0,
	Z -> 0.01,
	Clamp-> True
};


StartLG[path_, binpath_:"bin_win32/Release/lg.exe"] := 
	Module[{},
		If[LGlink != None, Uninstall[LGlink]; LGlink = None;];
		LGlink = Install[FileNameJoin[{path, binpath}]];
		LGpath = FileNameJoin[{path, binpath}];
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

		bytes = Check[SphericalPhase[scale,offset,w,h], Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


BlazeImage[image_, angle_,
	OptionsPattern[{
		Amplitude-> 1, 
		AmplitudeOffset-> 0, 
		Frequency-> 1,
		PhaseOffset-> 0,
		Size->{None, None}, 
		Start-> {0,0},
		Clamped-> False,
		Raw-> False
	}]] := 
	Module[ 
		{scale,w,h,x,y,data,window,ang,offset, amp, ampOffset,fc,cl},
		scale = N[OptionValue[Frequency]];
		data=If[ImageQ@image,ImageData[image], image];
		window=OptionValue[Size];
		w=IntegerPart@If[NumberQ@window[[1]],window[[1]], Dimensions[data][[2]]];
		h=IntegerPart@If[NumberQ@window[[2]],window[[2]], Dimensions[data][[1]]];
		amp=N[OptionValue[Amplitude]];
		ampOffset=N[OptionValue[AmplitudeOffset]];
		x =IntegerPart[OptionValue[Start][[1]]];
		y=IntegerPart[OptionValue[Start][[2]]];
		offset=N[OptionValue[PhaseOffset]];
		ang =N[angle];
		cl= If[OptionValue[Clamped], 1, 0];

		bytes=Check[Blaze[ang,scale,offset,amp, ampOffset,x,y,w,h,cl,data], Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


LGMag[l_, p_, 
	OptionsPattern[Join[
		{Raw-> False},
		DEFAULTS
	]]] := 
	Module[
		{scale,w0,w,h,cl,z,bytes,k,li,pi},
		scale = N[OptionValue[ScalingFactor]];
		w=IntegerPart[OptionValue[Width]];
		h=IntegerPart[OptionValue[Height]];
		cl=If[OptionValue[Clamp],1, 0];
		li =IntegerPart[l];
		pi=IntegerPart[p];

		k=N[OptionValue[WaveVector]];
		w0=N[OptionValue[BeamWaist]];
		z=Re[OptionValue[Z]];

		bytes = Check[
					LaguerreGuassianMag[w, h, scale, li, pi, k,w0,z,cl],
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
		{scale,w0,w,h,cl,z,bytes,k,li,pi},
		scale = N[OptionValue[ScalingFactor]];
		w=IntegerPart[OptionValue[Width]];
		h=IntegerPart[OptionValue[Height]];
		cl=If[OptionValue[Clamp],1, 0];
		li =IntegerPart[l];
		pi  =IntegerPart[p];

		k=N[OptionValue[WaveVector]];
		w0=N[OptionValue[BeamWaist]];
		z=N[OptionValue[Z]];

		bytes=Check[LaguerreGuassianPhase[w, h,scale, li, pi, k,w0,z,cl], Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


TwoLGPhase[l0_, p0_,l1_,p1_, 
	OptionsPattern[Join[
		{Raw-> False},
		DEFAULTS
	]]]:= 
	Module[
		{scale,w0,w,h,cl,z,bytes,k,l1i,l0i,p1i,p0i},
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

		bytes=Check[TwoLaguerreGuassianPhase[
		w, h, scale, l0i, p0i,l1i,p1i, k,w0,z,cl
		],Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


TwoLGMag[l0_, p0_,l1_,p1_, 
	OptionsPattern[Join[
		{Raw-> False},
		DEFAULTS
	]]]:= 
	Module[
		{scale,w0,w,h,cl,z,bytes,k,l1i,l0i,p1i,p0i},
		
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

		bytes = Check[TwoLaguerreGuassianMag[w, h, scale,l0i, p0i,l1i,p1i, k,w0,z,cl], Abort[]];

		If[OptionValue[Raw], bytes, Image[bytes]]
	];


End[];


EndPackage[];