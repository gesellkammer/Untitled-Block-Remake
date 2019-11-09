/*
PhasorLine: outputs sample indices to read a Buffer between
times x0 and x1

LinePhasor: like Line, only recurring and triggerable

TimedSchmitt: like Schmidt, but signal needs to stay a certian duration on a given state
*/


PhasorLine {
    *ar { arg x0=0, x1=1, dur=1, trig=0, mul=1, add=0;
		var impulse = Impulse.ar(1/dur);
		var trig2 = if( impulse > trig, impulse, trig );
        var sr = SampleRate.ir;
        var speed = (x1 - x0) / dur;
        ^Phasor.ar(trig2, speed/sr, x0*sr, x1*sr, x0*sr).madd(mul, add);
    }

	*kr {arg x0=0, x1=1, dur=1, trig=0, mul=1, add=0;
		var impulse = Impulse.kr(1/dur);
		var trig2 = if( impulse > trig, impulse, trig );
        var sr = ControlRate.ir;
        var speed = (x1 - x0) / dur;
        ^Phasor.kr(trig2, speed/sr, x0*sr, x1*sr, x0*sr).madd(mul, add);
    }
}


LinePhasor {
	*ar { arg start=0, end=1, dur=1, trig=0, resetPos=nil, mul=1, add=0;
		var sr = SampleRate.ir;
		var speed = (end - start) / dur;
		resetPos = resetPos ? start;
		^Phasor.ar(trig, speed/sr, start, end, resetPos).madd(mul, add);
	}

	*kr { arg start=0, end=1, dur=1, trig=0, resetPos=nil, mul=1, add=0;
		var sr = ControlRate.ir;
		var speed = (end - start) / dur;
		resetPos = resetPos ? start;
		^Phasor.kr(trig, speed/sr, start, end, resetPos).madd(mul, add);
	}
}


TrigSchmitt {
	/*
	A switch with a min amount of time open and closed

	in goes possitve      -> output is 1 for at least hitime seconds
	in goes non-possitive -> output is 0 for at least lotime seconds

	*/
	*ar { arg in, lotime=0.1, hitime=0.1;
		var up = in | Trig1.ar(in, hitime);
		var inv = 1 - in;
		var down = inv | Trig1.ar(inv, lotime);
		var state = Schmidt.ar(up + (down.neg), -0.999, 0.999);
		^state;
	}

	*kr { arg in, lotime=0.1, hitime=0.1;
		var up = in | Trig1.kr(in, hitime);
		var inv = 1 - in;
		var down = inv | Trig1.kr(inv, lotime);
		var state = Schmidt.kr(up + (down.neg), -0.999, 0.999);
		^state;
	}
}