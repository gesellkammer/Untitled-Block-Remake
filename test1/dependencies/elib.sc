Elib {
	classvar <registry;
	classvar <>a4;
	*initClass {
		registry = IdentityDictionary();
		a4 = 442;
	}
	*moses { |seq, func|
		var trues = List();
		var falses = List();
		seq.do {|x|
			if( func.(x) ) {
				trues.add(x);
			} /*else*/ {
				falses.add(x);
			}
		};
		^[trues.array, falses.array]
	}
	*seq_order_lower { |a, b|
		/*
		returns which of both sequences should be ordered first
		*/
		var n;
		var max_n = min(a.size, b.size);
		max_n.do {|n|
			if( a[n] < b[n] ) {
				^true;
			}{
				if( a[n] > b[n] ) {
					^false;
				};
			}
		};
		^true
	}
	*seq_order_greater { |a, b|
		/*
		returns which of both sequences should be ordered first
		*/
		var n;
		var max_n = min(a.size, b.size);
		max_n.do {|n|
			if( a[n] > b[n] ) {
				^true;
			}{
				if( a[n] < b[n] ) {
					^false;
				};
			}
		};
		^true
	}
	*midi_init {
		MIDIClient.init;
		MIDIIn.connectAll;
	}

	*decay2fback {|decaydur, delaytime|
		// the feedback needed for a comb to decrease its volume in 60 dB
		// in the given duration for the indicated delaytime
		// ^exp(log(0.01) * delaytime / decaydur);
		^exp(-4.605170185988091 * delaytime / decaydur);
	}

	*call_regularly {|func, interval, clock=\app, autostop=true|
		var routine;
		if(clock.isKindOf(Clock).not) {
			clock = (app: AppClock, system: SystemClock).at(clock);
			clock = clock ? SystemClock;
		};
		// call at an interval and when pressed CMD + "."
		CmdPeriod.doOnce({
			func.();
		});
		routine = Routine {
			loop {
				func.();
				interval.wait;
			};
		}.play(clock);
		if( autostop ) {
			CmdPeriod.doOnce { routine.stop };
		};
		^routine;
	}

	*unpickle {|path|
		path = path.absolutePath;
		if(File.exists(path)) {
			var file = File(path, "r");
			var str = file.readAllString;
			var obj = str.interpret;
			^obj;
		}{
			^nil
		};
	}

	*pickle {|obj, path|
		var file = File(path.absolutePath, "w");
		var str = obj.asCompileString;
		file.write(str);
		file.close;
	}

	*allpairs {|seq|
		^all{:[seq[x], seq[y]], x <- (0..seq.size-2), y<-(x..seq.size-1), x!=y};
	}

	*ringmod {|...freqs|
		var pairs = freqs.allPairs(freqs);
		var sidebands = pairs.collect {|pair|
			var freq0 = pair[0];
			var freq1 = pair[1];
			[freq0+freq1, abs(freq0-freq1)];
		}.flatten.as(Set).as(Array);
		^sidebands;
	}

	*audioTestOutput{|firstchan=0, lastchan=1, durperchan=0.5, gain=0.5, kind=0|
		var numchans = lastchan - firstchan + 1;
		^{
			var chan = LinePhasor.kr(firstchan, lastchan + 1, numchans * durperchan).floor;
			var pink = PinkNoise.ar * 0.5;
			var pulses = Impulse.ar(10) * 0.5;
			var sine = SinOsc.ar(1000) * 0.5;
			SelectX.ar(kind, [pink, pulses]) * gain !> Out.ar(chan, _);
			chan.pollchanged(label: "channel");
		}.play;
	}

	*asmidi{|x|
		if(x.isKindOf(String)) {
			^x.n2m;
		} {
			^x.asFloat;
		};
	}

	*setA4 {|freq|
		// this is used by m2f and f2m
		Priv_PitchConvertion.setA4(freq);
	}

	*difftoneSource{|difftone, interval, resolution=0.5|
		var diffFreq = Elib.asmidi(difftone).m2f;
		var ratio = interval.midiratio;
		var f0 = diffFreq / (ratio - 1);
		var f1 = f0 * ratio;
		^[f0.f2m.round(resolution), f1.f2m.round(resolution)];
	}

	*difftoneSourcesInRange{|difftone, minnote, maxnote, resolution=0.5, maxdiff=0.5, mininterval=0, maxinterval=7|
		// find all sources within the given range
		var desired = Elib.asmidi(difftone);
		var interval0 = max(resolution, mininterval);
		var intervals = (interval0,interval0+resolution..maxinterval);
		var pairs = intervals.collect{|interval| Elib.difftoneSource(desired, interval, resolution)};
		var midi0 = Elib.asmidi(minnote);
		var midi1 = Elib.asmidi(maxnote);
		pairs.postln;
		^all{:pair, pair <- pairs,
			(midi0 <= pair[0]) && (pair[1] <= midi1) && ((Elib.difftone(*pair) - desired).abs <= maxdiff)};
	}

	*difftone{|midi1, midi2|
		// returns the difftone between midi1 and midi2 as midinote
		^(Elib.asmidi(midi1).m2f - Elib.asmidi(midi2).m2f).abs.f2m;
	}

	*pitchbend2cents{|pitchbend, maxcents=200|
		^(((pitchbend/16383.0)*(maxcents*2.0)) - maxcents + 0.5).floor;
	}

	*cents2pitchbend{|cents, maxcents=200|
		^((cents+maxcents)/(maxcents*2.0) * 16383.0 + 0.5).floor;
	}
}

Difftone {
	*pair{|difftone, interval, resolution=0.5|
		var diffFreq = Elib.asmidi(difftone).m2f;
		var ratio = interval.midiratio;
		var f0 = diffFreq / (ratio - 1);
		var f1 = f0 * ratio;
		^[f0.f2m.round(resolution), f1.f2m.round(resolution)];
	}

	*findPairs{|difftone, minnote, maxnote, resolution=0.5, maxdiff=0.5, mininterval=0, maxinterval=7|
		// find all sources within the given range
		var desired = Elib.asmidi(difftone);
		var interval0 = max(resolution, mininterval);
		var intervals = (interval0,interval0+resolution..maxinterval);
		var pairs = intervals.collect{|interval| Difftone.pair(desired, interval, resolution)};
		var midi0 = Elib.asmidi(minnote);
		var midi1 = Elib.asmidi(maxnote);
		pairs.postln;
		^all{:pair, pair <- pairs,
			(midi0 <= pair[0]) && (pair[1] <= midi1) && ((Difftone.result(*pair) - desired).abs <= maxdiff)};
	}

	*result{|midi1, midi2|
		// returns the difftone between midi1 and midi2 as midinote
		^(Elib.asmidi(midi1).m2f - Elib.asmidi(midi2).m2f).abs.f2m;
	}
}

MovingAvg {
	var <wsize, <>avg, partialratio, oneratio;
	*new {|winsize=5, startvalue=0|
		^super.new.init(winsize, startvalue)
	}
	init {|winsize=5, startvalue=0|
		avg = startvalue;
		wsize = winsize;
		partialratio = (winsize - 1) / winsize;
		oneratio = 1 / winsize;
	}
	value  {|x| ^this.update(x); }
	update {|x|
		avg = (avg * partialratio) + (x * oneratio);
		^avg;
	}
}

Insert : UGen {
	/* like an insert in a mixing desk:
	Example: get a sound from hardware channel 0, apply some filtering, then sendit to channel 4
	for some hardware processing (or maybe some plugin via Jack?) and get it back via channel 4

	inch defaults to outch

	var source = SoundIn.ar(0) !> LPF.ar(_, 4000) !> Insert(_, 4, 4);
	*/
	*ar {|sig, outch, inch|
		inch = inch ? outch;
		Out.ar(outch, sig);
		^SoundIn.ar(inch);
	}
	*kr {|sig, outch, inch|
		inch = inch ? outch;
		Out.kr(outch, sig);
		^In.kr(inch);

	}
	*new {|sig, outch, inch|
		var rate = sig.asArray.collect(_.rate).unbubble;
		if( rate == 'audio' ) {
			^Insert.ar(sig, outch, inch);
		} {
			^Insert.kr(sig, outch, inch);
		};
	}
}

FFTL {
	// A simple wrapper around FFT to allocate a local buffer and perform the fft in one step
	// { FFTL(SoundIn.ar(0), 2048) !> ... !> IFFT(_) !> Out([0, 1]) }.play
	*new {|in, n=2048, hop=0.5, wintype=0, active=1, winsize=0|
		var nchannels = 1;  // TODO: could this be determined by the input??
		^FFT(LocalBuf(n, nchannels), in, hop, wintype, active, winsize);
	}
}

PrBinaryAgent {
	/*
	see the operator >< and !><

	Examples:

	SinOsc.ar(440) >< Pulse.ar(440) @0.5
	WhiteNoise.ar(0.2) !>< LPF.ar(_, 1000) @0.8

	*/

	var a, b, blendfunc;
	*new {|a, b, func|
		if(a.isNil) {
			"Attempted to create a BinaryAgent with nil as first argument".throw;
		};
		if(b.isNil) {
			"Attempted to create a BinaryAgent with nil as second argument".throw;
		};
		^super.new.init(a, b, func);
	}
	init {|a0, b0, func|
		a = a0;
		b = b0;
		blendfunc = func;
	}
	at {|which|
		^blendfunc.(which, a, b);
	}
	@ {|which|  // a synonim of at
		^this.at(which);
	}
	< {|which|
		// usage: SinOsc.ar(440) >@which< Saw.ar(440)   // performs a blend between SinOsc and Saw depending on which (0-1)
		^blendfunc.(which, a, b);
	}
}

+ SequenceableCollection {
	moses {|func|
		// returns 2 arrays: [truearray, falsearray]
		var trues = List();
		var falses = List();
		this.do {|x|
			if( func.(x) ) {
				trues.add(x);
			} /*else*/ {
				falses.add(x);
			}
		};
		^[trues.array, falses.array]
	}
	greater {|otherseq|
		// return of this sequence should be ordered later than otherseq.
		// returns false if they are equal.
		^Elib.seq_order_greater(this, otherseq);
	}
	allPairs {
		// ^Elib.allpairs(this);
		^all{:[this[x], this[y]], x <- (0..this.size-2), y<-(x..this.size-1), x!=y};
	}

}

+ Env {
	*fromxy { |xs, ys, curve='lin'|
		var deltas = (xs.size - 1).collect {|i|
			xs[i+1] - xs[i]
		};
		^Env.new(ys, deltas, curve);
	}
	*frommatrix { |matrix, curve='lin'|
		var xs, ys;
		# xs, ys = matrix.clump(2).flop;
		^Env.fromxy(xs, ys, curve);
	}
	*frompairs { |... pairs|
		var xs, ys;
		# xs, ys = pairs.flop;
		^Env.fromxy(xs, ys);
	}
	*bpf{ |...args|
		// the first arg is the curve, the rest is a seq
		var xs, ys;
		# xs, ys = args[1..].clump(2).flop;
		^Env.fromxy(xs, ys, curve:args[0]);
	}
	*bpfexpon { |...seq|
		var xs, ys;
		# xs, ys = seq.clump(2).flop;
		^Env.fromxy(xs, ys, curve:\exp);
	}
	*bpfstep { |...seq|
		var xs, ys;
		# xs, ys = seq.clump(2).flop;
		^Env.fromxy(xs, ys, curve:\step);
	}
	*bpfwelch { |...seq|
		var xs, ys;
		# xs, ys = seq.clump(2).flop;
		^Env.fromxy(xs, ys, curve:\welch);
	}
	*bpflinear { |...seq|
		var xs, ys;
		# xs, ys = seq.clump(2).flop;
		^Env.fromxy(xs, ys);
	}
	*bpfsine { |...seq|
		var xs, ys;
		# xs, ys = seq.clump(2).flop;
		^Env.fromxy(xs, ys, curve:\sine);
	}
	*bpfsquared { |...seq|
		var xs, ys;
		# xs, ys = seq.clump(2).flop;
		^Env.fromxy(xs, ys, curve:\squared);
	}

	*fade {|dur, fade=0.1, level=1, curve='lin'|
		^Env.bpf(curve, 0, 0, fade, level, dur-fade, level, dur, 0);
	}

	*afade {|dur, fade=0.1, level=1, curve='lin', done=2|
		^Env.fade(dur:dur, fade:fade, level:level, curve:curve).ar(doneAction:done);
	}

	asBuffer {|numpoints=8192, action, server|
		/*
		Env.bpfsine(0, 0, 1, 1).asBuffer(action:{|buf| ~mybuf = buf })
		*/
		assert{ numpoints.notNil && numpoints > 0 };
		server = server ? Server.default;
		^Buffer.sendCollection(server, this.discretize(numpoints), action:action)
	}
}


+ String {
	path_split {
		var a = this.split($/);
		var basename = a.last;
		var path = a[.. a.size-2].join($/);
		^[path, basename];
	}
	oscfunc { |func, port|
		^OSCFunc(func:func, path:this, recvPort:port);
	}
	n2m {
		^Priv_PitchConvertion.note_to_midi(this);
	}

	n2f {
		^this.n2m.m2f;
	}
}

+ Symbol {
	oscfunc { |func, port|
		^OSCFunc(func:func, path:this, recvPort:port);
	}
	n2m {
		^Priv_PitchConvertion.note_to_midi(this.asString);
	}

	n2f {
		^this.n2m.m2f;
	}
}

+ Buffer {
	*cue {|server, path|
		var channels = SoundFile(path).numChannels;
		^this.cueSoundFile(server, path, 0, channels);
	}
}

+ Synth {
	tail {|target|
		fork {
			this.server.sync;
			this.moveToTail(target)
		}
	}
	after {|target|
		fork {
			this.server.sync;
			this.moveAfter(target);
		}
	}
	before {|target|
		fork {
			this.server.sync;
			this.moveBefore(target);
		}
	}
	head {|target|
		fork {
			this.server.sync;
			this.moveToHead(target);
		}
	}
	controlNames {
		/*
		Example
		-------

		a = {|freq=440|
		SinOsc.ar(freq) * \amp.kr(0.2) !> Out(0)
		}.play;

		a.controlNames.postln;   // --> [ i_out, freq, amp ]

		*/
		^this.synthDef.allControlNames.collect(_.name);
	}


}

+ Array {
	nearest {|x|
		var x0, x1;
		var i = this.indexOfGreaterThan(x);
		if( i.isNil ) {
			^this.last
		};
		if( i == 0 ) {
			^this.first
		};
		x1 = this.at(i);
		x0 = this.at(i - 1);
		if( (x1 - x).abs < (x0 - x).abs ) {
			^x1
		} /* else */ {
			^x0
		};
	}
	zip {|other|
		^[this, other].flop;
	}
	>< {|other|
		if( other.isKindOf(SequenceableCollection) ) {
			{ other.size == this.size }.assert("><: only sequences of equal size can be blended");

			^PrBinaryAgent(this, other, func: {|by, a, b|
				a.size.collect {|i|
					var ai = a[i];
					var bi = b[i];
					ai.blend(bi, by);
				};
			});
		} {
			^PrBinaryAgent(this, other, func: {|by, as, b|
				as.collect {|a|
					a.blend(b, by);
				};
			});
		};
	}
	>@ {|which|
		^PrBinaryAgent(this, which, func:{|b, a, by| a.blend(b, by)});
	}
}

+ True {
	* {|other|
		^other;
	}
}

+ False {
	* {|other|
		^0;
	}
}


+ MIDIOut {
	mmcPlay {
		this.sysex(Int8Array(6).addAll([16rF0, 16r7F, 0, 6, 2, 16rF7]));
	}
	mmcStop {
		this.sysex(Int8Array(6).addAll([16rF0, 16r7F, 0, 6, 1, 16rF7]));
	}
	mmcRec {
		this.sysex(Int8Array(6).addAll([16rF0, 16r7F, 0, 6, 6, 16rF7]));
	}
	mmcRecStop {
		this.sysex(Int8Array(6).addAll([16rF0, 16r7F, 0, 7, 2, 16rF7]));
	}
	mmcPause {
		this.sysex(Int8Array(6).addAll([16rF0, 16r7F, 0, 6, 9, 16rF7]));
	}
	mmcFastForward {
		this.sysex(Int8Array(6).addAll([16rF0, 16r7F, 0, 6, 4, 16rF7]));
	}
	mmcRewind {
		this.sysex(Int8Array(6).addAll([16rF0, 16r7F, 0, 6, 5, 16rF7]));
	}
	*newMatchingRegex {|regex|
		var destinationIndex;
		if( MIDIClient.destinations.isNil ) {
			MIDIClient.init;
		};
		destinationIndex = MIDIClient.destinations.detectIndex {|dest|
			regex.matchRegexp(dest.device)
		};
		if( destinationIndex.notNil ) {
			^MIDIOut(destinationIndex, MIDIClient.destinations[destinationIndex].uid);
		}
	}
	*newMatching {|partOfName|
		var found = MIDIClient.destinations.select { |dest| dest.device.contains(partOfName) };
		var dev;
		if( found.size > 0 ) {
			dev = found[0];
			^MIDIOut.newByName(dev.device, dev.name);
		} {
			^nil;
		};
	}

}

+ MIDIClient {
	*getSourcesMatching { |device_pattern, port_pattern|
		/*
		- Returns an array of sources matching the pattern as a regex
		- Will return an empty array if no sources match

		Arguments:
		-- device_pattern: a regexp pattern
		-- port_pattern: a regexp pattern, or nil to match everything
		*/
		var sources;
		if( device_pattern.isNil ) {
			throw("MIDIClient.getSourcesMatching: device_pattern should be a string, not nil");
		};
		sources = MIDIClient.sources.select { |source| device_pattern.matchRegexp(source.device) };
		if ( (sources.size > 0) and: (port_pattern.notNil) ) {
			sources = sources.select {|source| port_pattern.matchRegexp(source.port) };
		};
		^sources;
	}

	*getDestinationsMatching {|device_pattern, port_pattern|
		var destinations = MIDIClient.destinations.select { |dest|
			device_pattern.matchRegexp(dest.device) };
		if ((destinations.size > 0) and: (port_pattern.notNil)) {
			destinations = destinations.select{|dest| port_pattern.matchRegexp(destinations.port)};
		};
		^destinations;
	}

	*getSource { |device, port|
		/* Get the midi source whose device and port match the given ones
		device and port can be a REGEX
		nil means ".*"
		if not found, returns nil
		if found, return always the first source found, even when there are more.
		*/
		var found;
		if( device.isNil ) {
			throw("MIDIClient.getSource: device should be a string, not nil");
		};
		found = MIDIClient.getSourcesMatching(device, port);
		^ if( found.size == 0, nil, {found[0]} );
	}

	*uidToName { |uid|
		/* convert a uid representing a midisource to its name */
		var found = MIDIClient.sources.select { |source| source.uid == uid };
		if( found.size > 0 ) {
			^found[0];
		}{
			^nil
		};
	}

	*monitor {
		var mon = ();
		if(MIDIClient.initialized.not) {
			MIDIClient.init;
			MIDIIn.connectAll;
		};
		mon.cc = MIDIFunc.cc({|val, cc, ch, uid| ["CC", MIDIClient.uidToName(uid), ch, cc, val].postln});
		mon.noteon = MIDIFunc.noteOn({|vel, note, ch, uid| ["NOTEON", MIDIClient.uidToName(uid), ch, note, vel].postln});
		^mon;
	}
}

+ Synth {
	synthDef {
		^SynthDefStorage.synthDefDict[defName][0]
	}
	getArgIndex { |argname|
		var synthdef = this.synthDef();
		var argname_assymbol = argname.asSymbol;
		synthdef.allControlNames.do {|control, control_index|
			if( control.name == argname_assymbol) {
				^control_index;
			};
		}
		^nil;
	}
}

+ SynthDef {
	getSpecs {
		^SynthDescLib.global[name.asSymbol].tryPerform(\metadata).tryPerform(\at, \specs);
	}
}

+ Integer {
	linspace {|start, end|
		var dx = (end - start) / (this - 1);
		^(start, start+dx .. end);
	}
}

+ Changed {
	*new { arg input, threshold = 0;
		var rate = input.asArray.collect(_.rate).unbubble;
		if( rate == 'audio' ) {
			^Changed.ar(input, threshold);
		} {
			^Changed.kr(input, threshold);
		};
	}
}

+ SoundIn {
	*new { |chan, mul=1, add=0|
		^SoundIn.ar(chan, mul, add);
	}
}

+ UGen {
	pollchanged {|label, trigid= -1|
		^Poll(Changed(this), this, label, trigid);
	}
	xfade {|other|
		/*
		USage:
		{ |which|
		SinOsc.ar(freq) xfade: Saw.ar(freq) @which   // which between 0-1
		}.play
		*/
		^PrBinaryAgent(this, other, func:{|by, a, b| SelectX2(by, a, b)});
	}
	>< {|other|
		// Usage:
		// { |which|
		//    SinOsc.ar(freq) >< Saw.ar(freq) at: which  // or SinOsc.ar(freq) >< Saw.ar(freq) @which
		// }.play
		//
		if(other.isNil) {
			"%s Cannot crossfade with nil".format(this.class.name).throw;
		};
		^PrBinaryAgent(this, other, func:{|by, a, b| a.blend(b, by)});
	}
	>@ {|which|
		// use it like this: SinOsc.ar(440) >@0.2< SinOsc.ar(880)
		// this is the same as SinOsc.ar(440) >< SinOsc.ar(880) at:0.2
		^PrBinaryAgent(this, which, func:{|b, a, by| a.blend(b, by)});
	}

	mul {|factor|
		// use case is to avoid parenthesis over big chunks:
		// SinOsc.ar(440) + (SinOsc.ar(880) * 0.2)   --> SinOsc.ar(440) + SinOsc.ar(880).mul(0.2)
		^(this * factor);
	}
}

+ PV_ChainUGen {
	copy {
		^PV_Copy(this, LocalBuf(this.numFrames));
	}
	copyphase {
		^PV_CopyPhase(this, LocalBuf(this.numFrames));
	}
	numFrames {
		var numframes = this.widthFirstAntecedents.collect {|ante|
			if( ante.isKindOf(LocalBuf) ) {
				ante.numFrames;
			}{
				0;
			};
		};
		numframes = numframes.select(_ > 0);
		if( numframes.size >= 1 ) {
			^numframes[0];
		}{
			Error("numframes: no antecedent reporting # of frames of this ChainUGen!").throw;
		}
	}
	numChannels {
		var numchannels = this.widthFirstAntecedents.collect {|ante|
			if( ante.isKindOf(LocalBuf) ) {
				ante.numChannels;
			}{
				0;
			};
		};
		numchannels = numchannels.select(_ > 0);
		if( numchannels.size >= 1 ) {
			^numchannels[0];
		}{
			Error("numchannels: no antecedent reporting # of channels of this ChainUGen!").throw;
		}
	}
}


+ SimpleNumber {
	m2n {
		^Priv_PitchConvertion.midi_to_note(this);
	}
	m2f {
		^(this.midicps * (Priv_PitchConvertion.a4/440));
	}
	f2m {
		^((this * (440/Priv_PitchConvertion.a4)).cpsmidi);
	}
	eps {
		^(this + (2e-16));
	}
	negeps {
		^(this - (2e-16));
	}

	r2i {
		^this.ratiomidi;
	}

	i2r {
		^this.midiratio;
	}

}

+ AbstractFunction {
	m2f {
		^(this.midicps * (Priv_PitchConvertion.a4/440));
	}
	f2m {
		^((this * (440/Priv_PitchConvertion.a4)).cpsmidi);
	}
}



+ Function {
	call_regularly {|interval, clock=\app, autostop=true|
		^Elib.call_regularly(this, interval, clock, autostop);
	}
}

+ Routine {
	autostop {
		"deprecated! Use untilCmdPeriod".warn;
		CmdPeriod.doOnce {this.stop};
		^this;
	}
	untilCmdPeriod {
		CmdPeriod.doOnce{this.stop};
		^this;
	}
}


+ OSCFunc {
	untilCmdPeriod {
		CmdPeriod.doOnce(this);
		^this;
	}
}
