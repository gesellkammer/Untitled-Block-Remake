/*

Remake of Untitled Block

Dependencies:

this patch has no dependencies

How the original worked:

* There are two audio streams
  * one very low (we will call this SUB) with
    frequencies between 10 and 200 Hz, in an almost gaussian distribution around 50 Hz
  * a wall of sound (we will call it WALL) with all other frequencies
* The SUB stream is sent to a compressor as its key stream (the stream controlling the amount
of compression), and a mix of the SUB and the WALL are sent to the compressor as the signal
being compressed. The compressor's ratio is set very low, so as to almost eliminate the sound
when the gate is open. The threshold of the compressor is tuned precisely to the level of the
key signal, with the result that it is constantly jumping up and down this threshold, opening and closing the compresser.
* The compressor is set to a vary high ratio with no gain compensation, creating in fact
an "inverted gate"
* Eventually the SUB is gained down, and without this stream present, the compressor is
closed and the WALL of sound is present constantly
* At certain moments, coordinated with the turning on of the SUB after beeing turned off,
a different eq curve is applied to the WALL, either accentuating certain regions (never in
a pronounced way which would result in the actual perception of pitch) or filtering out others.

Remake:

One of the main elements of this approach is the implementation of an amplitude follower
to the key signal, so that it is not necessary to adjust the threshold or the gain of the key
signal.
Like in the original there are two streams, sub and wall, passed in this case as
buffers (but this can be changed, allowing for control of the underlying audio)
The gate of the compressor has further control of its min. open and closed times, so as to
control the overall behaviour of the texture. This times are also subjected to a gaussian variation
centered on the given values, to avoid regularities.

*/

( // setup

fork {
	var patchDir = thisProcess.nowExecutingPath.dirname;

	[~sub, ~wall].do {|buf| buf.free; buf = nil };

	s.sync;

	~sub = Buffer.read(s, patchDir +/+ "snd/SUB3-44100.flac");
	~wall = Buffer.read(s, patchDir +/+ "snd/VENTWHITE-44100.flac");
}
)

(
var loopbuf2 = {|buf, marginLeft, marginRight|
	var bufsr = BufSampleRate.ir(buf);
	LoopBuf.ar(2, buf,
		startLoop:marginLeft*bufsr,
		endLoop:BufSamples.ir(buf)-(marginRight*bufsr));
};


var compress = {|sig, thresh= -24, knee=8, ratio=2, gain= -3, autogain=1, attack=0.0, release=0.1|
    /*
    -- the softknee begins at thresh - knee
    -- gain has only meaning when autogain is on
    -- a bigger ratio means more compression
    */
    var times = [attack, release];
    var slopeAbove = 1.0/ratio;
    var gain2 = Select.kr(autogain, [gain, (thresh.neg * ( 1 - slopeAbove )) + gain]); // autogain
    var amp = Amplitude.ar( sig, *times ).linlin(
        (thresh-knee).dbamp, (thresh+knee).dbamp, 0, 1
    ).clip(0,1).sqrt;
    knee = knee.max(0.0001);
    slopeAbove = amp.linlin(0, 1, slopeAbove, 1);
    thresh     = amp.linlin(0, 1, thresh, thresh-knee);
    Compander.ar( sig, sig, thresh.dbamp, 1, slopeAbove, *times ) * gain2.dbamp;
};


var guard = {|sig|
	Select.ar(CheckBadValues.ar(sig, 0, 0), [sig, DC.ar(0), DC.ar(0), sig]);
};


var limiter = {|sig, soft_thresh= -3, hard_thresh= -0.1, knee=3, ratio=4, attack=0.0008, release=0.1, lookahead=0|
    // a soft limiter
    var look = max(lookahead, ControlDur.ir);
    sig = compress.(sig, soft_thresh, knee, ratio, attack:attack, release:release, autogain:0);
    sig = Limiter.ar(sig, hard_thresh.dbamp, look);
    sig = guard.(sig);
    sig;
};

SynthDef(\blk, {
	|bufsub=0, bufwall=0, atk=0.005, rel=0.02, avgRate=10, keyEnvFactor=1,
	keyOffset=0, keyFactor=1, hiMinDur=0.02, loMinDur=0.05, gainlo=1, gainhi=0.05,
	subGain=2, wallGain=0.5, subKeyGain=1, wallKeyGain=0, durDev=0.15|
	/*
	bufsub: the buffer containing low freqs
	bufwall: the buffer containing the wall sound
	atk, rel: attack and release of the system
	avgRate: the rate at which the current level of the key signal is averaged. This
	is used to "follow" the key signal
	keyEnvFactor: the key signal is amp followed, and the up/down movements of this
	curve are derived as a factor of atk/rel. Use 1 to use the same attack/release times as
	the compressor
	keyOffset: a value to add to the key envelope to offset it artificialy from the threshold
	If set two high, the compressor remains always open. Normal settings: between -0.1 and 0
	hiMinDur: min. duration that the compressor gate will remain open (applying gainhi)
	loMinDur: min. duration that the compressor gate will remain closed (applying gainlo)
	gainhi: gain applied to the signal when the compressor is open
	gainlo: gain applied to the signal when the compressor is closed
	sugGain: gain of SUB passed to the output
	wallGain: gain of WALL passed to the output (signal = SUB*subGain + WALL*wallGain)
	subKeyGain: gain of SUB passed to the key input of the compressor
	wallKeyGain: gain of WALL passed to the key input of the compressor
	durDev: a deviation factor. A Gaussian distribution is calculated around hiMinDur and loMinDur,
	with this value as its deviation factor (0.1 is 10% deviation)

	*/
	var sub = loopbuf2.(bufsub, 1, 1);
	var wall = loopbuf2.(bufwall, 1, 1);
	var sig = (sub*subGain) + (wall*wallGain);
	var compkey = (sub[0] * subKeyGain) + (wall[0] * wallKeyGain);
	var keyenv = (
		Amplitude.ar(compkey, attackTime:atk*keyEnvFactor, releaseTime:rel*keyEnvFactor)
		* keyFactor
		+ keyOffset).clip(0, 1);
	var autothresh = (keyenv>0)*AverageOutput.ar(keyenv, Metro.kr(avgRate, 1)).lagud(1/avgRate, 0);
	var gainlosig = K2A.ar(gainlo);
	var gainhisig = K2A.ar(gainhi);
	var durChangeRate = 10;
	var loMinDurRnd = TGaussRand.kr(
		lo:loMinDur*(1-durDev),
		hi:loMinDur*(1+durDev),
		trig: Impulse.kr(durChangeRate, phase:0.5)).clip(0, 999);
	var hiMinDurRnd = TGaussRand.kr(
		lo:hiMinDur*(1-durDev),
		hi:hiMinDur*(1+durDev),
		trig: Impulse.kr(durChangeRate)).clip(0, 999);
	var compGate = (
		TrigSchmitt.ar(
			Schmidt.ar(keyenv, autothresh-keyOffset, autothresh+keyOffset),
			lotime:loMinDurRnd,
			hitime:hiMinDurRnd));
	var gain = Select.ar(compGate, [gainlosig, gainhisig]);
	sig = sig * gain.lagud(atk, rel);
	sig = CompanderD.ar(sig, thresh: -12.dbamp, slopeAbove:0.3, clampTime:atk, relaxTime:rel)*2;
	sig = limiter.(sig, attack:0.008, release:0.1);
	Out.ar(0, sig);
}).add;
)

x = Synth(\blk, args:[bufsub: ~sub, bufwall: ~wall])

x.set(\hiMinDur, 0.05, \loMinDur, 0.03, \keyGain, 2, \atk, 0.005, \rel, 0.02, \gainhi, 0.05, \gainlo, 1, \subGain, 1.5, \wallGain, 0.5, \keyFactor, 0.9, \avgRate, 8, \keyEnvFactor, 1, \durDev, 0.2)

