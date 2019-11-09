EM_Compress {
    *ar {|sig, thresh= -24, knee=8, ratio=2, gain= -3, autogain=1, attack=0.0, release=0.1|
        /*
        -- the softknee begins at thresh - knee
        -- gain has only meaning when autogain is on
        -- a bigger ratio means more compression compression
        */
        var times = [attack, release];
        var slopeAbove = 1.0/ratio;
        var gain2 = Select.kr(autogain, [gain, (thresh.neg * ( 1 - slopeAbove )) + gain]); // autogain
        var amp = Amplitude.ar( sig, *times ).linlin(
            (thresh-knee).dbamp, (thresh+knee).dbamp, 0, 1
        ).clip(0,1).sqrt;
        knee = knee.max(0.0001);
        // slopeAbove = 1.blend( slopeAbove, amp );
        // thresh = (thresh-knee).blend( thresh, amp );
        slopeAbove = amp.linlin(0, 1, slopeAbove, 1);
        thresh     = amp.linlin(0, 1, thresh, thresh-knee);
        ^Compander.ar( sig, sig, thresh.dbamp, 1, slopeAbove, *times ) * gain2.dbamp;
    }
    *new {|sig, thresh= -24, knee=8, ratio=2, gain= -3, autogain=1, attack=0.0, release=0.1|
        ^EM_Compress.ar(sig:sig, thresh:thresh, knee:knee, ratio:ratio, gain:gain, autogain:autogain, attack:attack, release:release);
    }
}

EM_Guard {
    *ar { |sig|
        ^Select.ar(CheckBadValues.ar(sig, 0, 0), [sig, DC.ar(0), DC.ar(0), sig]);
    }
    *new {|sig|
        ^EM_Guard.ar(sig);
    }
}

EM_Limit {
    // a soft limiter
    *ar {|sig, soft_thresh= -3, hard_thresh= -0.1, knee=3, ratio=4, attack=0.0008, release=0.5, lookahead=0|
        lookahead = max(lookahead, ControlDur.ir);
        ^(
            EM_Compress.ar(sig, soft_thresh, knee, ratio, attack:attack, release:release, autogain:0)
            !> Limiter.ar(_, hard_thresh.dbamp, lookahead)
            !> EM_Guard.ar(_)
        );
    }
    *new {|sig, soft_thresh= -3, hard_thresh= -0.1, knee=3, ratio=4, attack=0.0008, release=0.5, lookahead=0|
        ^EM_Limit.ar(sig, soft_thresh, hard_thresh, knee, ratio, attack, release, lookahead)
    }
}

EM_Denoise {
    *ar {|in, trig, smoothing=0.95, attendb=12, wet=0.7, i_fftsize=1024|
        /*
        Simple noise reduction algorithm.

        The signal `in` is monitored constantly. When trig is one, it is frozen and used as a noise sample to reduce from itself
        trig      : 0 to monitor for noise, 1 to freeze it (trigger when only noise is going through)
        smoothing : smooth the noise sample (0-1), as passed to PV_MagSmooth
        attendb   : a possitive db. How much attenuation should be applied
        wet       : mix the denoised signal with the original to mask artifacts.
                    A delay is applied to the original to compensate the delay of the FFT
        i_fftsize : the size of the fft. Bigger FFTs are more accurate but incur
                    in time smearing and most importantly in latency.
        */
        var source = in !> LeakDC.ar(_, 0.9);
        var orig  = FFTL(source, i_fftsize, hop:0.25);
        var smooth = PV_MagSmooth(orig.copy, smoothing);
        var smoothsaved = PV_Freeze(smooth, freeze:trig) !> PV_MagMulAdd(_, mul:attendb.dbamp);
        var reduced  = PV_MagSubtract(orig, smoothsaved, 1) !> IFFT(_);
        var fftdelay = i_fftsize/SampleRate.ir;
        var out = DelayN.ar(source, fftdelay, fftdelay) >< reduced at:wet;
        ^out;
    }
    *new {|in, trig, smoothing=0.95, attendb=12, wet=0.7, i_fftsize=1024|
        ^EM_Denoise.ar(in, trig, smoothing, attendb, wet, i_fftsize);
    }
    *split {|in, trig, smoothing=0.95, attendb=12, wet=0.7, i_fftsize=1024|
        /*
        Output three signals: the denoised signal, the noisy part which was removed, and the noise profile
        */
        var source = in !> LeakDC.ar(_, 0.9);
        var orig  = FFTL(source, i_fftsize, hop:0.25);
        var smooth = PV_MagSmooth(orig.copy, smoothing);
        var smoothsaved = PV_MagFreeze(smooth, freeze:trig) !> PV_MagMulAdd(_, mul:attendb.dbamp);
        var reduced  = PV_MagSubtract(orig.copy, smoothsaved, 1) !> IFFT(_);
        var fftdelay = i_fftsize/SampleRate.ir;
        var out = DelayN.ar(source, fftdelay, fftdelay) >< reduced at:wet;
        var noise = PV_CommonMag(orig, smoothsaved) !> IFFT(_);
        var noiseprofile = IFFT(smoothsaved);
        ^[out, noise, noiseprofile];
    }
}

EM_RemoveResidual {
    *ar {|in, threshold=0.15, wet=0.7, cutoff=12000, i_fftsize=1024|
        /*
        Reduce unstable components (noise) in the signal
        -------------------------------------------------

        threshold: 0-1. What to consider noise.
                   near to 0 cuts more noise, near to 1 lets more noise through
        wet      : mix the denoised signal with the original to mask artifacts.
                   A delay is applied to the original to compensate the delay of the FFT
        i_fftsize: the size of the fft. Bigger FFTs are more accurate but incur in
                   time-smearing and most importantly in latency.


        */
        var fftdelay = i_fftsize/SampleRate.ir;
        var source   = in !> BHiCut.ar(_, cutoff);
        var chain    = FFTL(in, i_fftsize);
        var orig     = PV_Copy(chain, LocalBuf(i_fftsize));
        var reduced  = (
            FFTL(in, i_fftsize)
            !> PV_PartialSynthF(_, threshold * (SampleRate.ir * 0.5 / i_fftsize))
            !> IFFT(_)
        );
        var out      = DelayN.ar(in, fftdelay, fftdelay) >< reduced at:wet;
        ^out;
    }
    *new {|in, attendb=12, cutoff=12000, wet=0.7, i_fftsize=1024|
        ^EM_RemoveResidual.ar(in, attendb, cutoff, wet, i_fftsize);
    }
}

EM_CubicDistortion {
    *ar {|in, drive=0.5, offset=0|
        /*
        drive should be smoothed. DC should be blocked

        cubicnl(drive,offset) = *(pregain) : +(offset) : clip(-1,1) : cubic
        with {
            pregain = pow(10.0,2*drive);
            clip(lo,hi) = min(hi) : max(lo);
            cubic(x) = x - x*x*x/3;
            postgain = max(1.0,1.0/pregain);
        };
        */
        var pregain  = pow(10, 2*drive);
        var postgain = max(1.0, 1.0/pregain);
        var cubic = {|x| x - (x*x*x/3)};
        var process = (in*pregain + offset).clip(-1, 1) !> cubic * postgain;
        ^process;
    }
}

EM_Hysteresis {
    // see EM_NoiseGate
    *ar {|in, thresh_att=0.03, thresh_rel=0.01, attack=0.01, hold=0.1, release=0.2, min_rest=0.01|
        var amp = Amplitude.ar(in);
        var t0 = amp > thresh_att !> Trig.ar(_, hold);
        var t1 = amp < thresh_rel * (1-t0) !> Trig.ar(_, min_rest);
        var trig_envelope = SetResetFF.ar(t0, t1);
        var out = EnvGen.ar(Env.adsr(attack, 0, 1, release), gate: trig_envelope);
        ^out;
    }
    *kr {|in, thresh_att=0.03, thresh_rel=0.01, attack=0.01, hold=0.1, release=0.2, min_rest=0.01|
        var amp = Amplitude.kr(in);
        var t0 = amp > thresh_att !> Trig.kr(_, hold);
        var t1 = amp < thresh_rel * (1-t0) !> Trig.kr(_, min_rest);
        var trig_envelope = SetResetFF.kr(t0, t1);
        var out = EnvGen.kr(Env.adsr(attack, 0, 1, release), gate:trig_envelope);
        ^out;
    }
}

EM_NoiseGate {
    /*
    in     : the signal to modify
    control: the signal to use as cotrol (in most cases it will be the same signal, but it could be a filtered signal, etc.)
    thresh_att: the control signal must be higher than this threshold (amp, NOT dB) for an action to begin
    thresh_rel: the release phase is entered when the control is softer than this value and we are past the hold-time
    attack : in seconds
    hold   : in seconds. The signal will be allowed at least hold time to be audible at full amplitude
    release: in seconds
    min_rest: when in silent mode, wait at least this time before any other action can be triggered
    */
    *ar {|in, control, thresh_att=0.06, thresh_rel= 0.01, attack=0.05, hold=0.1, release=0.2, min_rest=0.01|
        var hyst = EM_Hysteresis.ar(control, thresh_att, thresh_rel, attack, hold, release, min_rest);
        var out = in * hyst;
        ^out;
    }
    *kr {|in, control, thresh_att=0.06, thresh_rel= 0.01, attack=0.05, hold=0.1, release=0.2, min_rest=0.01|
        var hyst = EM_Hysteresis.kr(control, thresh_att, thresh_rel, attack, hold, release, min_rest);
        var out = in * hyst;
        ^out;
    }
}


+ UGen {
    oversamp2 {|func|
        /*
        Example
        -------

        { SoundIn.ar(0).oversamp2 { |sig| LOSER_WaveShapeDist.ar(sig) } !> Out([0, 1])
        }.play;

        {
            SoundIn.ar(0).oversamp2 {|sig|
                sig * pow(10, 2*drive) !>_.clip(-1, 1) !> {|sig| sig - (sig*sig*sig/3)}
            }
            !> Out([0, 1])
        )
        */
        var up = UpSample.ar(this);
        var post = func.(up);
        var down = DownSample.ar(post);
        ^down;
    }
}