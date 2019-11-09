ADC {
    *new {|chan|
        ^SoundIn.ar(chan);
    }
}

Stereo {
	*new {|source, chan=0|
		^Out.ar([chan, chan+1], source);
	}
}

+ Object {
    !> {|other, mod|
        if( mod.isNil ) {
            ^other.apply(this);
        };
	    ^switch( mod,
            // parallel, carry the original as first.
            // a !>.par func   --> [a, func(a)]
            // If you want it the other way around:
            // A) use parafter  --> a !>.parafter b  --> [b(a), a]
            // B) a !>.par b !>_.reverse

            \par, { [this, this !> other] },
            // a shortcut for the above: a !>.parafter func --> [b(a), a]
            \parafter, { [this !> other, this] },
            // side chain:  var source = SinOsc.ar(440) !>.side SendPeakRMS.ar(_, cmdName:'/peaks');
            // now source contains the SinOsc, the SendPeakRMS is collected in the UGen graph (like all UGens, which have no output)
            // this is similar to
            // A)  SinOsc.ar(440) !2 !< SendpeakRMS.ar(_) at: 1
            // B)  SinOsc.ar(440) !>.par SendPeakRMS.ar(_) at: 0
            \side,   { other.apply(this); this },
            // a !>.xfade b(_, 1000) at:0.8  ---> a xfade: b(a, 1000) @0.8
            // this could also be achieved: WhiteNoise.ar !>.par LPF.ar(_, 1000) !> SelectX.ar(0.8, _)
            \xfade, { this xfade: ( this !> other ) },
            // a !>.blend b at:c  ---->  a >< (a !> b) @c
            // also:   a !>< b at:c -->  a >< (a !> b) @c
            \blend, { this >< (this !> other) },
            // a_ugen + b_ugen !>.beam \sum !> keep acting on this
            // \sum.unbeam.printosc("sum")  // monitor this part of the signal
            \beam, { this.beam(other) },
            // shortcut to this !> other !> _[index], now can be written this !>.0 other(_)
            // see !><
            0, { other.apply(this)[0] },
            1, { other.apply(this)[1] }
        )
	}
	!< {|other|
	    ^other.applyOneToMany(this);
	}

    <! {|other|
        ^( other !> this );
    }

    >! {|other|
        ^( other !< this );
    }

    !>< {|other, mod|
        if( mod.isNil ) {
            ^( this >< (this !> other) );
        };
        // for the case where other expands to stereo --> sig !><.0 GVerb.ar(_) // this will output the channel 0 of GVerb
        ^switch(mod,
            0, { this >< ((this !> other)[0]) },
            1, { this >< ((this !> other)[1]) }
        );
    }

    mix {
	    ^(Mix(this));
    }

	arity {
		^0
	}
	apply {|x|
		^this.value(x)
	}

	par {|other|
		^[this, other]
	}

    ++ {|other|
        if(isKindOf(other, SequenceableCollection)) {
            ^[this] ++ other;
        }
        ^[this, other];
    }
    beam {|other|  // other should be a symbol --> see Symbol::unbeam
        // a temporary storage to send and receive signals
        /*
        {
        SoundIn.ar(0).beam("orig") !> LPF.ar(_, 440) !> ... do something with this

        Out.ar(0, \orig.unbeam);  // the original signal is here.

        See Also !>.beam
        */
        Elib.registry[other.asSymbol] = this;
        ^this;
    }
}

+ Symbol {
    unbeam { ^Elib.registry[this] }  // see Object::beam
}

+ SequenceableCollection {
	!> { |other, which|
		var out;
		if( which.isNil ) {
			^other.apply(this)
		} {
			out = other.apply(this[which]);
			^(this.copy.slotPut(which, out))
		}
	}
	!< { |other| ^other.apply_distribute(*this) }

	apply_distribute{ |xs|
	    "error".postln
	}

	applyOneToMany{ |obj|
	    if( obj.isKindOf(SequenceableCollection) ){
	       ^this.apply_distribute(obj);
	    }{  // else
	       ^this.collect{ |x| x.apply(obj) };
	    };
	}
}

+BinaryOpFunction {
	arity {
		^(this.slotAt(this.slotIndex(\a)).arity +
		  this.slotAt(this.slotIndex(\b)).arity
		)
	}
}

+UnaryOpFunction {
	arity {
		^(this.slotAt(this.slotIndex(\a)).arity)
	}
}

+ Function {
	arity {
		^this.def.argNames.size
	}
	doesChannelExpansion {
		if( this.def.selectors.isNil, {^false});
		^this.def.selectors.any
			{|x| x.asClass.superclasses.includes(UGen)}
	}
}

+ AbstractFunction {
    apply {|obj|
        switch( this.arity,
            1, { ^(this.value(obj)) },
            2, { ^{|a| this.value(obj, a)} },
            3, { ^{|a, b| this.value(obj, a, b)} },
            4, { ^{|a, b, c| this.value(obj, a, b, c)} },
            5, { ^{|a, b, c, d| this.value(obj, a, b, c, d)} },
        );
        ^{ |...args| this.values(*([obj] ++ args)) };
    }

	apply_distribute {|...seq|
		var tmp, arity;
		arity = this.arity;
		if( seq.size < arity ) {
			^{|...args| this.value(*(seq ++ args))}
		};
		if( seq.size == arity ) {
			^this.value(*seq);
		} /* else: seq.size > arity */ {
			// ^this.value(seq);
			^( [this.value(*seq[..arity])] ++ seq[arity..] );
		};
	}

	applyOneToMany { |obj|
	    ^this.value(*Array.fill(this.arity, obj));
	}
}

+ Mix {
    *apply {|array|
        ^Mix(array);

    }
}

+ Out {
    *new {|chan|
        if( chan.isKindOf(Bus) ) {
            if( chan.rate == 'audio' ) {
                ^Out.ar(chan.index, _);
            } /* else */ {
                ^Out.kr(chan.index, _);
            }
        } /* else */ {
			^{|x|
				if( x.rate == 'audio') {
					Out.ar(chan, x);
				} {
					Out.kr(chan, x);
				}
			}
			// ^Out.ar(chan, _);
        }
    }
}

+ Symbol {
	ndef {|func|
		^Ndef(this, func)
	}
}

+ ListPattern {
	* apply {|seq|
        ^this.new(seq);
    }
}

// TODO: fix channel expansion in UGens
