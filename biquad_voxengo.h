/**
 * @file biquad_voxengo.h
 *
 * @version 1.2
 *
 * @brief Perfect biquad filter design code.
 *
 * Email: aleksey.vaneev@gmail.com or info@voxengo.com
 *
 * LICENSE:
 *
 * Copyright (c) 2026 Aleksey Vaneev
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef BIQUAD_VOXENGO
#define BIQUAD_VOXENGO

#include <math.h>

typedef struct {
    double b0, b1, b2, a0, a1, a2;
} biquad_t;

enum biquadType { BT_PEQ = 0, BT_BPF };

/**
 * @brief Designs a bell, notch, or band-pass filter.
 *
 * @param type Filter type.
 * @param SampleRate Sampling rate (e.g., 48000).
 * @param Freq Center frequency in Hz.
 * @param Gain Bell or notch filter's gain in linear scale (e.g., 2 for +6 dB);
 * ignored for band-pass filter.
 * @param BW Filter's bandwidth in octaves. `BW==2/ln(2)*asinh(1/(2*Q))`.
 * @param[out] f Resulting biquad filter coefficients.
 */

static inline void cookBiquadVoxengo( const int type, const double SampleRate,
    const double Freq, const double Gain, const double BW, biquad_t* const f )
{
    double Fp, Fb, xp, xb, y, v2, w, r1, r2, den, A2, B2, A, B;
    double g0, gb, gp, gn, G0w, Gn, r, t, u;

    // ---- normalized frequencies (tan blows up without clamp) ----
    Fp = Freq / SampleRate;
    if( Fp < 1e-9 ) Fp = 1e-9; else if( Fp > 0.4999999 ) Fp = 0.4999999;
    Fb = Fp * pow( 0.5, BW * 0.5 );

    // ---- normalized band-edge detuning ----
    // r = 1 exactly when the octave band's upper edge reaches Nyquist.
    // y = r^2 is the response-skew DOF: it fixes the Nyquist anchor gain.
    // `rs` is a shift parameter with 2.0 yielding an intended design;
    // 1.7 shifts Gn and produces a better analog prototype match.
    const double rs = 2.0;
    r = ( Fp*Fp - Fb*Fb ) / ( rs * Fb * ( 0.25 - Fp*Fp ));
    y = r * r;

    // ---- family anchors, squared gains ----
    if( type == BT_PEQ )
    {
        if( fabs( Gain - 1.0 ) < 1e-9 )
        {
            f->a0 = f->b0 = 1.0; f->a1 = f->a2 = f->b1 = f->b2 = 0.0;
            return;
        }

        gn = ( 1.0 + Gain*y ) / ( 1.0 + y/Gain );
        g0 = 1.0;
        gb = Gain;
        gp = Gain * Gain;
        v2 = Gain / ( Gain + y ); // ( gp - gn ) / ( gp - g0 ) equivalent.
    }
    else
    {
        gn = y / ( 1.0 + y );
        g0 = 0.0;
        gb = 0.5;
        gp = 1.0;
        v2 = 1.0 / ( 1.0 + y );
    }

    // ---- warped frequency axis ----
    const double C_PI = 3.14159265358979324;
    xp = tan( C_PI * Fp ); xp *= xp;
    xb = tan( C_PI * Fb ); xb *= xb;

    w = xp * sqrt( v2 );
    Gn = sqrt( gn );
    G0w = sqrt( g0 ) * w;

    // ---- 2x2 linear solve for A^2, B^2 ----
    t = w - xp; u = G0w - Gn*xp;
    r1 = ( gp*t*t - u*u ) / xp;

    t = w - xb; u = G0w - Gn*xb;
    r2 = ( gb*t*t - u*u ) / xb;

    den = gb - gp;

    A2 = ( r1 - r2 ) / den;
    B2 = ( gb*r1 - gp*r2 ) / den;

    A = sqrt( A2 );
    B = sqrt( B2 );

    // ---- emit biquad ----
    f->a0 = 1.0 + w + A;
    f->a1 = 2.0 * ( w - 1.0 );
    f->a2 = 1.0 + w - A;
    f->b0 = Gn + G0w + B;
    f->b1 = 2.0 * ( G0w - Gn );
    f->b2 = Gn + G0w - B;
}

#endif // BIQUAD_VOXENGO
