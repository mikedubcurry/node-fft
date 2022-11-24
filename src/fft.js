/*===========================================================================*\
 * Fast Fourier Transform (Cooley-Tukey Method)
 *
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Fast Fourier Transform.
\*===========================================================================*/

//------------------------------------------------
// Note: Some of this code is not optimized and is
// primarily designed as an educational and testing
// tool.
// To get high performace would require transforming
// the recursive calls into a loop and then loop
// unrolling. All of this is best accomplished
// in C or assembly.
//-------------------------------------------------

//-------------------------------------------------
// The following code assumes a complex number is
// an array: [real, imaginary]
//-------------------------------------------------
const complex = require("./complex");
const fftUtil = require("./fftutil");
const twiddle = require("bit-twiddle");

module.exports = {
  //-------------------------------------------------
  // Calculate FFT for vector where vector.length
  // is assumed to be a power of 2.
  //-------------------------------------------------
  fft: (vector) => {
    const X = [];
    const N = vector.length;

    // Base case is X = x + 0i since our input is assumed to be real only.
    if (N == 1) {
      if (Array.isArray(vector[0]))
        //If input vector contains complex numbers
        return [[vector[0][0], vector[0][1]]];
      else return [[vector[0], 0]];
    }

    // Recurse: all even samples
    const X_evens = this.fft(vector.filter(even));
    // Recurse: all odd samples
    const X_odds = this.fft(vector.filter(odd));

    // Now, perform N/2 operations!
    for (let k = 0; k < N / 2; k++) {
      // t is a complex number!
      const t = X_evens[k],
        e = complex.multiply(fftUtil.exponent(k, N), X_odds[k]);

      X[k] = complex.add(t, e);
      X[k + N / 2] = complex.subtract(t, e);
    }

    function even(__, ix) {
      return ix % 2 == 0;
    }

    function odd(__, ix) {
      return ix % 2 == 1;
    }

    return X;
  },
  //-------------------------------------------------
  // Calculate FFT for vector where vector.length
  // is assumed to be a power of 2.  This is the in-
  // place implementation, to avoid the memory
  // footprint used by recursion.
  //-------------------------------------------------
  fftInPlace: (vector) => {
    const N = vector.length;

    const trailingZeros = twiddle.countTrailingZeros(N); //Once reversed, this will be leading zeros

    // Reverse bits
    for (let k = 0; k < N; k++) {
      const p = twiddle.reverse(k) >>> (twiddle.INT_BITS - trailingZeros);
      if (p > k) {
        const complexTemp = [vector[k], 0];
        vector[k] = vector[p];
        vector[p] = complexTemp;
      } else {
        vector[p] = [vector[p], 0];
      }
    }

    //Do the DIT now in-place
    for (let len = 2; len <= N; len += len) {
      for (let i = 0; i < len / 2; i++) {
        const w = fftUtil.exponent(i, len);
        for (let j = 0; j < N / len; j++) {
          const t = complex.multiply(w, vector[j * len + i + len / 2]);
          vector[j * len + i + len / 2] = complex.subtract(
            vector[j * len + i],
            t
          );
          vector[j * len + i] = complex.add(vector[j * len + i], t);
        }
      }
    }
  },
};
