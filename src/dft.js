/*===========================================================================*\
 * Discrete Fourier Transform (O(n^2) brute-force method)
 *
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Fast Fourier Transform.
\*===========================================================================*/

//------------------------------------------------
// Note: this code is not optimized and is
// primarily designed as an educational and testing
// tool.
//------------------------------------------------
const complex = require("./complex");
const fftUtil = require("./fftutil");

//-------------------------------------------------
// Calculate brute-force O(n^2) DFT for vector.
//-------------------------------------------------
const dft = (vector) => {
  const X = [],
    N = vector.length;

  for (let k = 0; k < N; k++) {
    X[k] = [0, 0]; //Initialize to a 0-valued complex number.

    for (let i = 0; i < N; i++) {
      const exp = fftUtil.exponent(k * i, N);
      let term;
      if (Array.isArray(vector[i])) term = complex.multiply(vector[i], exp);
      //If input vector contains complex numbers
      else term = complex.multiply([vector[i], 0], exp); //Complex mult of the signal with the exponential term.
      X[k] = complex.add(X[k], term); //Complex summation of X[k] and exponential
    }
  }

  return X;
};

module.exports = dft;
