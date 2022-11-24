/*===========================================================================*\
 * Fast Fourier Transform Frequency/Magnitude passes
 *
 * (c) Vail Systems. Joshua Jung and Ben Bryan. 2015
 *
 * This code is not designed to be highly optimized but as an educational
 * tool to understand the Fast Fourier Transform.
\*===========================================================================*/

//-------------------------------------------------
// The following code assumes a complex number is
// an array: [real, imaginary]
//-------------------------------------------------
const complex = require("./complex");

//-------------------------------------------------
// By Eulers Formula:
//
// e^(i*x) = cos(x) + i*sin(x)
//
// and in DFT:
//
// x = -2*PI*(k/N)
//-------------------------------------------------
const mapExponent = {};
const exponent = (k, N) => {
  const x = -2 * Math.PI * (k / N);

  mapExponent[N] = mapExponent[N] || {};
  mapExponent[N][k] = mapExponent[N][k] || [Math.cos(x), Math.sin(x)]; // [Real, Imaginary]

  return mapExponent[N][k];
};

//-------------------------------------------------
// Calculate FFT Magnitude for complex numbers.
//-------------------------------------------------
const fftMag = (fftBins) => {
  const ret = fftBins.map(complex.magnitude);
  return ret.slice(0, ret.length / 2);
};

//-------------------------------------------------
// Calculate Frequency Bins
//
// Returns an array of the frequencies (in hertz) of
// each FFT bin provided, assuming the sampleRate is
// samples taken per second.
//-------------------------------------------------
const fftFreq = (fftBins, sampleRate) => {
  const stepFreq = sampleRate / fftBins.length;
  const ret = fftBins.slice(0, fftBins.length / 2);

  return ret.map((__, ix) => {
    return ix * stepFreq;
  });
};

//-------------------------------------------------
// Exports
//-------------------------------------------------
module.exports = {
  fftMag,
  fftFreq,
  exponent,
};
