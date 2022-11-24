//-------------------------------------------------
// Add two complex numbers
//-------------------------------------------------
const add = (a, b) => {
  return [a[0] + b[0], a[1] + b[1]];
};

//-------------------------------------------------
// Subtract two complex numbers
//-------------------------------------------------
const subtract = (a, b) => {
  return [a[0] - b[0], a[1] - b[1]];
};

//-------------------------------------------------
// Multiply two complex numbers
//
// (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
//-------------------------------------------------
const multiply = (a, b) => {
  return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
};

//-------------------------------------------------
// Calculate |a + bi|
//
// sqrt(a*a + b*b)
//-------------------------------------------------
const magnitude = (c) => {
  return Math.sqrt(c[0] * c[0] + c[1] * c[1]);
};

//-------------------------------------------------
// Exports
//-------------------------------------------------
module.exports = {
  add,
  subtract,
  multiply,
  magnitude,
};
