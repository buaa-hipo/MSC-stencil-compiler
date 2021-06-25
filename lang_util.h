
// Regular binary ops:
// Operations that take two oprands, and returns a single operand with the same
// type

enum class BinaryOpType : int {
  mul,
  add,
  sub,
  truediv, // will be lower into div in type checking
  floordiv, // will be lower into div in type checking
  div,
  mod,
  max,
  min,
  bit_and,
  bit_or,
  bit_xor,
  cmp_lt,
  cmp_le,
  cmp_gt,
  cmp_ge,
  cmp_eq,
  cmp_ne,
  atan2,
  undefined
};

enum class UnaryOpType : int {
  neg,
  sqrt,
  floor,
  ceil,
  cast,
  abs,
  sgn,
  sin,
  asin,
  cos,
  acos,
  tan,
  tanh,
  inv,
  rcp,
  exp,
  log,
  rsqrt,
  bit_not,
  logic_not,
  undefined
};


enum class DataType : int {
  i1,
  i8,
  i16,
  i32,
  i64,
  f16,
  f32,
  f64,
  u8,
  u16,
  u32,
  u64,
  ptr,
  none,  // "void"
  unknown
};




