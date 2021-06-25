
#undef DEFINE_EXPRESSION_OP_BINARY
#undef DEFINE_EXPRESSION_OP_UNARY
#undef DEFINE_EXPRESSION_FUNC

#define DEFINE_EXPRESSION_OP_UNARY(opname)                                                  \
  ExprPtr opname(ExprPtr expr) {                                                            \
    auto _Una = std::make_shared<UnaryOpExpr>(UnaryOpType::opname, expr );                  \
    ExprPtr Una = std::dynamic_pointer_cast<Expr>(_Una);                                    \
    return Una;                                                                             \
  }                                                                        


#define DEFINE_EXPRESSION_OP_BINARY(op, opname)                                             \
  ExprPtr operator op(ExprPtr lhs, ExprPtr rhs) {                                           \
    auto _Bin = std::make_shared<BinaryOpExpr>(BinaryOpType::opname, lhs, rhs);             \
    ExprPtr Bin = std::dynamic_pointer_cast<Expr>(_Bin);                                    \
    return Bin;                                                                             \
  }                                                                        

#define DEFINE_EXPRESSION_FUNC(opname)                                                      \
  ExprPtr opname(ExprPtr lhs, ExprPtr rhs) {                                                \
    auto _Fun = std::make_shared<FuncOpExpr>(BinaryOpType::opname, lhs, rhs);             \
    ExprPtr Fun = std::dynamic_pointer_cast<Expr>(_Fun);                                    \
    return Fun;                                                                             \
  }                                                                        

#define DEFINE_EXPRESSION_OP_BINARY_CONST_RHS(op, opname)                                   \
  ExprPtr operator op(ExprPtr lhs, int i) {                                                 \
    ExprPtr rhs = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>((int)i));  \
    auto _Bin = std::make_shared<BinaryOpExpr>(BinaryOpType::opname, lhs, rhs);             \
    ExprPtr Bin = std::dynamic_pointer_cast<Expr>(_Bin);                                    \
    return Bin;                                                                             \
  }                                                                        

#define DEFINE_EXPRESSION_OP_BINARY_CONST_LHS(op, opname)                                   \
  ExprPtr operator op(int i, ExprPtr rhs) {                                                 \
    ExprPtr lhs = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>((int)i));  \
    auto _Bin = std::make_shared<BinaryOpExpr>(BinaryOpType::opname, lhs, rhs);             \
    ExprPtr Bin = std::dynamic_pointer_cast<Expr>(_Bin);                                    \
    return Bin;                                                                             \
  }                                                                        



DEFINE_EXPRESSION_OP_UNARY(sqrt)
DEFINE_EXPRESSION_OP_UNARY(floor)
DEFINE_EXPRESSION_OP_UNARY(abs)
DEFINE_EXPRESSION_OP_UNARY(sin)
DEFINE_EXPRESSION_OP_UNARY(asin)
DEFINE_EXPRESSION_OP_UNARY(cos)
DEFINE_EXPRESSION_OP_UNARY(acos)
DEFINE_EXPRESSION_OP_UNARY(tan)
DEFINE_EXPRESSION_OP_UNARY(tanh)
DEFINE_EXPRESSION_OP_UNARY(inv)
DEFINE_EXPRESSION_OP_UNARY(rcp)
DEFINE_EXPRESSION_OP_UNARY(rsqrt)
DEFINE_EXPRESSION_OP_UNARY(exp)
DEFINE_EXPRESSION_OP_UNARY(log)

//
DEFINE_EXPRESSION_OP_BINARY_CONST_RHS(+, add)
DEFINE_EXPRESSION_OP_BINARY_CONST_LHS(+, add)
DEFINE_EXPRESSION_OP_BINARY_CONST_RHS(-, sub)
DEFINE_EXPRESSION_OP_BINARY_CONST_LHS(-, sub)
DEFINE_EXPRESSION_OP_BINARY_CONST_RHS(*, mul)
DEFINE_EXPRESSION_OP_BINARY_CONST_LHS(*, mul)
//

DEFINE_EXPRESSION_OP_BINARY(+, add)
DEFINE_EXPRESSION_OP_BINARY(-, sub)
DEFINE_EXPRESSION_OP_BINARY(*, mul)
DEFINE_EXPRESSION_OP_BINARY(/, div)
DEFINE_EXPRESSION_OP_BINARY(%, mod)
DEFINE_EXPRESSION_OP_BINARY(&&, bit_and)
DEFINE_EXPRESSION_OP_BINARY(||, bit_or)
// DEFINE_EXPRESSION_OP_BINARY(&, bit_and)
// DEFINE_EXPRESSION_OP_BINARY(|, bit_or)
DEFINE_EXPRESSION_OP_BINARY (^, bit_xor)
DEFINE_EXPRESSION_OP_BINARY(<, cmp_lt)
DEFINE_EXPRESSION_OP_BINARY(<=, cmp_le)
DEFINE_EXPRESSION_OP_BINARY(>, cmp_gt)
DEFINE_EXPRESSION_OP_BINARY(>=, cmp_ge)
DEFINE_EXPRESSION_OP_BINARY(==, cmp_eq)
DEFINE_EXPRESSION_OP_BINARY(!=, cmp_ne)

DEFINE_EXPRESSION_FUNC(min);
DEFINE_EXPRESSION_FUNC(max);
DEFINE_EXPRESSION_FUNC(atan2);
DEFINE_EXPRESSION_FUNC(truediv);
DEFINE_EXPRESSION_FUNC(floordiv);

#undef DEFINE_EXPRESSION_FUNC
#undef DEFINE_EXPRESSION_OP_UNARY
#undef DEFINE_EXPRESSION_OP_BINARY
