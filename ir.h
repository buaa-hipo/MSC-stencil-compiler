#include <cassert>
#include <cstdint>
#include <map>
#include <queue>
#include <algorithm>
#include <fstream>
#include "common.h"
#include "lang_util.h"
#include "util.h"
#include "CodeGen.h"
#include "IO.h"
#include "runtime.h"

class IRBuilder;
class IRNode;
class Block;
class Stmt;

class Expression;
class Expr;
class ExprGroup;

class VariableExpr;

class FuncOpExpr;
class AssignOpExpr;
class BinaryOpExpr;
class UnaryOpExpr;
class SpNode;
class SpNodeIndicesExpr; 
class TeNodeExpr; 
class TeNodeIndicesExpr;
class FrontendForStmt;
class RangeForStmt;
class FunctionStmt;
class CallFuncExpr;
class CallNodeExpr;
class FrontendAllocaStmt;
class FrontendForStmt;
class FunctionStmt;
class PlaceHolderExpr; 
class InputExpr;
class OutputExpr;
class MPIExpr;
class ScheduleForStmt;
class AthreadStmt; 

class IRVisitor {
 public:

  IRVisitor() {}

  // default visitor
  virtual void visit(Stmt *stmt) {}
  virtual void visit(FuncOpExpr *op) {}
  virtual void visit(AssignOpExpr *op) {}
  virtual void visit(BinaryOpExpr *op) {}
  virtual void visit(UnaryOpExpr *op) {}
  virtual void visit(VariableExpr *op) {}
  virtual void visit(SpNodeIndicesExpr *op) {}
  virtual void visit(TeNodeExpr * op) {}
  virtual void visit(TeNodeIndicesExpr * op) {}
  virtual void visit(CallFuncExpr * op) {}
  virtual void visit(CallNodeExpr * op) {}
  virtual void visit(FrontendAllocaStmt * op) {}
  virtual void visit(FrontendForStmt * op) {}
  virtual void visit(RangeForStmt * op) {}
  virtual void visit(FunctionStmt * op) {}
  virtual void visit(PlaceHolderExpr * op) {}
  virtual void visit(InputExpr * op) {}
  virtual void visit(OutputExpr * op) {}
  virtual void visit(MPIExpr * op) {}
  virtual void visit(ScheduleForStmt * op) {}
  virtual void visit(AthreadStmt * op) {}

#define DEFINE_VISIT(T)            \
  virtual void visit(T *stmt) {    \
    visit((Stmt *)stmt);       \
  }

  DEFINE_VISIT(Block);

  /*
#define DEFINE_VISIT(T, t)            \
  virtual void visit(T * t) {}        \
  virtual VecInt visit_shape(T * t) {}

  DEFINE_VISIT(Stmt , stmt);
  DEFINE_VISIT(FuncOpExpr, op);
  DEFINE_VISIT(BinaryOpExpr, op);
  DEFINE_VISIT(UnaryOpExpr, op);
  DEFINE_VISIT(VariableExpr, op);
  DEFINE_VISIT(SpNodeIndicesExpr, op);
  DEFINE_VISIT(TeNodeExpr ,op);
  */
};


#define DEFINE_ACCEPT                        \
  void accept(IRVisitor *visitor) override { \
    visitor->visit(this);                    \
  }

class IRNode {
 public: 
  Block * parent;
  IRNode * get_ir_root();
  virtual void accept(IRVisitor *visitor) {}
  virtual ~IRNode() = default;
};

class Stmt : public IRNode {
 public:
  //Block * parent;
  //std::shared_ptr<Block> parent;


  Stmt() {
    parent = nullptr;
  }

  static uint64 operand_hash(Stmt *stmt) {
    return uint64(1) << ((uint64(stmt) >> 4) % 64);
  }

  //IRNode * get_ir_root();

  virtual ~Stmt() override = default;
};



class Block : public IRNode {
 public:
  //std::shared_ptr<Block> parent;
  //Block * parent;
  //std::vector<Stmt *> statements;
  //std::vector<std::shared_ptr<Stmt> > statements;
  std::vector<std::shared_ptr<IRNode> > statements;

  Block() {
    parent = nullptr;
  }

  //int locate(std::shared_ptr<Stmt> stmt) {
  int locate(std::shared_ptr<IRNode> stmt) {
    for (int i = 0; i < (int)statements.size(); i++) {
      if (statements[i] == stmt) {
        return i;
      }
    }
    return -1;
  }

  void erase(int location) {
    statements.erase(statements.begin() + location);
  }
  
  //void erase(std::shared_ptr<Stmt> stmt) {
  void erase(std::shared_ptr<IRNode> stmt) {
    for (int i = 0; i < (int)statements.size(); i++) {
      if (statements[i] == stmt) {
        erase(i);
        break;
      }
    }
  }

  //void insert(std::shared_ptr<Stmt> stmt, int location = -1 ) {
  void insert(std::shared_ptr<IRNode> stmt, int location = -1 ) {
    stmt->parent = this;

    //shared_from_this() should be work at https://stackoverflow.com/questions/11711034/stdshared-ptr-of-this
    //but we make sure that the stmt has been set parent before insert if we use std::shared_ptr<Block> parent !!!

    if (location == -1) {
      statements.push_back(std::move(stmt));
    } else {
      statements.insert(statements.begin() + location, std::move(stmt));
    }
  }


  //void replace_statements_in_range(int start, int end, VecStatement &&stmts);

  //void set_statements( std::vector<std::shared_ptr<Stmt> > stmts ) {
  void set_statements( std::vector<std::shared_ptr<IRNode> > stmts ) {
    statements.clear();
    for (int i = 0; i < (int)stmts.size(); i++) {
      insert(std::move(stmts[i]), i);
    }
  }


  //void replace_with(std::shared_ptr<Stmt> old_statement,
  //                         std::shared_ptr<Stmt> new_statement) {
  void replace_with(std::shared_ptr<IRNode> old_statement,
                           std::shared_ptr<IRNode> new_statement) {
    int location = -1;
    for (int i = 0; i < (int)statements.size(); i++) {
      if (old_statement == statements[i]) {
        location = i;
        break;
      }
    }
    assert(location>-1);
    statements.erase(statements.begin() + location);
    insert(std::move(new_statement), location);
  }

  //void insert_before( std::shared_ptr<Stmt> old_statement, std::shared_ptr<Stmt> new_statement) {
  void insert_before( std::shared_ptr<IRNode> old_statement, std::shared_ptr<IRNode> new_statement) {
    int location = -1;
    for (int i = 0; i < (int)statements.size(); i++) {
      if (old_statement == statements[i]) {
        location = i;
        break;
      }
    }
    assert(location != -1);
    insert(std::move(new_statement), location);
  }

  //void replace_with(Stmt *old_statement,
  //                  VecStatement &new_statements,
  //                  bool replace_usages = true) {
  //  int location = -1;
  //  for (int i = 0; i < (int)statements.size(); i++) {
  //    if (old_statement == statements[i].get()) {
  //      location = i;
  //      break;
  //    }
  //  }
  //  TC_ASSERT(location != -1);
  //  if (replace_usages)
  //    old_statement->replace_with(new_statements.back().get());
  //  trash_bin.push_back(std::move(statements[location]));
  //  statements.erase(statements.begin() + location);
  //  for (int i = (int)new_statements.size() - 1; i >= 0; i--) {
  //    insert(std::move(new_statements[i]), location);
  //  }
  //}

  //Stmt *lookup_var(Ident ident) const;

  //std::shared_ptr<Stmt> mask();
  std::shared_ptr<IRNode> mask();

  //std::shared_ptr<Stmt> back() const {
  std::shared_ptr<IRNode> back() const {
    return statements.back();
  }

  DEFINE_ACCEPT
};

  //IRNode *get_ir_root() {
  //  auto block = parent;
  //  while (block->parent)
  //    block = block->parent;
  //  return dynamic_cast<IRNode *>(block);
  //}

//IRNode * Stmt::get_ir_root() {
IRNode * IRNode::get_ir_root() {
  auto block = parent;
  while (block->parent)
    block = block->parent;
  return dynamic_cast<IRNode *>(block);
}




class FrontendContext {
 private:
  std::shared_ptr<IRBuilder> current_builder;
  std::shared_ptr<Block> root_node;

 public:
  FrontendContext();

  IRBuilder &builder() {
    return *current_builder;
  }

  IRNode *root();

  std::shared_ptr<Block> get_root() {
    return std::move(root_node);
  }
};

FrontendContext::FrontendContext() {
  root_node = std::make_shared<Block>();
  current_builder = std::make_shared<IRBuilder>(root_node);
}

IRNode *FrontendContext::root() {
  auto ptr =std::dynamic_pointer_cast<IRNode>(root_node);
  return ptr.get();
}

//std::shared_ptr<FrontendContext> context = std::make_shared<FrontendContext>();
std::shared_ptr<FrontendContext> context;
std::shared_ptr<FrontendContext> context_slave_sunway;
std::vector<std::pair<int, DataType> > call_extern_func_slave;
//for temporary...
std::map<std::pair<int, int>, std::shared_ptr<Expr> > found;

void init(){
  context = std::make_shared<FrontendContext>();
  context_slave_sunway = std::make_shared<FrontendContext>();
  found.clear();
  call_extern_func_slave.clear();
}


IRBuilder &current_ast_builder() {
  return context->builder();
}

IRBuilder &slave_sunway_ast_builder() {
  return context_slave_sunway->builder();
}

class IRBuilder {
 private:
  //std::vector<Block *> stack;
  std::vector<std::shared_ptr<Block> > stack;

 public:
  IRBuilder(std::shared_ptr<Block> initial) {
    stack.push_back(initial);
  }

  //void insert(std::shared_ptr<Stmt> stmt, int location = -1);
  void insert(std::shared_ptr<IRNode> stmt, int location = -1);

  std::shared_ptr<Block> current_block() {
    if (stack.empty())
      return nullptr;
    else
      return stack.back();
  }

  //std::shared_ptr<Stmt> get_last_stmt();
  std::shared_ptr<IRNode> get_last_stmt();

};


//void IRBuilder::insert(std::shared_ptr<Stmt> stmt, int location) {
void IRBuilder::insert(std::shared_ptr<IRNode> stmt, int location) {
  assert(!stack.empty());
  stack.back()->insert(std::move(stmt), location);
}

//std::shared_ptr<Stmt> IRBuilder::get_last_stmt() {
std::shared_ptr<IRNode> IRBuilder::get_last_stmt() {
  return stack.back()->back();
}

class FrontendAllocaStmt : public Stmt {
  public:
    VecInt _shape;
    DataType _dt;
    int id;
    std::shared_ptr<SpNode> SpNodePtr;
    std::shared_ptr<VariableExpr> VariableExprPtr;

    FrontendAllocaStmt(VecInt shape, DataType dt, int id) {
      this->_shape = shape;
      this->_dt = dt;
      this->id = id;
      this->SpNodePtr = nullptr;
      this->VariableExprPtr = nullptr;
    }

    FrontendAllocaStmt(VecInt shape, DataType dt, int id, std::shared_ptr<SpNode> SpNodePtr) {
      this->_shape = shape;
      this->_dt = dt;
      this->id = id;
      this->SpNodePtr = SpNodePtr;
      this->VariableExprPtr = nullptr;
    }

    FrontendAllocaStmt(VecInt shape, DataType dt, int id, std::shared_ptr<VariableExpr> VariableExprPtr) {
      this->_shape = shape;
      this->_dt = dt;
      this->id = id;
      this->SpNodePtr = nullptr;
      this->VariableExprPtr = VariableExprPtr;
    }

    bool alloc_spnode(){
      return this->SpNodePtr != nullptr;
    }

    DataType dt(){
      return this->_dt;
    }
    VecInt shape(){
      return this->_shape;
    }

    ~FrontendAllocaStmt() {}

  DEFINE_ACCEPT
};

/*
class AllocaStmt : public Stmt {
 public:
  AllocaStmt(DataType type) {
    ret_type = VectorType(1, type);
  }

  AllocaStmt(int width, DataType type) {
    ret_type = VectorType(width, type);
  }

  virtual bool has_global_side_effect() const override {
    return false;
  }

  DEFINE_ACCEPT
};
*/

// expr.* file

class Expr : public IRNode {
  public:
    DataType _dt;
    VecInt _shape;
    DataType dt(){
      return this->_dt;
    }
    void set_datatype(DataType dt){
      this->_dt = dt;
    }
    VecInt shape(){
      return this->_shape;
    }
    void set_shape(VecInt shape){
      this->_shape = shape;
    }
    Expr() {}
    ~Expr() {}
};

typedef std::shared_ptr<Expr> ExprPtr;
#include "./expression.h"

class ExprGroup {
 public:
  std::vector<std::shared_ptr<Expr> > exprs;

  ExprGroup() {
  }

  ExprGroup(const ExprPtr a) {
    exprs.push_back(a);
  }

  ExprGroup(const ExprPtr a, const ExprPtr b) {
    exprs.push_back(a);
    exprs.push_back(b);
  }

  ExprGroup(ExprGroup a, const ExprPtr b) {
    exprs = a.exprs;
    exprs.push_back(b);
  }

  ExprGroup(const ExprPtr a, ExprGroup b) {
    exprs = b.exprs;
    exprs.insert(exprs.begin(), a);
  }

  void push_back(const ExprPtr expr) {
    exprs.emplace_back(expr);
  }

  std::size_t size() const {
    return exprs.size();
  }

  bool empty() const {
    return exprs.empty();
  }

  const ExprPtr operator[](int i) const {
    assert(0<=i && i< exprs.size());
    return exprs[i];
  }

  ExprPtr operator[](int i) {
    assert(0<=i && i< exprs.size());
    return exprs[i];
  }

  /*
  std::string serialize() {
    std::string ret;
    for (int i = 0; i < (int)exprs.size(); i++) {
      ret += exprs[i].serialize();
      if (i + 1 < (int)exprs.size()) {
        ret += ", ";
      }
    }
    return ret;
  }
  */

  ExprGroup loaded() const;
};


inline ExprGroup operator,(const ExprPtr a, const ExprPtr b) {
  fmt::print(" from operator , (ExprPtr, ExprPtr)\n");
  return ExprGroup(a, b);
}

inline ExprGroup operator,(const ExprGroup &a, const ExprPtr b) {
  fmt::print(" from operator , (ExprGroup, ExprPtr)\n");
  return ExprGroup(a, b);
}


class VariableExpr : public Expr {
  //private:
  public:
    static int counter;
    int id;
    bool is_variable;
    bool is_const_var;
    bool is_immediate;
    bool is_assigned;
    float64 value_f64;
    float32 value_f32;
    int32 value_i32;

  public:

    VariableExpr(){
      set_datatype(DataType::none);
      set_shape(SHAPE_SCALAR);
      is_variable=false;
      is_const_var=false;
      is_immediate=false;
      is_assigned=false;
      this->id = this->counter++;
    }
    VariableExpr(DataType dt){
      set_datatype(dt);
      set_shape(SHAPE_SCALAR);
      is_variable=true;
      is_const_var=false;
      is_immediate=false;
      is_assigned=false;
      this->id = this->counter++;
    }
    //
    VariableExpr(DataType dt, float64 immediate){
      assert(dt == DataType::f64);
      set_datatype(dt);
      set_shape(SHAPE_SCALAR);
      is_variable=false;
      is_const_var=false;
      is_immediate=false;
      is_assigned=true;
      value_f64=immediate;
      this->id = this->counter++;
    }
    VariableExpr(DataType dt, float32 immediate){
      assert(dt == DataType::f32);
      set_datatype(dt);
      set_shape(SHAPE_SCALAR);
      is_variable=false;
      is_const_var=false;
      is_immediate=false;
      is_assigned=true;
      value_f32=immediate;
      this->id = this->counter++;
    }
    VariableExpr(DataType dt, int32 immediate){
      assert(dt == DataType::i32);
      set_datatype(dt);
      set_shape(SHAPE_SCALAR);
      is_variable=false;
      is_const_var=false;
      is_immediate=false;
      is_assigned=true;
      value_i32=immediate;
      this->id = this->counter++;
    }
    //
    VariableExpr(float64 immediate){
      set_datatype(DataType::f64);
      set_shape(SHAPE_SCALAR);
      is_variable=false;
      is_const_var=false;
      is_immediate=true;
      is_assigned=false;
      value_f64=immediate;
      this->id = this->counter++;
    }
    VariableExpr(float32 immediate){
      set_datatype(DataType::f32);
      set_shape(SHAPE_SCALAR);
      is_variable=false;
      is_const_var=false;
      is_immediate=true;
      is_assigned=false;
      value_f32=immediate;
      this->id = this->counter++;
    }
    VariableExpr(int32 immediate){
      set_datatype(DataType::i32);
      set_shape(SHAPE_SCALAR);
      is_variable=false;
      is_const_var=false;
      is_immediate=true;
      is_assigned=false;
      value_i32=immediate;
      this->id = this->counter++;
    }

    ~VariableExpr(){
      set_datatype(DataType::none);
      set_shape(SHAPE_SCALAR);
      is_variable=false;
      is_const_var=false;
      is_immediate=false;
      is_assigned=false;
    }

  DEFINE_ACCEPT
};

int VariableExpr::counter = 0;

class SpNode{
  public:
    static const int default_temp_size;
    int _temp_size;
    size_t _ndim;
    VecInt _shape;
    //VecInt _halo_size;
    int _halo_size;
    int _time_window_size;
    static int counter;
    int id;
    DataType _dt;
    SpNode(){
      _ndim=0;
      _temp_size = default_temp_size;
    }
    SpNode(size_t ndim, VecInt shape, DataType dt){
      _ndim=ndim;
      _shape=shape;
      _dt=dt;
      _temp_size = default_temp_size;
      this->id = this->counter++;
    }
    int halo_size(){
      return this->_halo_size;
    }
    void set_halo_size(int halo_size){
      assert( 0 <= halo_size);
      this->_halo_size = halo_size;
    }
    void set_time_window_size(int time_window_size) {
      this->_time_window_size = time_window_size;
    }
    int time_window_size(){
      return this->_time_window_size;
    }
    bool use_time_window(){
      return this->_time_window_size > 0 ;
    }
    int temp_size(){
      return this->_temp_size;
    }
    void set_temp_size(int temp_size){
      this->_temp_size = temp_size;
    }
    size_t ndim(){
      return this->_ndim;
    }
    DataType dt(){
      return this->_dt;
    }
    VecInt shape(){
      return this->_shape;
    }
    ~SpNode(){
      this->_ndim=0;
      _temp_size = -1;
    }
};

int SpNode::counter = 0;
const int SpNode::default_temp_size = 1000;

class SpNodeIndicesExpr : public Expr {
  public:
    std::shared_ptr<SpNode> SpNodePtr;
    //std::shared_ptr<ExprGroup> ExprGroupPtr
    ExprGroup eg;
    std::shared_ptr<Expr> ArgTempPtr;

    void set_arg_temp(std::shared_ptr<Expr> ArgTempPtr) {
      this->ArgTempPtr = ArgTempPtr;
    }
    void plus_arg_temp(std::shared_ptr<Expr> _ArgTempPtr) {
      std::cout<<"from plus_arg_temp"<<std::endl;
      this->ArgTempPtr = _ArgTempPtr + this->ArgTempPtr;

      //auto _Bin = std::make_shared<BinaryOpExpr>(BinaryOpType::add, _ArgTempPtr, this->ArgTempPtr);             
      //ExprPtr Bin = std::dynamic_pointer_cast<Expr>(_Bin);                                    
      //this->ArgTempPtr=Bin;
    }

    SpNodeIndicesExpr(){}

    //SpNodeIndicesExpr(std::shared_ptr<SpNode> SpNodePtr, std::shared_ptr<ExprGroup> ExprGroupPtr): SpNodePtr(SpNodePtr),
    //  ExprGroupPtr(ExprGroupPtr){}
    SpNodeIndicesExpr(std::shared_ptr<SpNode> SpNodePtr, ExprGroup eg): SpNodePtr(SpNodePtr),eg(eg){
      set_datatype(SpNodePtr->dt());
      set_shape(SpNodePtr->shape());
    }

    ~SpNodeIndicesExpr(){}

  DEFINE_ACCEPT
};

class FrontendTensor {
  public:
    std::shared_ptr<SpNode> SpNodePtr;
    //std::shared_ptr<ExprGroup> ExprGroupPtr;

    FrontendTensor(){}

    FrontendTensor(VecInt shape, DataType dt){
      this->SpNodePtr = std::make_shared<SpNode>(shape.size(), shape, dt);
    }
    
    std::shared_ptr<Expr> operator[](ExprGroup eg){
      assert(this->SpNodePtr->ndim() == eg.size());
      std::shared_ptr<SpNodeIndicesExpr> SpNodeIndicesExprPtr = std::make_shared<SpNodeIndicesExpr>(this->SpNodePtr, eg);
      return std::dynamic_pointer_cast<Expr>(SpNodeIndicesExprPtr);
    }
    
    ~FrontendTensor(){}
};

class TeNode{
  public:
    //the value of temporal is only for analysis of 
    //data dependancy in time dimension. 
    int _temporal;
    size_t _ndim;
    VecInt _shape;
    static int counter;
    int id;
    DataType _dt;
    TeNode(){
      _temporal=0;
      _ndim=0;
    }
    TeNode(int temporal, size_t ndim, VecInt shape, DataType dt){
      _temporal=temporal;
      _ndim=ndim;
      _shape=shape;
      _dt=dt;
      this->id = this->counter++;
    }
    size_t ndim(){
      return _ndim;
    }
    void set_temporal(int temporal){
      _temporal = temporal;
    }
    int temporal(){
      return _temporal;
    }
    VecInt shape(){
      return this->_shape;
    }
    DataType dt(){
      return this->_dt;
    }
    ~TeNode(){
      _temporal=0;
      _ndim=0;
    }
};

int TeNode::counter = 0;

class TeNodeExpr : public Expr {
  public:
    std::shared_ptr<TeNode> TeNodePtr;
    // CallFunc to kernel
    std::shared_ptr<FunctionStmt> FunctionStmtPtr;
    //

    TeNodeExpr(){}

    TeNodeExpr(int temporal, VecInt shape, DataType dt){
      this->TeNodePtr = std::make_shared<TeNode>(temporal, shape.size(), shape, dt);
      set_datatype(TeNodePtr->dt());
      set_shape(TeNodePtr->shape());
    }
    TeNodeExpr(std::shared_ptr<TeNode> TeNodePtr): TeNodePtr(TeNodePtr){
      set_datatype(TeNodePtr->dt());
      set_shape(TeNodePtr->shape());
    }

    ~TeNodeExpr(){}

  DEFINE_ACCEPT
};

class TeNodeIndicesExpr : public Expr {
  public:
    std::shared_ptr<TeNode> TeNodePtr;
    ExprGroup eg;

    TeNodeIndicesExpr(){}

    TeNodeIndicesExpr(std::shared_ptr<TeNode> TeNodePtr, ExprGroup eg): TeNodePtr(TeNodePtr),eg(eg){
      set_datatype(TeNodePtr->dt());
      set_shape(TeNodePtr->shape());
    }

    ~TeNodeIndicesExpr(){}

  DEFINE_ACCEPT
};

typedef std::shared_ptr<VariableExpr> VariableExprPtr;

class FuncOpExpr : public Expr {
 public:
  BinaryOpType type;
  ExprPtr lhs;
  ExprPtr rhs;

  FuncOpExpr(const BinaryOpType &type, ExprPtr lhs, ExprPtr rhs)
      : type(type), lhs(lhs), rhs(rhs) {
  }

  DEFINE_ACCEPT
};

class AssignOpExpr : public Expr {
  public:
    ExprPtr lhs;
    ExprPtr rhs;

    AssignOpExpr(ExprPtr lhs, ExprPtr rhs)
        : lhs(lhs), rhs(rhs) {
    }

    DEFINE_ACCEPT
};

// AssignNode
#define EQUALS <<
ExprPtr operator << (ExprPtr lhs, ExprPtr rhs) {                                           
  return std::dynamic_pointer_cast<Expr>( std::make_shared<AssignOpExpr>(lhs, rhs) );             
}                                                                        
// initialization
ExprPtr operator <<= (ExprPtr lhs, ExprPtr rhs) {                                           
  return std::dynamic_pointer_cast<Expr>( std::make_shared<AssignOpExpr>(lhs, rhs) );             
}                                                                        

class BinaryOpExpr : public Expr {
 public:
  BinaryOpType type;
  ExprPtr lhs;
  ExprPtr rhs;

  BinaryOpExpr(const BinaryOpType &type, ExprPtr lhs, ExprPtr rhs)
      : type(type), lhs(lhs), rhs(rhs) {
  }

  DEFINE_ACCEPT
};


class UnaryOpExpr : public Expr {
 public:
  UnaryOpType type;
  ExprPtr operand;

  UnaryOpExpr(const UnaryOpType &type, ExprPtr operand): type(type), operand(operand) {
  }

  DEFINE_ACCEPT
};

class InputExpr : public Expr {
  public:
    std::vector<std::shared_ptr<Expr> > args;
    std::shared_ptr<SpNode> SpNodePtr;
    std::string name_file;
    bool use_mpi;
    std::vector<int> shape_mpi;
    InputExpr(){}
    InputExpr(std::vector<std::shared_ptr<Expr> > args, std::shared_ptr<SpNode> SpNodePtr, std::string name_file, VecInt shape, DataType dt): 
      args(args), SpNodePtr(SpNodePtr), name_file(name_file) {
        this->set_shape(shape);
        this->set_datatype(dt);
        this->use_mpi = false;
    }
    InputExpr(std::vector<std::shared_ptr<Expr> > args, std::shared_ptr<SpNode> SpNodePtr, std::string name_file, std::vector<int> shape_mpi, VecInt shape, DataType dt): 
      args(args), SpNodePtr(SpNodePtr), name_file(name_file), shape_mpi(shape_mpi) {
        this->set_shape(shape);
        this->set_datatype(dt);
        this->use_mpi = true;
    }

    ~InputExpr(){}

  DEFINE_ACCEPT
};

class OutputExpr : public Expr {
  public:
    std::vector<std::shared_ptr<Expr> > args;
    std::shared_ptr<SpNode> SpNodePtr;
    std::string name_file;
    OutputExpr(){}
    OutputExpr(std::vector<std::shared_ptr<Expr> > args, std::shared_ptr<SpNode> SpNodePtr, std::string name_file, VecInt shape, DataType dt): 
      args(args), SpNodePtr(SpNodePtr), name_file(name_file) {
        this->set_shape(shape);
        this->set_datatype(dt);
    }

    ~OutputExpr(){}

  DEFINE_ACCEPT
};

class MPIExpr : public Expr {
  public:
    std::vector<std::shared_ptr<Expr> > args;
    std::shared_ptr<SpNode> SpNodePtr;
    std::vector<int> shape_mpi;
    //1:asyn_start, 2:asyn_end, 3:syn
    int type_call;
    MPIExpr(){}
    MPIExpr(std::vector<std::shared_ptr<Expr> > args, std::shared_ptr<SpNode> SpNodePtr, std::vector<int> shape_mpi, int type_call, VecInt shape, DataType dt): 
      args(args), SpNodePtr(SpNodePtr), shape_mpi(shape_mpi), type_call(type_call) {
        this->set_shape(shape);
        this->set_datatype(dt);
    }

    ~MPIExpr(){}

  DEFINE_ACCEPT
};

class CallFuncExpr : public Expr {
  public:
    int id_func;
    std::vector<std::shared_ptr<Expr> > args;
    CallFuncExpr(){}
    CallFuncExpr(std::vector<std::shared_ptr<Expr> > args, int id_func, VecInt shape, DataType dt): 
      args(args), id_func(id_func) {
        this->set_shape(shape);
        this->set_datatype(dt);
    }

    ~CallFuncExpr(){}

  DEFINE_ACCEPT
};

class CallNodeExpr : public Expr {
  public:
    std::shared_ptr<TeNode> TeNodePtr;
    CallNodeExpr(){}
    CallNodeExpr( std::shared_ptr<TeNode> TeNodePtr, VecInt shape, DataType dt ): TeNodePtr(TeNodePtr) {
        this->set_shape(shape);
        this->set_datatype(dt);
    }

    ~CallNodeExpr(){}

  DEFINE_ACCEPT
};

class RangeForStmt : public Stmt {
  public:
    std::shared_ptr<Expr> loop_var;
    //int begin;
    std::shared_ptr<Expr> begin;
    //int end;
    std::shared_ptr<Expr> end;
    std::shared_ptr<Block> body;
    int parallelize;
    int vectorize;

    RangeForStmt(){}
    RangeForStmt(std::shared_ptr<Expr> loop_var, 
                  std::shared_ptr<Expr> begin,
                  std::shared_ptr<Expr> end,
                  std::shared_ptr<Block> body,
                  int parallelize,
                  int vectorize)
      : loop_var(loop_var),
        begin(begin),
        end(end),
        body(std::move(body)),
        parallelize(parallelize),
        vectorize(vectorize){}

    ~RangeForStmt(){}

  DEFINE_ACCEPT
};

class FrontendForStmt : public Stmt {
  public:
    ExprGroup eg;
    //VecInt begin_ends;
    VecInt begins;
    VecInt ends;
    VecInt strides;
    std::shared_ptr<Block> body;
    //parallelize;
    int num_threads;

    FrontendForStmt(){}
    //FrontendForStmt(ExprGroup eg, VecInt begin_ends, std::shared_ptr<Block> body): eg(eg), begin_ends(begin_ends), body(body) {}
    //default stride is one
    FrontendForStmt(ExprGroup eg, VecInt begins, VecInt ends, int num_threads, std::shared_ptr<Block> body): eg(eg), begins(begins), ends(ends), num_threads(num_threads), body(body) {
      assert(begins.size() == ends.size());
      strides.clear();
      for(size_t i=0;i<begins.size();i++){
        strides.push_back(1);
      }
    }

    ~FrontendForStmt(){}

  DEFINE_ACCEPT
};

class Axis {
  public:
    //the original id of variable and for loop, which is not splited.
    int id_var;
    //the order of nested for loop with splited or not.
    int order;
    //0:not splited, 1:splited and outer, 2:splited and inner
    int split_flag;
    //
    int begin;
    int end;
    int stride;
    //
    int num_threads;
    std::string name;

    Axis(){
      this->id_var = -1;
      this->order = -1;
      this->split_flag = 0;
      this->name = "";
      this->begin = 0;
      this->end = 0;
      this->stride = 1;
      this->num_threads = 0;
    }
    Axis(int id_var, int order, int split_flag, int begin, int end, int stride, int num_threads, std::string name): id_var(id_var), order(order), split_flag(split_flag), begin(begin), end(end), stride(stride), num_threads(num_threads), name(name) {}

    void set(int id_var, int order, int split_flag, int begin=0, int end=0, int stride=1, int num_threads=0, std::string name=""){
      this->id_var = id_var;
      this->order = order;
      this->split_flag = split_flag;
      this->begin = begin;
      this->end = end;
      this->stride = stride;
      this->num_threads = num_threads;
      this->name = name;
    }
    void set_num_threads(int num_threads){
      assert(num_threads>=0);
      this->num_threads = num_threads;
    }
    //operator<
    bool operator<(const Axis &a)const{
      return this->order < a.order;
      //if(this->id_var == a.id_var)return this->split_flag < a.split_flag;
      //return this->id_var < a.id_var;
    }
    void operator=(const Axis &a){
      this->id_var      = a.id_var;
      this->order       = a.order;
      this->split_flag  = a.split_flag;
      this->begin       = a.begin;
      this->end         = a.end;
      this->stride      = a.stride;
      this->num_threads = a.num_threads;
      this->name        = a.name;
    }
    bool operator==(const Axis &a)const{
      return 
      this->id_var      == a.id_var       &&
      this->order       == a.order        &&
      this->split_flag  == a.split_flag   &&
      this->begin       == a.begin        &&
      this->end         == a.end          &&
      this->stride      == a.stride       &&
      this->num_threads == a.num_threads  &&
      this->name        == a.name;
    }
    bool operator!=(const Axis &a)const{
      return !((*this) == a) ;
    }
    bool is_valid(){
      return this->id_var != -1 ;
    }
    void clear(){
      this->id_var = -1;
    }
    ~Axis(){
      this->id_var = -1;
      this->order = -1;
      this->split_flag = 0;
      this->name = "";
      this->end = 0;
      this->stride = 1;
      this->num_threads = 0;
    }
};

class CacheRead {
  public:
    std::shared_ptr<SpNode> SpNodePtr;
    Axis axis;
    std::string allocation;

    CacheRead(){
      this->SpNodePtr = nullptr;
    }
    CacheRead(std::shared_ptr<SpNode> SpNodePtr, std::string allocation="global"): SpNodePtr (SpNodePtr), allocation(allocation) {}

    void set(std::shared_ptr<SpNode> SpNodePtr){
      this->SpNodePtr = SpNodePtr;
    }
    void set(const Axis &x){
      this->axis = x;
    }
    void clear(){
      this->SpNodePtr = nullptr;
      this->axis.clear();
    }
    void operator=(const CacheRead &a){
      this->SpNodePtr = a.SpNodePtr;
      this->axis = a.axis;
    }
    bool operator==(const CacheRead &a)const{
      return (this->SpNodePtr == a.SpNodePtr && this->axis == a.axis);
    }
    bool use_cache_read(){
      return this->SpNodePtr != nullptr;
    }

    ~CacheRead(){
      this->SpNodePtr = nullptr;
    }
};

class CacheWrite {
  public:
    std::shared_ptr<TeNode> TeNodePtr;
    Axis axis;
    std::string allocation;

    CacheWrite(){
      this->TeNodePtr = nullptr;
    }
    CacheWrite(std::shared_ptr<TeNode> TeNodePtr, std::string allocation="global"): TeNodePtr (TeNodePtr), allocation(allocation) {}

    void set(std::shared_ptr<TeNode> TeNodePtr){
      this->TeNodePtr = TeNodePtr;
    }
    void set(const Axis &x){
      this->axis = x;
    }
    void clear(){
      this->TeNodePtr = nullptr;
      this->axis.clear();
    }
    void operator=(const CacheWrite &a){
      this->TeNodePtr = a.TeNodePtr;
      this->axis = a.axis;
    }
    bool operator==(const CacheWrite &a)const{
      return (this->TeNodePtr == a.TeNodePtr && this->axis == a.axis);
    }
    bool use_cache_write(){
      return this->TeNodePtr != nullptr;
    }

    ~CacheWrite(){
      this->TeNodePtr = nullptr;
    }
};

class ScheduleForStmt : public Stmt {
  public:
    ExprGroup eg;
    std::vector<Axis> axes;
    std::shared_ptr<Block> body;
    CacheRead  cacheread;
    CacheWrite cachewrite;

    ScheduleForStmt(){}
    ScheduleForStmt(ExprGroup eg, std::vector<Axis> axes, CacheRead  cacheread, CacheWrite cachewrite, std::shared_ptr<Block> body): eg(eg), axes(axes), cacheread(cacheread), cachewrite(cachewrite), body(body) {
    }

    ~ScheduleForStmt(){}

  DEFINE_ACCEPT
};

class AthreadStmt : public Stmt {
  public:
    int id_func_slave;
    std::shared_ptr<Expr> ArgTempPtr;
    std::shared_ptr<SpNode> SpNodePtr;
    std::shared_ptr<TeNode> TeNodePtr;

    AthreadStmt(){}
    AthreadStmt(int id_func_slave, std::shared_ptr<Expr> ArgTempPtr, std::shared_ptr<SpNode> SpNodePtr, std::shared_ptr<TeNode> TeNodePtr): id_func_slave(id_func_slave), ArgTempPtr(ArgTempPtr), SpNodePtr(SpNodePtr), TeNodePtr(TeNodePtr) {
    }

    ~AthreadStmt(){}

  DEFINE_ACCEPT
};

class FunctionStmt : public Stmt {
  public:
    static int counter;
    int id;
    //std::vector<std::shared_ptr<Expr> > args;
    std::vector<std::pair<DataType, std::shared_ptr<Expr>> > args;
    std::shared_ptr<Block> body;
    FunctionStmt(){}
    //FunctionStmt(std::vector<std::shared_ptr<Expr> > args, std::shared_ptr<Block> body): args(args), body(body) {
    FunctionStmt(std::vector<std::pair<DataType, std::shared_ptr<Expr>> > args, std::shared_ptr<Block> body): args(args), body(body) {
      this->id=this->counter++;
    }

    ~FunctionStmt(){}

  DEFINE_ACCEPT
};

int FunctionStmt::counter =0;

class PlaceHolderExpr : public Expr {
  public:
    std::shared_ptr<FrontendAllocaStmt> FrontendAllocaStmtPtr;
    std::shared_ptr<CallFuncExpr> CallFuncExprPtr;
    std::shared_ptr<TeNodeIndicesExpr> TeNodeIndicesExprPtr;

    PlaceHolderExpr(){}
    PlaceHolderExpr(std::shared_ptr<FrontendAllocaStmt> FrontendAllocaStmtPtr, std::shared_ptr<CallFuncExpr> CallFuncExprPtr, std::shared_ptr<TeNodeIndicesExpr> TeNodeIndicesExprPtr):FrontendAllocaStmtPtr(FrontendAllocaStmtPtr),CallFuncExprPtr(CallFuncExprPtr), TeNodeIndicesExprPtr(TeNodeIndicesExprPtr){
      set_datatype(TeNodeIndicesExprPtr->TeNodePtr->dt());
      set_shape(TeNodeIndicesExprPtr->TeNodePtr->shape());
    }

    ~PlaceHolderExpr(){}

  DEFINE_ACCEPT
};


#define DefVar(x, dt)  \
  auto __##x##__ = std::make_shared<VariableExpr>(DataType::dt); \
  ExprPtr x = std::dynamic_pointer_cast<Expr>(__##x##__); \
  current_ast_builder().insert(std::move(std::dynamic_pointer_cast<IRNode>(std::make_shared<FrontendAllocaStmt>(SHAPE_SCALAR , DataType::dt, __##x##__->id))));

#define DefVar_Value(x, dt, value)  \
  std::shared_ptr<VariableExpr> __##x##__; \
  if(DataType::dt == DataType::f32 ){   \
    float _value = ((float)value); \
    __##x##__ = std::make_shared<VariableExpr>(DataType::dt, _value); \
  }   \
  else if(DataType::dt == DataType::f64 ){   \
    double _value = ((double)value);    \
    __##x##__ = std::make_shared<VariableExpr>(DataType::dt, _value); \
  }   \
  else if(DataType::dt == DataType::i32 ){   \
    int _value = ((int)value);   \
    __##x##__ = std::make_shared<VariableExpr>(DataType::dt, _value); \
  }   \
  else{   \
    fmt::print("not correct type!");    \
    exit(-1);   \
  }   \
  ExprPtr x = std::dynamic_pointer_cast<Expr>(__##x##__); \
  current_ast_builder().insert(std::move(std::dynamic_pointer_cast<IRNode>(std::make_shared<FrontendAllocaStmt>(SHAPE_SCALAR , DataType::dt, __##x##__->id, __##x##__ ))));

#define DefTensor1D(x, dh, dt, d1)  \
  assert(dh >=0 && d1>=0 && d1>=dh*2);              \
  FrontendTensor x(VecInt {d1}, DataType::dt);  \
  (x.SpNodePtr)->set_halo_size( dh );            \
  (x.SpNodePtr)->set_time_window_size( -1 );            \
  current_ast_builder().insert(std::move(std::dynamic_pointer_cast<IRNode>(std::make_shared<FrontendAllocaStmt>(VecInt {d1+dh*2}, DataType::dt, (x.SpNodePtr)->id, (x.SpNodePtr) ))));

#define DefTensor2D(x, dh, dt, d1, d2)  \
  assert(dh >=0 && d1>=0 && d2>=0 && d1>=dh*2 && d2>=dh*2);              \
  FrontendTensor x(VecInt {d1, d2}, DataType::dt); \
  (x.SpNodePtr)->set_halo_size( dh );               \
  (x.SpNodePtr)->set_time_window_size( -1 );            \
  current_ast_builder().insert(std::move(std::dynamic_pointer_cast<IRNode>(std::make_shared<FrontendAllocaStmt>(VecInt {d1+dh*2, d2+dh*2}, DataType::dt, (x.SpNodePtr)->id, (x.SpNodePtr) ))));

#define DefTensor3D(x, dh, dt, d1, d2, d3)  \
  assert(dh >=0 && d1>=0 && d2>=0 && d3>=0 && d1>=dh*2 && d2>=dh*2 && d3>=dh*2);              \
  FrontendTensor x(VecInt {d1, d2, d3}, DataType::dt); \
  (x.SpNodePtr)->set_halo_size( dh );                    \
  (x.SpNodePtr)->set_time_window_size( -1 );            \
  current_ast_builder().insert(std::move(std::dynamic_pointer_cast<IRNode>(std::make_shared<FrontendAllocaStmt>(VecInt {d1+dh*2, d2+dh*2, d3+dh*2}, DataType::dt, (x.SpNodePtr)->id, (x.SpNodePtr) ))));

// DefTensorxD_TimeWin will optimize memory usage

#define DefTensor1D_TimeWin(x, win_size, dh, dt, d1)  \
  assert(win_size>0 && dh >=0 && d1>=0 && d1>=dh*2);              \
  FrontendTensor x(VecInt {d1}, DataType::dt);  \
  (x.SpNodePtr)->set_halo_size( dh );            \
  (x.SpNodePtr)->set_time_window_size( win_size );            \
  (x.SpNodePtr)->set_temp_size( win_size );            \
  current_ast_builder().insert(std::move(std::dynamic_pointer_cast<IRNode>(std::make_shared<FrontendAllocaStmt>(VecInt {d1+dh*2}, DataType::dt, (x.SpNodePtr)->id, (x.SpNodePtr) ))));

#define DefTensor2D_TimeWin(x, win_size, dh, dt, d1, d2)  \
  assert(win_size>0 && dh >=0 && d1>=0 && d2>=0 && d1>=dh*2 && d2>=dh*2);              \
  FrontendTensor x(VecInt {d1, d2}, DataType::dt); \
  (x.SpNodePtr)->set_halo_size( dh );               \
  (x.SpNodePtr)->set_time_window_size( win_size );            \
  (x.SpNodePtr)->set_temp_size( win_size );            \
  current_ast_builder().insert(std::move(std::dynamic_pointer_cast<IRNode>(std::make_shared<FrontendAllocaStmt>(VecInt {d1+dh*2, d2+dh*2}, DataType::dt, (x.SpNodePtr)->id, (x.SpNodePtr) ))));

#define DefTensor3D_TimeWin(x, win_size, dh, dt, d1, d2, d3)  \
  assert(win_size>0 && dh >=0 && d1>=0 && d2>=0 && d3>=0 && d1>=dh*2 && d2>=dh*2 && d3>=dh*2);              \
  FrontendTensor x(VecInt {d1, d2, d3}, DataType::dt); \
  (x.SpNodePtr)->set_halo_size( dh );                    \
  (x.SpNodePtr)->set_time_window_size( win_size );            \
  (x.SpNodePtr)->set_temp_size( win_size );            \
  current_ast_builder().insert(std::move(std::dynamic_pointer_cast<IRNode>(std::make_shared<FrontendAllocaStmt>(VecInt {d1+dh*2, d2+dh*2, d3+dh*2}, DataType::dt, (x.SpNodePtr)->id, (x.SpNodePtr) ))));

#define DefShapeMPI1D(x, D1) \
  assert(D1>0); \
  std::vector<int> x {D1};

#define DefShapeMPI2D(x, D1, D2) \
  assert(D1>0 && D2>0); \
  std::vector<int> x {D1, D2};

#define DefShapeMPI3D(x, D1, D2, D3) \
  assert(D1>0 && D2>0 && D3>0); \
  std::vector<int> x {D1, D2, D3};

const std::string schedule = "schedule";

#define Var ExprPtr

#define Tensor FrontendTensor

/*
ExprPtr operator+(ExprPtr lhs, ExprPtr rhs){
  auto _Bin = std::make_shared<BinaryOpExpr>(BinaryOpType::add, lhs.get(), rhs.get());
  ExprPtr Bin = std::dynamic_pointer_cast<Expr>(_Bin);
  return Bin;
}

ExprPtr operator*(ExprPtr lhs, ExprPtr rhs){
  auto _Bin = std::make_shared<BinaryOpExpr>(BinaryOpType::mul, lhs.get(), rhs.get());
  ExprPtr Bin = std::dynamic_pointer_cast<Expr>(_Bin);
  return Bin;
}
*/


#define DEFINE_CASE_OP_BINARY(op, opname)                                             \
  case BinaryOpType::opname:                                                          \
    fmt::print(#op);                                                                  \
    break;

#define DEFINE_CASE_OP_FUNC(opname)                                                   \
  case BinaryOpType::opname:                                                          \
    fmt::print(#opname);                                                              \
    break;

#define DEFINE_CASE_OP_UNARY(opname)                                                  \
  case UnaryOpType::opname:                                                           \
    fmt::print(#opname);                                                              \
    break;                                                                            

#define DEFINE_CASE_OP_BINARY_CODEGEN(op, opname)                                     \
  case BinaryOpType::opname:                                                          \
    add( #op );                                                                       \
    break;

#define DEFINE_CASE_OP_FUNC_CODEGEN(opname)                                           \
  case BinaryOpType::opname:                                                          \
    add( #opname );                                                                       \
    break;

#define DEFINE_CASE_OP_UNARY_CODEGEN(opname)                                          \
  case UnaryOpType::opname:                                                           \
    add( #opname );                                                                       \
    break;                                                                            


class TypeCheck : public IRVisitor {
  public:

    void visit(AssignOpExpr * op) override {
      op->lhs->accept(this);
      op->rhs->accept(this);
      VecInt shape_lhs = op->lhs->shape();
      VecInt shape_rhs = op->rhs->shape();
      if(same_vec(shape_lhs, shape_rhs)){
        op->set_shape(shape_lhs);
      }  
      else{
        std::cout<<"err: the shape of lhs is diff from rhs in AssignOp"<<std::endl;
        std::exit(-1);
      }
      //to check datatype
      DataType dt_lhs = op->lhs->dt();
      DataType dt_rhs = op->rhs->dt();
      if(is_uint(dt_lhs) || is_uint(dt_rhs) ){
        std::cout<<"err: DataType uint has not been supported in BinaryOp"<<std::endl;
        std::exit(-1);
      }
      DataType dt_up = convert_up(dt_lhs, dt_rhs);
      op->set_datatype(dt_up);
    }

    static void run(IRNode *node) {
      auto p = TypeCheck();
      fmt::print("==========\n");
      fmt::print("TypeCheck {{\n");
      node->accept(&p);
      fmt::print("}}\n");
      fmt::print("==========\n");
    }

    void visit(BinaryOpExpr *op) override {
      fmt::print("TypeCheck BinaryOp ");
      op->lhs->accept(this);
      op->rhs->accept(this);
      VecInt shape_lhs = op->lhs->shape();
      VecInt shape_rhs = op->rhs->shape();
      if(same_vec(shape_lhs, shape_rhs)){
        op->set_shape(shape_lhs);
      }  
      else if(same_vec(shape_lhs, SHAPE_SCALAR)){
        std::cout<<"lhs is a scalar in BinaryOp"<<std::endl;
        op->set_shape(shape_rhs);
      }
      else if(same_vec(shape_rhs, SHAPE_SCALAR)){
        std::cout<<"rhs is a scalar in BinaryOp"<<std::endl;
        op->set_shape(shape_lhs);
      }
      else{
        std::cout<<"err: the shape of lhs is diff from rhs and none of them are scalar in BinaryOp"<<std::endl;
        std::exit(-1);
      }
      //to check datatype
      DataType dt_lhs = op->lhs->dt();
      DataType dt_rhs = op->rhs->dt();
      if(is_uint(dt_lhs) || is_uint(dt_rhs) ){
        std::cout<<"err: DataType uint has not been supported in BinaryOp"<<std::endl;
        std::exit(-1);
      }
      DataType dt_up = convert_up(dt_lhs, dt_rhs);
      op->set_datatype(dt_up);
    }

    void visit(FuncOpExpr *op) override {
      fmt::print("TypeCheck FuncOp ");
      op->lhs->accept(this);
      op->rhs->accept(this);
      VecInt shape_lhs = op->lhs->shape();
      VecInt shape_rhs = op->rhs->shape();
      if(same_vec(shape_lhs, shape_rhs)){
        op->set_shape(shape_lhs);
      }  
      else if(same_vec(shape_lhs, SHAPE_SCALAR)){
        std::cout<<"lhs is a scalar in FuncOp"<<std::endl;
        op->set_shape(shape_rhs);
      }
      else if(same_vec(shape_rhs, SHAPE_SCALAR)){
        std::cout<<"rhs is a scalar in FuncOp"<<std::endl;
        op->set_shape(shape_lhs);
      }
      else{
        std::cout<<"err: the shape of lhs is diff from rhs and none of them are scalar in FuncOp"<<std::endl;
        std::exit(-1);
      }
      //to check datatype
      DataType dt_lhs = op->lhs->dt();
      DataType dt_rhs = op->rhs->dt();
      if(is_uint(dt_lhs) || is_uint(dt_rhs) ){
        std::cout<<"err: DataType uint has not been supported in FuncOp"<<std::endl;
        std::exit(-1);
      }
      DataType dt_up = convert_up(dt_lhs, dt_rhs);
      op->set_datatype(dt_up);
    }

    void visit(UnaryOpExpr *op) override{
      fmt::print("TypeCheck UnaryOp ");
      op->operand->accept(this);
      VecInt shape = op->operand->shape();
      DataType dt = op->operand->dt();
      op->set_shape(shape);
      op->set_datatype(dt);
    }

    void visit(VariableExpr *op) override {
      fmt::print("TypeCheck Variable ");
    }
    void visit(SpNodeIndicesExpr *op) override {
      fmt::print("TypeCheck SpNodeIndices ");
    }
};

class KernelIRVisitor : public IRVisitor {
  public:
    std::shared_ptr<Expr> ArgTempPtr;
    void set_arg_temp(std::shared_ptr<Expr> ArgTempPtr) {
      this->ArgTempPtr = ArgTempPtr;
    }
    static void run(IRNode *node, std::shared_ptr<Expr> ArgTempPtr) {
      auto p = KernelIRVisitor();
      p.set_arg_temp(ArgTempPtr);
      fmt::print("==========\n");
      fmt::print("KernelIRVisitor {{\n");
      node->accept(&p);
      fmt::print("}}\n");
      fmt::print("==========\n");
    }
    void visit(AssignOpExpr * op) override {
      op->lhs->accept(this);
      op->rhs->accept(this);
    }
    void visit(BinaryOpExpr *op) override {
      op->lhs->accept(this);
      op->rhs->accept(this);
    }
    void visit(FuncOpExpr *op) override {
      op->lhs->accept(this);
      op->rhs->accept(this);
    }
    void visit(UnaryOpExpr *op) override{
      op->operand->accept(this);
    }
    void visit(VariableExpr *op) override {
    }
    void visit(SpNodeIndicesExpr *op) override {
      op->set_arg_temp( this->ArgTempPtr ); 
    }
};

class VariableVisitor : public IRVisitor {
  public:
    static int id;
    static void set_id(int _id) {
      id = _id;
    }
    static int get_id(){
      return id;
    }
    static void run(IRNode *node) {
      auto p = VariableVisitor();
      node->accept(&p);
    }
    void visit(VariableExpr *op) override {
      set_id(op->id);
    }
};

int VariableVisitor::id ;

class GraphIRVisitorPlus : public IRVisitor {
  public:
    std::shared_ptr<Expr> ArgTempPtr;
    void set_arg_temp(std::shared_ptr<Expr> ArgTempPtr) {
      this->ArgTempPtr = ArgTempPtr;
    }
    static void run(IRNode *node, std::shared_ptr<Expr> ArgTempPtr) {
      auto p = GraphIRVisitorPlus();
      p.set_arg_temp(ArgTempPtr);
      fmt::print("==========\n");
      fmt::print("GraphIRVisitorPlus {{\n");
      node->accept(&p);
      fmt::print("}}\n");
      fmt::print("==========\n");
    }
    void visit(AssignOpExpr * op) override {
      op->lhs->accept(this);
      //op->rhs->accept(this);
    }
    void visit(SpNodeIndicesExpr *op) override {
      op->plus_arg_temp( this->ArgTempPtr );
      for(size_t i=0; i<(op->eg).size(); i++){
        ((op->eg)[i])->accept(this);
      }
    }
};

class GraphIRVisitor : public IRVisitor {
  public:
    static std::vector<std::shared_ptr<PlaceHolderExpr> > PlaceHolderExprPtrVecs;
    static void run(IRNode *node) {
      PlaceHolderExprPtrVecs.clear();
      auto p = GraphIRVisitor();
      fmt::print("==========\n");
      fmt::print("GraphIR {{\n");
      node->accept(&p);
      fmt::print("}}\n");
      fmt::print("==========\n");
    }
    void visit(BinaryOpExpr *op) override {
      op->lhs->accept(this);
      op->rhs->accept(this);
    }
    void visit(FuncOpExpr *op) override {
      op->lhs->accept(this);
      op->rhs->accept(this);
    }
    void visit(UnaryOpExpr *op) override{
      op->operand->accept(this);
    }
    void visit(VariableExpr *op) override {
    }
    void visit(SpNodeIndicesExpr *op) override {
      op->ArgTempPtr->accept(this);
      for(size_t i=0; i<(op->eg).size(); i++){
        ((op->eg)[i])->accept(this);
      }
    }
    void visit(TeNodeExpr * op) override {
    }
    void visit(AssignOpExpr * op) override {
      op->lhs->accept(this);
      op->rhs->accept(this);
    }

    void visit(PlaceHolderExpr * op) override {
      std::cout<<"from graphirvisitor in visiting placeholderexpr "<<std::endl;
      std::shared_ptr<PlaceHolderExpr> PlaceHolderExprPtr = 
              std::make_shared<PlaceHolderExpr>(op->FrontendAllocaStmtPtr, op->CallFuncExprPtr, op->TeNodeIndicesExprPtr);
      PlaceHolderExprPtrVecs.push_back(PlaceHolderExprPtr);
    }

    void visit(FrontendForStmt * op) override {
      op->body->accept(this);
    }

    void visit(RangeForStmt * op) override {
      op->body->accept(this);
    }

    void visit(Block * op) override {
      for(auto state: op->statements){
        state->accept(this);
      }
    }

    void visit(FunctionStmt * op) override {
      op->body->accept(this);
    }

};

std::vector<std::shared_ptr<PlaceHolderExpr> > GraphIRVisitor::PlaceHolderExprPtrVecs={};

class CodeGenCPU : public IRVisitor {
  public:
    static std::string Code_C;
    static void add(std::string input){
      Code_C += input;
    }
    static void init(){
      Code_C="";
    }
    static std::string run(IRNode *node, std::vector<int> vec_run_id , bool use_mpi=false, int timing_id=-1, int start=-1, int end=-1) {
      init();
      auto p = CodeGenCPU();
      add("\n//==========\n");
      add("//CodeGenCPU {{\n");
      add(headers_cpu);
      add(headers_omp);
      //mpi
      if(use_mpi){
        add(headers_mpi);
      }
      node->accept(&p);
      add("int main(){\n");
      //mpi
      if(use_mpi){
        add(" MPI_Init(NULL, NULL);\n");
        add(" MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);\n");
        add(" MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);\n");
      }
      for(auto v: vec_run_id){
        if(timing_id != -1 && timing_id == v){
          if(use_mpi){
            add(timing_start_mpi);
          }
          else{
            add(timing_start);
          }
        }
        add(std::string(" func")+std::to_string(v)+"();\n");
        if(timing_id != -1 && timing_id == v){
          add(std::string(" int timesteps = ")+std::to_string(end-start+1)+";\n");
          if(use_mpi){
            add(timing_end_mpi);
          }
          else{
            add(timing_end);
          }
        }
      }
      //mpi
      if(use_mpi){
        add(" MPI_Finalize();\n");
      }
      add(" return 0;\n}\n");
      add("//}}\n");
      add("//==========\n");
      return Code_C;
    }

    // print rules: one Expr as one atomic operation
    void visit(BinaryOpExpr *op) override {
      add(" (");
      op->lhs->accept(this);

      switch( op->type ){
        DEFINE_CASE_OP_BINARY_CODEGEN(+, add)
        DEFINE_CASE_OP_BINARY_CODEGEN(-, sub)
        DEFINE_CASE_OP_BINARY_CODEGEN(*, mul)
        DEFINE_CASE_OP_BINARY_CODEGEN(/, div)
        DEFINE_CASE_OP_BINARY_CODEGEN(%, mod)
        DEFINE_CASE_OP_BINARY_CODEGEN(&&, bit_and)
        DEFINE_CASE_OP_BINARY_CODEGEN(||, bit_or)
        DEFINE_CASE_OP_BINARY_CODEGEN(^, bit_xor)
        DEFINE_CASE_OP_BINARY_CODEGEN(<, cmp_lt)
        DEFINE_CASE_OP_BINARY_CODEGEN(<=, cmp_le)
        DEFINE_CASE_OP_BINARY_CODEGEN(>, cmp_gt)
        DEFINE_CASE_OP_BINARY_CODEGEN(>=, cmp_ge)
        DEFINE_CASE_OP_BINARY_CODEGEN(==, cmp_eq)
        DEFINE_CASE_OP_BINARY_CODEGEN(!=, cmp_ne)

        default:
          add("err : not correct Binary Op!");
          break;
      }

      op->rhs->accept(this);
      add(") ");

    }

    void visit(FuncOpExpr *op) override {
      switch( op->type ){
        DEFINE_CASE_OP_FUNC_CODEGEN(min)
        DEFINE_CASE_OP_FUNC_CODEGEN(max)
        DEFINE_CASE_OP_FUNC_CODEGEN(atan2)
        DEFINE_CASE_OP_FUNC_CODEGEN(truediv)
        DEFINE_CASE_OP_FUNC_CODEGEN(floordiv)

        default:
          add("err : not correct Func Op!");
          break;
      }

      add("(");
      op->lhs->accept(this);
      add(",");
      op->rhs->accept(this);
      add(") ");
    }

    void visit(UnaryOpExpr *op) override{
      switch( op->type ){

        DEFINE_CASE_OP_UNARY_CODEGEN(sqrt)
        DEFINE_CASE_OP_UNARY_CODEGEN(floor)
        DEFINE_CASE_OP_UNARY_CODEGEN(abs)
        DEFINE_CASE_OP_UNARY_CODEGEN(sin)
        DEFINE_CASE_OP_UNARY_CODEGEN(asin)
        DEFINE_CASE_OP_UNARY_CODEGEN(cos)
        DEFINE_CASE_OP_UNARY_CODEGEN(acos)
        DEFINE_CASE_OP_UNARY_CODEGEN(tan)
        DEFINE_CASE_OP_UNARY_CODEGEN(tanh)
        DEFINE_CASE_OP_UNARY_CODEGEN(inv)
        DEFINE_CASE_OP_UNARY_CODEGEN(rcp)
        DEFINE_CASE_OP_UNARY_CODEGEN(rsqrt)
        DEFINE_CASE_OP_UNARY_CODEGEN(exp)
        DEFINE_CASE_OP_UNARY_CODEGEN(log)

        default:
          add("err : not correct Unary Op!");
          break;
      }
      add("(");
      op->operand->accept(this);
      add(") ");

    }

    void visit(VariableExpr *op) override {
      //print("{}{} = {} {} {}", op->type_hint(), op->name(),
      //      binary_op_type_name(op->op_type), op->lhs->name(),
      //      op->rhs->name());
      add(" (");
      if(op->is_variable){
        add(std::string("var")+std::to_string(op->id));
      }
      else if(op->is_const_var){
        add(std::string("const_var")+std::to_string(op->id));
      }
      else if(op->is_immediate || op->is_assigned){
        if(op->dt() == DataType::f64)
          add(std::to_string(op->value_f64));
        else if(op->dt() == DataType::f32)
          add(std::to_string(op->value_f32));
        else if(op->dt() == DataType::i32)
          add(std::to_string(op->value_i32));
        else
          add("err : This DataType has not been implemented!");
      }
      else{ add("err : not found correct Variable!"); }
      add(") ");
    }

    void visit(SpNodeIndicesExpr *op) override {
      add(std::string(" SpNode")+std::to_string(op->SpNodePtr->id));
      if(op->SpNodePtr->time_window_size()>0){
        add("[ (");
        op->ArgTempPtr->accept(this);
        add(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ) ] ");
      }
      else if(op->SpNodePtr->time_window_size() < 0){
        add("[");
        op->ArgTempPtr->accept(this);
        add("] ");
      }
      for(size_t i=0; i<(op->eg).size(); i++){
        add("[");
        ((op->eg)[i])->accept(this);
        add(std::string(" + ")+std::to_string(op->SpNodePtr->halo_size()));
        add("] ");
        //if(i!=(op->eg).size()-1){
        //  add(",");
        //}
      }
    }

    void visit(MPIExpr *op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }
      //
      std::string type_call = "";
      if(op->type_call == 1)      type_call = "_start";
      else if(op->type_call == 2) type_call = "_end";
      else if(op->type_call == 3) type_call = "";
      //
      add(std::string(" exchange_halo_")+std::to_string(op->SpNodePtr->ndim())+"D_"+type+type_call+" ( ");
      add(std::string(" &SpNode")+std::to_string(op->SpNodePtr->id)+"[");
      for(auto v:op->args){
        v->accept(this);
      }
      // in specific time step
      if(op->SpNodePtr->use_time_window()){
        add(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ]");
      }
      else{
        add("]");
      }
      //
      for(int i=0;i<(op->SpNodePtr->shape()).size();i++){
        add("[0]");
      }
      //
      add(", ");
      add("0, ");
      for(int i=0;i<(op->SpNodePtr->shape()).size();i++){
        add(std::to_string((op->SpNodePtr->shape())[i])+", ");
      }
      add(std::to_string(op->SpNodePtr->halo_size())+", ");
      add("mpi_rank, ");
      for(int i=0;i<(op->shape_mpi).size();i++){
        add(std::to_string((op->shape_mpi)[i])+", ");
      }
      add("mpi_size ) ;\n");
    }

    void visit(OutputExpr *op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }
      //
      add(std::string(" OutputData_")+type+" ( ");
      add(std::string(" &SpNode")+std::to_string(op->SpNodePtr->id)+"[");
      for(auto v:op->args){
        v->accept(this);
      }
      // in specific time step
      if(op->SpNodePtr->use_time_window()){
        add(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ]");
      }
      else{
        add("]");
      }
      //
      for(int i=0;i<(op->shape()).size();i++){
        add("[0]");
      }
      add(", ");
      int size=1;
      for(auto v: op->shape()){
        size = size * (v + 2*(op->SpNodePtr->halo_size()));
      }
      add(std::to_string(size)+", ");
      add(std::string("\"")+op->name_file+"\" );\n");
    }

    void visit(InputExpr *op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }
      //
      add(std::string(" InputData_")+type+" ( ");
      add(std::string(" &SpNode")+std::to_string(op->SpNodePtr->id)+"[");
      add("0]");
      /*
      for(auto v:op->args){
        v->accept(this);
      }
      // in specific time step
      if(op->SpNodePtr->use_time_window()){
        add(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ]");
      }
      else{
        add("]");
      }
      */
      //
      for(int i=0;i<(op->shape()).size();i++){
        add("[0]");
      }
      add(", ");
      //
      int size = op->SpNodePtr->temp_size();
      if(op->SpNodePtr->use_time_window()){
        size = op->SpNodePtr->time_window_size();
      }
      for(auto v: op->shape()){
        size = size * (v + 2*(op->SpNodePtr->halo_size()));
      }
      add(std::to_string(size)+", ");
      add(std::string("\"")+op->name_file+"\" );\n");
    }

    void visit(TeNodeIndicesExpr *op) override {
      add(std::string(" TeNode")+std::to_string(op->TeNodePtr->id));
      for(size_t i=0; i<(op->eg).size(); i++){
        add("[");
        ((op->eg)[i])->accept(this);
        add("] ");
      }
    }

    void visit(TeNodeExpr * op) override {
      add(std::string(" TeNode")+std::to_string(op->TeNodePtr->id));
      for(size_t i=0; i<(op->shape()).size(); i++){
        add(std::string("[")+ std::to_string((op->shape())[i]) +"]");
      }
      std::cout<<" ";
    }

    void visit(PlaceHolderExpr * op) override {
      op->TeNodeIndicesExprPtr->accept(this);
    }
    
    void visit(AssignOpExpr * op) override {
      //add(" (");
      op->lhs->accept(this);
      add("=");
      op->rhs->accept(this);
      //add(") ");
      add(";\n");
    }

    void visit(CallFuncExpr * op) override {
      add(std::string(" func") + std::to_string(op->id_func) + "(");
      for(int i=0;i< (op->args).size();i++){
        ((op->args)[i])->accept(this);
        if(i!=(op->args).size()-1){
          add(",");   
        }
      }
      add(");\n");
    }

    void visit(CallNodeExpr * op) override {
      add(std::string(" TeNode")+std::to_string(op->TeNodePtr->id)+" ");
    }

    void visit(ScheduleForStmt * op) override {
      assert((op->axes).size()>0);
      /*
      for(auto v: op->axes){
        std::string split_flag="";
        if(v.split_flag==1)split_flag="_outer";
        if(v.split_flag==2)split_flag="_inner";
        add( std::string(" int var") + std::to_string(v.id_var)+ split_flag+" ;\n" );
      }
      */
      if( ((op->axes)[0]).num_threads > 0 ){
        std::string num_threads = std::to_string( ((op->axes)[0]).num_threads );
        add(std::string("# pragma omp parallel for num_threads (")+num_threads+")\n");
      }
      size_t num_loops = (op->axes).size();
      for(int i=0; i< num_loops; i++){
        auto v = (op->axes)[i];
        std::string split_flag="";
        if(v.split_flag==1)split_flag="_outer";
        if(v.split_flag==2)split_flag="_inner";
        //
        std::string name_var = std::string("var")+ std::to_string(v.id_var) + split_flag;
        add(std::string("for( int ")+ name_var +" = " +std::to_string(v.begin) + "; ");
        add(name_var + std::string(" < ")+ std::to_string(v.end)+"; ");
        add(name_var + std::string(" += ")+ std::to_string(v.stride)+" ) {\n");
      }
      //
      std::queue<Axis> que;
      for( auto v: op->axes ){
        que.push(v);
      }
      while(!que.empty()){
        Axis axis = que.front();
        que.pop();
        if(axis.split_flag==0)continue;
        bool found =false;
        std::queue<Axis> temp_que;
        Axis temp_axis;
        while(!que.empty()){
          temp_axis = que.front();
          que.pop();
          if(axis.id_var == temp_axis.id_var){
            found=true;
            break;
          }
          temp_que.push(temp_axis);
        }
        if(found==false){
          fmt::print("error : not found the corresponding outer and inner axis!");
          exit(-1);
        }
        while(!temp_que.empty()){
          que.push(temp_que.front());
          temp_que.pop();
        }
        if(axis.split_flag==2){
          std::swap(axis, temp_axis);
        }
        //
        std::string name_var = std::string("var")+ std::to_string(axis.id_var) ;
        add(std::string(" int ") + name_var + " = " + name_var+"_outer * "+std::to_string(temp_axis.end)+" + "+name_var+"_inner ;\n");
      }
      //
      op->body->accept(this);
      for(int i=0; i< num_loops; i++){
        add("}\n");
      }
    }

    void visit(FrontendForStmt * op) override {
      //
      size_t num_loops = (op->begins).size();
      //std::string num_threads = std::to_string(max_num_threads);
      if(op->num_threads>0){
        std::string num_threads = std::to_string(op->num_threads);
        add(std::string("# pragma omp parallel for num_threads (")+num_threads+")\n");
      }
      for(int i=0; i< num_loops; i++){
        add("for( int ");
        ((op->eg)[i])->accept(this);
        add(std::string("= ")+ std::to_string((op->begins)[i])+"; ");
        ((op->eg)[i])->accept(this);
        add(std::string("< ")+ std::to_string((op->ends)[i])+"; ");
        ((op->eg)[i])->accept(this);
        //add("+=1 ) {\n");
        add(std::string("+=")+ std::to_string((op->strides)[i])+" ) {\n");
      }
      op->body->accept(this);
      for(int i=0; i< num_loops; i++){
        add("}\n");
      }

      /*
      add("StructFor( ");
      for(int i=0;i<(op->begin_ends).size();i++){
        ((op->eg)[i])->accept(this);
        add(std::string(": 0 -> ") + std::to_string((op->begin_ends)[i]));
        if(i!=(op->begin_ends).size()-1){
          add(",");
        }
      }
      add(") {\n");
      op->body->accept(this);
      add("}\n");
      */
    }

    void visit(RangeForStmt * op) override {
      add("for( int ");
      op->loop_var->accept(this);
      add("= ");
      op->begin->accept(this);
      add("; ");
      op->loop_var->accept(this);
      add("<= ");
      op->end->accept(this);
      add("; ");
      op->loop_var->accept(this);
      add("+=1 ");
      add(") {\n");
      op->body->accept(this);
      add("}\n");
    }

    void visit(FrontendAllocaStmt * op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }

      if( same_vec(op->shape(), SHAPE_SCALAR) ){
        if( op->VariableExprPtr == nullptr ){
          add(type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->VariableExprPtr->is_assigned){
          std::string value ;
          if(type == "float")   value = std::to_string(op->VariableExprPtr->value_f32);
          if(type == "double")  value = std::to_string(op->VariableExprPtr->value_f64);
          if(type == "int")     value = std::to_string(op->VariableExprPtr->value_i32);
          add(type+" var"+std::to_string(op->id)+" = "+ value + " ;");
        }
        else{
          add("err : Variable is not assigned and try to been allocated!");
        }
        /*
        if(op->is_variable){
          add(type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->is_const_var){
          add(std::string("const ")+type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->is_immediate){
          add(std::to_string(op->value_f64));
        }
        else{ add("err : not found correct Variable!"); }
        */
      }
      else if(op->alloc_spnode()){
        add(type+" SpNode"+std::to_string(op->id));
        add(std::string("[")+ std::to_string(op->SpNodePtr->temp_size()) +"]");
        for(size_t i=0; i<(op->shape()).size(); i++){
          add(std::string("[")+ std::to_string((op->shape())[i]) +"]");
        }
        add(";");
      }
      else {//alloc TeNode
        //
        add(type+" TeNode"+std::to_string(op->id));
        for(size_t i=0; i<(op->shape()).size(); i++){
          add(std::string("[")+ std::to_string((op->shape())[i]) +"]");
        }
        add(";");
      }
      add("\n");
    }

    void visit(Block * op) override {
      for(auto state: op->statements){
        state->accept(this);
      }
    }

    void visit(FunctionStmt * op) override {
      add(std::string("void func")+ std::to_string(op->id)+"(");
      for(int i=0; i < (op->args).size(); i++){
        auto arg = (op->args)[i];
        std::string type;
        auto dt = arg.first;
        if(dt == DataType::f64){
          type="double";
        }
        else if(dt == DataType::f32){
          type="float";
        }
        else if(dt == DataType::i32){
          type="int";
        }
        else{
          add("err : This DataType has not been implemented!");
        }
        add(type);
        (arg.second)->accept(this);
        if(i!=(op->args).size()-1){
          add(", ");
        }
      }
      add(") {\n");
      op->body->accept(this);
      add("}\n");
    }

};

std::string CodeGenCPU::Code_C="";

class CodeGenSunwaySlave : public IRVisitor {
  public:
    static std::string Code_C;
    static void add(std::string input){
      Code_C += input;
    }
    static void init(){
      Code_C="";
    }
    static std::string run(IRNode *node) {
      init();
      auto p = CodeGenSunwaySlave();
      add("\n//==========\n");
      add("//CodeGenSunwaySlave {{\n");
      add(headers_slave);
      add(headers_argument_host_slave);
      add(headers_dma_get_put_slave); 
      //
      node->accept(&p);
      add("//}}\n");
      add("//==========\n");
      return Code_C;
    }

    // print rules: one Expr as one atomic operation
    void visit(BinaryOpExpr *op) override {
      add(" (");
      op->lhs->accept(this);

      switch( op->type ){
        DEFINE_CASE_OP_BINARY_CODEGEN(+, add)
        DEFINE_CASE_OP_BINARY_CODEGEN(-, sub)
        DEFINE_CASE_OP_BINARY_CODEGEN(*, mul)
        DEFINE_CASE_OP_BINARY_CODEGEN(/, div)
        DEFINE_CASE_OP_BINARY_CODEGEN(%, mod)
        DEFINE_CASE_OP_BINARY_CODEGEN(&&, bit_and)
        DEFINE_CASE_OP_BINARY_CODEGEN(||, bit_or)
        DEFINE_CASE_OP_BINARY_CODEGEN(^, bit_xor)
        DEFINE_CASE_OP_BINARY_CODEGEN(<, cmp_lt)
        DEFINE_CASE_OP_BINARY_CODEGEN(<=, cmp_le)
        DEFINE_CASE_OP_BINARY_CODEGEN(>, cmp_gt)
        DEFINE_CASE_OP_BINARY_CODEGEN(>=, cmp_ge)
        DEFINE_CASE_OP_BINARY_CODEGEN(==, cmp_eq)
        DEFINE_CASE_OP_BINARY_CODEGEN(!=, cmp_ne)

        default:
          add("err : not correct Binary Op!");
          break;
      }

      op->rhs->accept(this);
      add(") ");

    }

    void visit(FuncOpExpr *op) override {
      switch( op->type ){
        DEFINE_CASE_OP_FUNC_CODEGEN(min)
        DEFINE_CASE_OP_FUNC_CODEGEN(max)
        DEFINE_CASE_OP_FUNC_CODEGEN(atan2)
        DEFINE_CASE_OP_FUNC_CODEGEN(truediv)
        DEFINE_CASE_OP_FUNC_CODEGEN(floordiv)

        default:
          add("err : not correct Func Op!");
          break;
      }

      add("(");
      op->lhs->accept(this);
      add(",");
      op->rhs->accept(this);
      add(") ");
    }

    void visit(UnaryOpExpr *op) override{
      switch( op->type ){

        DEFINE_CASE_OP_UNARY_CODEGEN(sqrt)
        DEFINE_CASE_OP_UNARY_CODEGEN(floor)
        DEFINE_CASE_OP_UNARY_CODEGEN(abs)
        DEFINE_CASE_OP_UNARY_CODEGEN(sin)
        DEFINE_CASE_OP_UNARY_CODEGEN(asin)
        DEFINE_CASE_OP_UNARY_CODEGEN(cos)
        DEFINE_CASE_OP_UNARY_CODEGEN(acos)
        DEFINE_CASE_OP_UNARY_CODEGEN(tan)
        DEFINE_CASE_OP_UNARY_CODEGEN(tanh)
        DEFINE_CASE_OP_UNARY_CODEGEN(inv)
        DEFINE_CASE_OP_UNARY_CODEGEN(rcp)
        DEFINE_CASE_OP_UNARY_CODEGEN(rsqrt)
        DEFINE_CASE_OP_UNARY_CODEGEN(exp)
        DEFINE_CASE_OP_UNARY_CODEGEN(log)

        default:
          add("err : not correct Unary Op!");
          break;
      }
      add("(");
      op->operand->accept(this);
      add(") ");

    }

    void visit(VariableExpr *op) override {
      //print("{}{} = {} {} {}", op->type_hint(), op->name(),
      //      binary_op_type_name(op->op_type), op->lhs->name(),
      //      op->rhs->name());
      add(" (");
      if(op->is_variable){
        add(std::string("var")+std::to_string(op->id));
      }
      else if(op->is_const_var){
        add(std::string("const_var")+std::to_string(op->id));
      }
      else if(op->is_immediate || op->is_assigned){
        if(op->dt() == DataType::f64)
          add(std::to_string(op->value_f64));
        else if(op->dt() == DataType::f32)
          add(std::to_string(op->value_f32));
        else if(op->dt() == DataType::i32)
          add(std::to_string(op->value_i32));
        else
          add("err : This DataType has not been implemented!");
      }
      else{ add("err : not found correct Variable!"); }
      add(") ");
    }

    void visit(SpNodeIndicesExpr *op) override {
      add(std::string(" cache_read"));
      for(size_t i=0; i<(op->eg).size(); i++){
        add("[");
        ((op->eg)[i])->accept(this);
        add(std::string(" + ")+std::to_string(op->SpNodePtr->halo_size()));
        add("] ");
      }
    }

    void visit(TeNodeIndicesExpr *op) override {
      add(std::string(" cache_write"));
      for(size_t i=0; i<(op->eg).size(); i++){
        add("[");
        ((op->eg)[i])->accept(this);
        add("] ");
      }
    }

    void visit(TeNodeExpr * op) override {
      add(std::string(" TeNode")+std::to_string(op->TeNodePtr->id));
      for(size_t i=0; i<(op->shape()).size(); i++){
        add(std::string("[")+ std::to_string((op->shape())[i]) +"]");
      }
      std::cout<<" ";
    }

    void visit(PlaceHolderExpr * op) override {
      op->TeNodeIndicesExprPtr->accept(this);
    }
    
    void visit(AssignOpExpr * op) override {
      //add(" (");
      op->lhs->accept(this);
      add("=");
      op->rhs->accept(this);
      //add(") ");
      add(";\n");
    }

    void visit(CallFuncExpr * op) override {
      add(std::string(" func") + std::to_string(op->id_func) + "(");
      for(int i=0;i< (op->args).size();i++){
        ((op->args)[i])->accept(this);
        if(i!=(op->args).size()-1){
          add(",");   
        }
      }
      add(");\n");
    }

    void visit(CallNodeExpr * op) override {
      add(std::string(" TeNode")+std::to_string(op->TeNodePtr->id)+" ");
    }

    void visit(ScheduleForStmt * op) override {
      // assert((op->axes).size()>=2 && (op->axes)[0].id_var==(op->axes)[1].id_var && (op->axes)[0].split_flag==1 && (op->axes)[1].split_flag==2 && (op->axes)[0].order==0 && (op->axes)[1].order==1 );
      assert((op->cacheread).use_cache_read());
      std::string type;
      auto dt = ((op->cacheread).SpNodePtr)->dt();
      if(dt == DataType::f64){
        type="double";
      }
      else if(dt == DataType::f32){
        type="float";
      }
      else if(dt == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }

      /*
      for(auto v: op->axes){
        std::string split_flag="";
        if(v.split_flag==1)split_flag="_outer";
        if(v.split_flag==2)split_flag="_inner";
        add( std::string(" int var") + std::to_string(v.id_var)+ split_flag+" ;\n" );
      }
      */
      //ExprGroup eg;
      //std::vector<Axis> axes;
      //std::shared_ptr<Block> body;
      //CacheRead  cacheread;
      //CacheWrite cachewrite;

      //default 64 cores in CPE.
      add(" my_id = athread_get_id(-1);\n");
      //sorted axes
      //std::vector<Axis> sorted_axes = op->axes;
      //std::sort(sorted_axes.begin(), sorted_axes.end());
      for(auto v: op->axes){
        std::string split_flag="";
        if(v.split_flag==1)split_flag="_outer";
        if(v.split_flag==2)split_flag="_inner";
        std::string name_var = std::string("var")+ std::to_string(v.id_var) + split_flag;
        add(std::string(" int ")+name_var+" =0;\n");
      }
      //
      size_t num_loops = (op->axes).size();
      auto halo_size = ((op->cacheread).SpNodePtr)->halo_size();
      // auto tile_size = (op->axes)[1].end;
      /*
       *  2021/1/23 11:10 Modified by Bangduo Chen
       *  此处处理方案要求compute_at外层order必须为 {x,y,z}o, 即不能出现 {x,y,z}i
       *  且顺序必须为为xo, yo, zo,
       *  compute_at内层order不受限制
       *  readCache和writeCache应该位于同一循环内
       */
      // 找到要插入cacheRead和cacheWrite的循环位置, i的取值范围为问题维度范围, 此处受原设计限制搜索了划分后的所有loop
      int compute_at_index = -1;
      for (int i = 0; i < num_loops; i++) {
        auto v = (op->axes)[i];
        auto u = (op->cacheread).axis;
        if (v.split_flag == 2)
          break;
        // if ((op->axes)[i] == (op->cacheread).axis) {
        if (v.id_var == u.id_var && v.split_flag == u.split_flag) {
          compute_at_index = i;
          break;
        }
      }
      // NOTE: 此处compute_at_index不应该为-1, 否则是没找到相应的compute_at位置或compute_at外层的循环出现了inner
      assert(compute_at_index != -1 && "Find compute_at_index failed or inner loop appear out of compute_at position\n");

      /* 获取各维度tile_size, compute_at之外必须reorder后的顺序(xo, yo, zo), 获取reorder后顺序的的tile值*/      
      std::vector<int> tile_size;
      tile_size.resize(3);
      for (int axes_i = 0; axes_i < num_loops-1; axes_i++) {
        if ((op->axes)[axes_i].split_flag == 1) {
          // outer所在的循环, 搜索相应的inner循环
          for (int axes_j = axes_i+1; axes_j < num_loops; axes_j++) {
            if ((op->axes)[axes_i].id_var == (op->axes)[axes_j].id_var && (op->axes)[axes_j].split_flag == 2) {
              // tile_size.push_back((op->axes)[axes_j].end);
              tile_size[(op->axes)[axes_i].id_var] = (op->axes)[axes_j].end;
            }
          }
        }
      }

      // 计算并生成输入数组和cacheRead数组的大小字符串
      std::vector<int> shape_spnode = ((op->cacheread).SpNodePtr)->shape();
      std::vector<int> shape_tenode = ((op->cachewrite).TeNodePtr)->shape();
      std::string trans_shape_spnode = "";
      std::string index_shape_spnode = "";
      for(int i=0; i<shape_spnode.size(); i++){
        // if(i!=0)
        // trans_shape_spnode = trans_shape_spnode + "[" + std::to_string(shape_spnode[i]+2*halo_size) + "]";
        // //
        // int size = shape_spnode[i];
        // if(i==0)
        //   size = tile_size;
        // index_shape_spnode = index_shape_spnode + "[" + std::to_string(size+2*halo_size) + "]";
        int size = shape_spnode[i];
        // 输入数组的大小
        if (i != 0)
          trans_shape_spnode = trans_shape_spnode + "[" + std::to_string(shape_spnode[i]+2*halo_size) + "]";

        if (i <= compute_at_index)
          size = tile_size[i];
        // cacheRead部分大小
        index_shape_spnode = index_shape_spnode + "[" + std::to_string(size+2*halo_size) + "]";
      }

      // 计算并生成输出数组和cacheWrite数组的大小字符串
      std::string trans_shape_tenode = "";
      std::string index_shape_tenode = "";
      for(int i=0; i<shape_tenode.size(); i++){
        // if(i!=0)
        // trans_shape_tenode = trans_shape_tenode + "[" + std::to_string(shape_tenode[i]) + "]";
        // //
        // int size = shape_tenode[i];
        // if(i==0)
        //   size = tile_size;
        // index_shape_tenode = index_shape_tenode + "[" + std::to_string(size) + "]";
        int size = shape_tenode[i];
        // 输出数组大小
        if (i != 0)
          trans_shape_tenode = trans_shape_tenode + "[" + std::to_string(shape_tenode[i]) + "]";

        if (i <= compute_at_index)
          size = tile_size[i];   
        // cacheWrite部分大小
        index_shape_tenode = index_shape_tenode + "[" + std::to_string(size) + "]";
      }
      add(type+" cache_read "+index_shape_spnode+";\n");
      add(type+" cache_write"+index_shape_tenode+";\n");
      //add(type+" * SpNode = arg->SpNodeTemp;\n");
      //add(type+" * TeNode = arg->TeNode;\n");
      add(type+" (*SpNode)  "+trans_shape_spnode+ "= ( "+type+" (*) "+trans_shape_spnode+") arg->SpNodeTemp;\n");
      add(type+" (*TeNode) "+trans_shape_tenode+ "= ( "+type+" (*) "+trans_shape_tenode+") arg->TeNode;\n");
      
      // 创建嵌套循环for语句部分，并在指定的位置插入DMA get语句
      std::vector<std::string> name_var_list;
      for(int i=0; i< num_loops; i++){
        auto v = (op->axes)[i];
        std::string split_flag="";
        if(v.split_flag==1)split_flag="_outer";
        if(v.split_flag==2)split_flag="_inner";
        //
        std::string name_var = std::string("var")+ std::to_string(v.id_var) + split_flag;
        name_var_list.push_back(name_var);
        std::string begin_slave = std::to_string(v.begin);
        if(i==0) begin_slave = "my_id";
        add(std::string("for( ")+ name_var +" = " + begin_slave + "; ");
        add(name_var + std::string(" < ")+ std::to_string(v.end)+"; ");
        std::string stride_slave = std::to_string(v.stride);
        if(i==0) stride_slave = "64";
        add(name_var + std::string(" += ")+ stride_slave +" ) {\n");
        // if(i==0){
        //   add("//DMA_get\n");
        //   std::string index_spnode = "";
        //   std::string index_semi_spnode = "";
        //   int shape_size_spnode =1;
        //   for(int j=0; j<shape_spnode.size(); j++){
        //     auto v = shape_spnode[j];
        //     index_spnode = index_spnode +"[0]";
        //     if(j!=0)
        //     index_semi_spnode = index_semi_spnode +"[0]";
        //     if(j==0)
        //       shape_size_spnode = shape_size_spnode * (tile_size + 2 * halo_size);
        //     else
        //       shape_size_spnode = shape_size_spnode * (v + 2 * halo_size);
        //   }
        //   add("Get( &SpNode["+name_var+" * "+std::to_string((op->axes)[i+1].end)+" ]"+index_semi_spnode+", &cache_read"+index_spnode+ ", "+std::to_string(shape_size_spnode)+" * sizeof("+type+") );\n");
        // }
        if (i == compute_at_index) {
          add("//DMA_get\n");
          int dim = shape_spnode.size();
          // DMA 采用按面加载的方式， 一维被认为是1x一维宽度的面， 除三维外，其余维度均为load一个面， 三维则需要加载tile_z+2*halo_size个面
          int load_times = 1;
          if (dim == 3)
            load_times = tile_size[0] + 2*halo_size;
          add(std::string(" int load_times = ") + std::to_string(load_times)+";\n");
          std::string index_spnode = "";
          std::string index_semi_spnode = "";
          // 生成输入数组和cacheRead数组的维度索引
          int j;
          for (j = 0; j < dim; j++) {
            std::string z_iter_string = ((j == 0 && dim == 3) ? "+z_iter" : "");
            // 输入数组
            if (j <= compute_at_index)
              index_semi_spnode += ("[" + name_var_list[j] + "*" + std::to_string(tile_size[j]) + z_iter_string + "]");
            else 
              index_semi_spnode += "[0]";
            
            // cacheRead数组
            if (j == 0 && dim == 3) 
              index_spnode += "[z_iter]";
            else 
              index_spnode += "[0]";
          }
          // 计算Get语句参数
          int cnt = 1;
          int stride = 0;
          int bsize = 0;
          for (j = 0; j < dim; j++) {
            int dim_size = (j <= compute_at_index) ? (tile_size[j] + 2*halo_size) : (shape_spnode[j] + 2*halo_size);
            if (dim == 3 && j == 0) {
              // 三维情况按面加载, 不计算最高维度大小
              dim_size = 1;
            }
            cnt *= dim_size;
          }

          if (compute_at_index != 0 && compute_at_index == (dim-1)) {
            stride = shape_spnode[dim-1] - tile_size[dim-1];
            // bsize = (tile_size_x + 2*halo_size)
            bsize = tile_size[dim-1]+2*halo_size;
          }
          // 生成Get语句
          add("Get_Stride(&SpNode" + index_semi_spnode + ", &cache_read" 
              + index_spnode + ", " + "load_times" + ", " + std::to_string(cnt) + "*sizeof(" + type + ")"
              + ", " + std::to_string(stride) + "*sizeof(" + type + ")" + ", " + std::to_string(bsize) + "*sizeof(" + type + ")" + ")\n");
        }
      }
      //que for axes
      // std::queue<Axis> que;
      // for( auto v: op->axes ){
      //   que.push(v);
      // }
      // //the outer most axis and corresponding inner axis
      // Axis axis_outer = que.front();
      // que.pop();
      // Axis axis_inner = que.front();
      // que.pop();
      // assert(axis_outer.split_flag==1 && axis_inner.split_flag==2);
      // std::string name_var = std::string("var")+ std::to_string(axis_outer.id_var) ;
      // //add(std::string(" int ") + name_var + " = " + name_var+"_outer * "+std::to_string(axis_inner.end)+" + "+name_var+"_inner ;\n");
      // add(std::string(" int ") + name_var + " = " + name_var+"_inner ;\n");
      // //find the corresponding inner axis for each outer axis. 
      // while(!que.empty()){
      //   Axis axis = que.front();
      //   que.pop();
      //   if(axis.split_flag==0)continue;
      //   bool found =false;
      //   std::queue<Axis> temp_que;
      //   Axis temp_axis;
      //   while(!que.empty()){
      //     temp_axis = que.front();
      //     que.pop();
      //     if(axis.id_var == temp_axis.id_var){
      //       found=true;
      //       break;
      //     }
      //     temp_que.push(temp_axis);
      //   }
      //   if(found==false){
      //     fmt::print("error : not found the corresponding outer and inner axis!");
      //     exit(-1);
      //   }
      //   while(!temp_que.empty()){
      //     que.push(temp_que.front());
      //     temp_que.pop();
      //   }
      //   if(axis.split_flag==2){
      //     std::swap(axis, temp_axis);
      //   }
      //   //
      //   std::string name_var = std::string("var")+ std::to_string(axis.id_var) ;
      //   add(std::string(" int ") + name_var + " = " + name_var+"_outer * "+std::to_string(temp_axis.end)+" + "+name_var+"_inner ;\n");
      // }
      //

      // 遍历嵌套循环, 生成计算过程更新点坐标
      int tile_size_iter = 0;
      for (int axes_iter = 0; axes_iter < num_loops; axes_iter++) {
        // 跳过inner循环
        if ((op->axes)[axes_iter].split_flag != 1)
          continue;

        std::string name_var = std::string("var") + std::to_string((op->axes)[axes_iter].id_var);
        // 待处理的循环位于compute_at外侧, 相应维度坐标值直接设置为inner的值 
        if (axes_iter <= compute_at_index) {
          add(std::string(" int ") + name_var + " = " + name_var + "_inner;\n");
        } else {
          // 待处理的循环位于compute_at内侧, 相应维度的坐标值为outer*tile_size+inner
          add(std::string(" int ") + name_var + " = " + name_var+"_outer * "+std::to_string(tile_size[tile_size_iter])+" + "+name_var+"_inner ;\n");
        }
        tile_size_iter++;
      }
      
      op->body->accept(this);
      for(int i=0; i< num_loops; i++){
        // if(i==num_loops-1){
        //   add("//DMA_put\n");
        //   std::string index_tenode = "";
        //   std::string index_semi_tenode = "";
        //   int shape_size_tenode =1;
        //   for(int j=0; j<shape_tenode.size(); j++){
        //     auto v = shape_tenode[j];
        //     index_tenode = index_tenode +"[0]";
        //     if(j!=0)
        //     index_semi_tenode = index_semi_tenode +"[0]";
        //     if(j==0)
        //       shape_size_tenode *= tile_size ;
        //     else
        //       shape_size_tenode *= v ;
        //   }
        //   //name_var of outer most loop.
        //   int outer_most_loop = 0;
        //   auto v = (op->axes)[outer_most_loop];
        //   std::string split_flag="";
        //   if(v.split_flag==1)split_flag="_outer";
        //   if(v.split_flag==2)split_flag="_inner";
        //   //
        //   std::string name_var = std::string("var")+ std::to_string(v.id_var) + split_flag;
        //   //
        //   add(" Put( &TeNode["+name_var+" * "+std::to_string((op->axes)[outer_most_loop+1].end)+" ]"+index_semi_tenode+", &cache_write"+index_tenode+ ", "+std::to_string(shape_size_tenode)+" * sizeof("+type+") );\n");
        // }
        if (i == compute_at_index+1) {
          add("// DMA_put\n");
          int dim = shape_tenode.size();
          // DMA采用按面写回的方式，一维被认为是1x一维宽度的面，除三维外，其余维度均为store一个面，三维则需写回tile_z个面
          int store_times = 1;
          if (dim == 3)
            store_times = tile_size[0];
          add(std::string(" int store_times = ") + std::to_string(store_times)+";\n");

          std::string index_tenode = "";
          std::string index_semi_tenode = "";
          // 生成输出数组和cacheWrite数组的维度索引
          int j;
          for (j = 0; j < dim; j++) {
            std::string z_iter_string = ((j == 0 && dim == 3) ? "+z_iter" : "");
            // 输出数组
            if (j <= compute_at_index)
              index_semi_tenode += ("[" + name_var_list[j] + "*" + std::to_string(tile_size[j]) + z_iter_string + "]");
            else
              index_semi_tenode += "[0]";
            
            // cacheWrite数组
            if (j == 0 && dim == 3)
              index_tenode += "[z_iter]";
            else 
              index_tenode += "[0]";
          }
          // 计算put语句参数
          int cnt = 1;
          int stride = 0;
          int bsize = 0;
          for (j = 0; j < dim; j++) {
            int dim_size = (j <= compute_at_index) ? tile_size[j] : shape_tenode[j];
            if (dim == 3 && j == 0) {
              // 三维情况按面加载, 不计算最高维度大小
              dim_size = 1;
            }
            cnt *= dim_size;
          }

          if (compute_at_index != 0 && compute_at_index == (dim - 1)) {
            stride = shape_tenode[dim-1] - tile_size[dim-1];
            // bsize = tile_size_x
            bsize = tile_size[dim-1];
          }

          // 生成Put语句
          add("Put_Stride(&TeNode"+index_semi_tenode + ", &cache_write"
              + index_tenode + ", " + "store_times" + ", " + std::to_string(cnt) + "*sizeof(" + type + ")"
              + ", " + std::to_string(stride) + "*sizeof(" + type + ")" + ", " + std::to_string(bsize) + "*sizeof(" + type + ")" + ") \n");
        }
        add("}\n");
      }
    }

    void visit(FrontendForStmt * op) override {
      //
      size_t num_loops = (op->begins).size();
      //std::string num_threads = std::to_string(max_num_threads);
      if(op->num_threads>0){
        std::string num_threads = std::to_string(op->num_threads);
        add(std::string("# pragma omp parallel for num_threads (")+num_threads+")\n");
      }
      for(int i=0; i< num_loops; i++){
        add("for( int ");
        ((op->eg)[i])->accept(this);
        add(std::string("= ")+ std::to_string((op->begins)[i])+"; ");
        ((op->eg)[i])->accept(this);
        add(std::string("< ")+ std::to_string((op->ends)[i])+"; ");
        ((op->eg)[i])->accept(this);
        //add("+=1 ) {\n");
        add(std::string("+=")+ std::to_string((op->strides)[i])+" ) {\n");
      }
      op->body->accept(this);
      for(int i=0; i< num_loops; i++){
        add("}\n");
      }
    }

    void visit(RangeForStmt * op) override {
      add("for( int ");
      op->loop_var->accept(this);
      add("= ");
      op->begin->accept(this);
      add("; ");
      op->loop_var->accept(this);
      add("<= ");
      op->end->accept(this);
      add("; ");
      op->loop_var->accept(this);
      add("+=1 ");
      add(") {\n");
      op->body->accept(this);
      add("}\n");
    }

    void visit(FrontendAllocaStmt * op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }

      if( same_vec(op->shape(), SHAPE_SCALAR) ){
        if( op->VariableExprPtr == nullptr ){
          add(type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->VariableExprPtr->is_assigned){
          std::string value ;
          if(type == "float")   value = std::to_string(op->VariableExprPtr->value_f32);
          if(type == "double")  value = std::to_string(op->VariableExprPtr->value_f64);
          if(type == "int")     value = std::to_string(op->VariableExprPtr->value_i32);
          add(type+" var"+std::to_string(op->id)+" = "+ value + " ;");
        }
        else{
          add("err : Variable is not assigned and try to been allocated!");
        }
        /*
        if(op->is_variable){
          add(type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->is_const_var){
          add(std::string("const ")+type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->is_immediate){
          add(std::to_string(op->value_f64));
        }
        else{ add("err : not found correct Variable!"); }
        */
      }
      else if(op->alloc_spnode()){
        add(type+" SpNode"+std::to_string(op->id));
        add(std::string("[")+ std::to_string(op->SpNodePtr->temp_size()) +"]");
        for(size_t i=0; i<(op->shape()).size(); i++){
          add(std::string("[")+ std::to_string((op->shape())[i]) +"]");
        }
        add(";");
      }
      else {//alloc TeNode
        //
        add(type+" TeNode"+std::to_string(op->id));
        for(size_t i=0; i<(op->shape()).size(); i++){
          add(std::string("[")+ std::to_string((op->shape())[i]) +"]");
        }
        add(";");
      }
      add("\n");
    }

    void visit(Block * op) override {
      for(auto state: op->statements){
        state->accept(this);
      }
    }

    void visit(FunctionStmt * op) override {
      add(std::string("void func")+ std::to_string(op->id)+"(");
      for(int i=0; i < (op->args).size(); i++){
        auto arg = (op->args)[i];
        std::string type;
        auto dt = arg.first;
        if(dt == DataType::f64){
          type="double";
        }
        else if(dt == DataType::f32){
          type="float";
        }
        else if(dt == DataType::i32){
          type="int";
        }
        else{
          add("err : This DataType has not been implemented!");
        }
        add(std::string(" Arg_")+type+" * arg ");
        //(arg.second)->accept(this);
        //if(i!=(op->args).size()-1){
        //  add(", ");
        //}
      }
      add(") {\n");
      op->body->accept(this);
      add("}\n");
    }

};

std::string CodeGenSunwaySlave::Code_C="";


class CodeGenSunway : public IRVisitor {
  public:
    static std::string Code_C;
    static void add(std::string input){
      Code_C += input;
    }
    static void init(){
      Code_C="";
    }
    static std::string run(IRNode *node, std::vector<int> vec_run_id , bool use_mpi=false, int timing_id =-1, int start=-1, int end=-1) {
      init();
      auto p = CodeGenSunway();
      add("\n//==========\n");
      add("//CodeGenSunway {{\n");
      add(headers_athread);
      add(headers_cpu);
      add(headers_argument_host_slave);
      //mpi
      if(use_mpi){
        add(headers_mpi);
      }
      //extern SLAVE_FUN
      for(auto v: call_extern_func_slave){
        auto dt = v.second;
        std::string type;
        if(dt == DataType::f64){
          type="double";
        }
        else if(dt == DataType::f32){
          type="float";
        }
        else if(dt == DataType::i32){
          type="int";
        }
        else{
          add("err : This DataType has not been implemented!");
        }
        add(std::string("extern SLAVE_FUN(func")+std::to_string(v.first)+")(Arg_"+type+" * arg);\n");
      }
      add("\n");
      //
      node->accept(&p);
      add("int main(){\n");
      //mpi
      if(use_mpi){
        add(" MPI_Init(NULL, NULL);\n");
        add(" MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);\n");
        add(" MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);\n");
      }
      add(" athread_init();\n");
      for(auto v: vec_run_id){
        if(timing_id != -1 && timing_id == v){
          if(use_mpi){
            add(timing_start_mpi);
          }
          else{
            add(timing_start);
          }

        }
        add(std::string(" func")+std::to_string(v)+"();\n");
        if(timing_id != -1 && timing_id == v){
          add(std::string(" int timesteps = ")+std::to_string(end-start+1)+";\n");
          if(use_mpi){
            add(timing_end_mpi);
          }
          else{
            add(timing_end);
          }

        }
      }
      add(" athread_halt();\n");
      //mpi
      if(use_mpi){
        add(" MPI_Finalize();\n");
      }
      add(" return 0;\n}\n");
      add("//}}\n");
      add("//==========\n");
      return Code_C;
    }

    // print rules: one Expr as one atomic operation
    void visit(BinaryOpExpr *op) override {
      add(" (");
      op->lhs->accept(this);

      switch( op->type ){
        DEFINE_CASE_OP_BINARY_CODEGEN(+, add)
        DEFINE_CASE_OP_BINARY_CODEGEN(-, sub)
        DEFINE_CASE_OP_BINARY_CODEGEN(*, mul)
        DEFINE_CASE_OP_BINARY_CODEGEN(/, div)
        DEFINE_CASE_OP_BINARY_CODEGEN(%, mod)
        DEFINE_CASE_OP_BINARY_CODEGEN(&&, bit_and)
        DEFINE_CASE_OP_BINARY_CODEGEN(||, bit_or)
        DEFINE_CASE_OP_BINARY_CODEGEN(^, bit_xor)
        DEFINE_CASE_OP_BINARY_CODEGEN(<, cmp_lt)
        DEFINE_CASE_OP_BINARY_CODEGEN(<=, cmp_le)
        DEFINE_CASE_OP_BINARY_CODEGEN(>, cmp_gt)
        DEFINE_CASE_OP_BINARY_CODEGEN(>=, cmp_ge)
        DEFINE_CASE_OP_BINARY_CODEGEN(==, cmp_eq)
        DEFINE_CASE_OP_BINARY_CODEGEN(!=, cmp_ne)

        default:
          add("err : not correct Binary Op!");
          break;
      }

      op->rhs->accept(this);
      add(") ");

    }

    void visit(FuncOpExpr *op) override {
      switch( op->type ){
        DEFINE_CASE_OP_FUNC_CODEGEN(min)
        DEFINE_CASE_OP_FUNC_CODEGEN(max)
        DEFINE_CASE_OP_FUNC_CODEGEN(atan2)
        DEFINE_CASE_OP_FUNC_CODEGEN(truediv)
        DEFINE_CASE_OP_FUNC_CODEGEN(floordiv)

        default:
          add("err : not correct Func Op!");
          break;
      }

      add("(");
      op->lhs->accept(this);
      add(",");
      op->rhs->accept(this);
      add(") ");
    }

    void visit(UnaryOpExpr *op) override{
      switch( op->type ){

        DEFINE_CASE_OP_UNARY_CODEGEN(sqrt)
        DEFINE_CASE_OP_UNARY_CODEGEN(floor)
        DEFINE_CASE_OP_UNARY_CODEGEN(abs)
        DEFINE_CASE_OP_UNARY_CODEGEN(sin)
        DEFINE_CASE_OP_UNARY_CODEGEN(asin)
        DEFINE_CASE_OP_UNARY_CODEGEN(cos)
        DEFINE_CASE_OP_UNARY_CODEGEN(acos)
        DEFINE_CASE_OP_UNARY_CODEGEN(tan)
        DEFINE_CASE_OP_UNARY_CODEGEN(tanh)
        DEFINE_CASE_OP_UNARY_CODEGEN(inv)
        DEFINE_CASE_OP_UNARY_CODEGEN(rcp)
        DEFINE_CASE_OP_UNARY_CODEGEN(rsqrt)
        DEFINE_CASE_OP_UNARY_CODEGEN(exp)
        DEFINE_CASE_OP_UNARY_CODEGEN(log)

        default:
          add("err : not correct Unary Op!");
          break;
      }
      add("(");
      op->operand->accept(this);
      add(") ");

    }

    void visit(VariableExpr *op) override {
      //print("{}{} = {} {} {}", op->type_hint(), op->name(),
      //      binary_op_type_name(op->op_type), op->lhs->name(),
      //      op->rhs->name());
      add(" (");
      if(op->is_variable){
        add(std::string("var")+std::to_string(op->id));
      }
      else if(op->is_const_var){
        add(std::string("const_var")+std::to_string(op->id));
      }
      else if(op->is_immediate || op->is_assigned){
        if(op->dt() == DataType::f64)
          add(std::to_string(op->value_f64));
        else if(op->dt() == DataType::f32)
          add(std::to_string(op->value_f32));
        else if(op->dt() == DataType::i32)
          add(std::to_string(op->value_i32));
        else
          add("err : This DataType has not been implemented!");
      }
      else{ add("err : not found correct Variable!"); }
      add(") ");
    }

    void visit(SpNodeIndicesExpr *op) override {
      add(std::string(" SpNode")+std::to_string(op->SpNodePtr->id));
      if(op->SpNodePtr->time_window_size()>0){
        add("[ (");
        op->ArgTempPtr->accept(this);
        add(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ) ] ");
      }
      else if(op->SpNodePtr->time_window_size() < 0){
        add("[");
        op->ArgTempPtr->accept(this);
        add("] ");
      }
      for(size_t i=0; i<(op->eg).size(); i++){
        add("[");
        ((op->eg)[i])->accept(this);
        add(std::string(" + ")+std::to_string(op->SpNodePtr->halo_size()));
        add("] ");
        //if(i!=(op->eg).size()-1){
        //  add(",");
        //}
      }
    }

    void visit(MPIExpr *op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }
      //
      std::string type_call = "";
      if(op->type_call == 1)      type_call = "_start";
      else if(op->type_call == 2) type_call = "_end";
      else if(op->type_call == 3) type_call = "";
      //
      add(std::string(" exchange_halo_")+std::to_string(op->SpNodePtr->ndim())+"D_"+type+type_call+" ( ");
      add(std::string(" &SpNode")+std::to_string(op->SpNodePtr->id)+"[");
      for(auto v:op->args){
        v->accept(this);
      }
      // in specific time step
      if(op->SpNodePtr->use_time_window()){
        add(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ]");
      }
      else{
        add("]");
      }
      //
      for(int i=0;i<(op->SpNodePtr->shape()).size();i++){
        add("[0]");
      }
      //
      add(", ");
      add("0, ");
      for(int i=0;i<(op->SpNodePtr->shape()).size();i++){
        add(std::to_string((op->SpNodePtr->shape())[i])+", ");
      }
      add(std::to_string(op->SpNodePtr->halo_size())+", ");
      add("mpi_rank, ");
      for(int i=0;i<(op->shape_mpi).size();i++){
        add(std::to_string((op->shape_mpi)[i])+", ");
      }
      add("mpi_size ) ;\n");
    }

    void visit(OutputExpr *op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }
      //
      add(std::string(" OutputData_")+type+" ( ");
      add(std::string(" &SpNode")+std::to_string(op->SpNodePtr->id)+"[");
      for(auto v:op->args){
        v->accept(this);
      }
      // in specific time step
      if(op->SpNodePtr->use_time_window()){
        add(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ]");
      }
      else{
        add("]");
      }
      //
      for(int i=0;i<(op->shape()).size();i++){
        add("[0]");
      }
      add(", ");
      int size=1;
      for(auto v: op->shape()){
        size = size * (v + 2*(op->SpNodePtr->halo_size()));
      }
      add(std::to_string(size)+", ");
      add(std::string("\"")+op->name_file+"\" );\n");
    }

    void visit(InputExpr *op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }
      //
      add(std::string(" InputData_")+type+" ( ");
      add(std::string(" &SpNode")+std::to_string(op->SpNodePtr->id)+"[");
      add("0]");
      /*
      for(auto v:op->args){
        v->accept(this);
      }
      // in specific time step
      if(op->SpNodePtr->use_time_window()){
        add(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ]");
      }
      else{
        add("]");
      }
      */
      //
      for(int i=0;i<(op->shape()).size();i++){
        add("[0]");
      }
      add(", ");
      //
      int size = op->SpNodePtr->temp_size();
      if(op->SpNodePtr->use_time_window()){
        size = op->SpNodePtr->time_window_size();
      }
      for(auto v: op->shape()){
        size = size * (v + 2*(op->SpNodePtr->halo_size()));
      }
      add(std::to_string(size)+", ");
      add(std::string("\"")+op->name_file+"\" );\n");
    }

    void visit(TeNodeIndicesExpr *op) override {
      add(std::string(" TeNode")+std::to_string(op->TeNodePtr->id));
      for(size_t i=0; i<(op->eg).size(); i++){
        add("[");
        ((op->eg)[i])->accept(this);
        add("] ");
      }
    }

    void visit(TeNodeExpr * op) override {
      add(std::string(" TeNode")+std::to_string(op->TeNodePtr->id));
      for(size_t i=0; i<(op->shape()).size(); i++){
        add(std::string("[")+ std::to_string((op->shape())[i]) +"]");
      }
      std::cout<<" ";
    }

    void visit(PlaceHolderExpr * op) override {
      op->TeNodeIndicesExprPtr->accept(this);
    }
    
    void visit(AssignOpExpr * op) override {
      //add(" (");
      op->lhs->accept(this);
      add("=");
      op->rhs->accept(this);
      //add(") ");
      add(";\n");
    }

    void visit(CallFuncExpr * op) override {
      add(std::string(" func") + std::to_string(op->id_func) + "(");
      for(int i=0;i< (op->args).size();i++){
        ((op->args)[i])->accept(this);
        if(i!=(op->args).size()-1){
          add(",");   
        }
      }
      add(");\n");
    }

    void visit(CallNodeExpr * op) override {
      add(std::string(" TeNode")+std::to_string(op->TeNodePtr->id)+" ");
    }

    void visit(AthreadStmt * op) override {
      //int id_func_slave;
      //std::shared_ptr<SpNode> SpNodePtr;
      //std::shared_ptr<TeNode> TeNodePtr;
      DataType dt = op->SpNodePtr->dt();
      std::string type;
      if(dt == DataType::f64){
        type="double";
      }
      else if(dt == DataType::f32){
        type="float";
      }
      else if(dt == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }
      std::vector<int> shape_spnode = op->SpNodePtr->shape();
      std::vector<int> shape_tenode = op->TeNodePtr->shape();
      //TeNode
      add(std::string("Arg_")+type+" arg = { &TeNode"+std::to_string(op->TeNodePtr->id));
      for(int i=0; i<shape_tenode.size(); i++){
        add("[0]");
      }
      //SpNode
      add(std::string(", &SpNode")+std::to_string(op->SpNodePtr->id));
      if(op->SpNodePtr->time_window_size()>0){
        add("[ (");
        op->ArgTempPtr->accept(this);
        add(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ) ] ");
      }
      else if(op->SpNodePtr->time_window_size() < 0){
        add("[");
        op->ArgTempPtr->accept(this);
        add("] ");
      }
      for(size_t i=0; i<shape_spnode.size(); i++){
        add("[0]");
      }
      add(" };\n");
      //athread
      add(std::string("athread_spawn( func")+std::to_string(op->id_func_slave)+", &arg );\n");
      add("athread_join();\n");
    }

    void visit(ScheduleForStmt * op) override {
      assert((op->axes).size()>0);
      /*
      for(auto v: op->axes){
        std::string split_flag="";
        if(v.split_flag==1)split_flag="_outer";
        if(v.split_flag==2)split_flag="_inner";
        add( std::string(" int var") + std::to_string(v.id_var)+ split_flag+" ;\n" );
      }
      */
      if( ((op->axes)[0]).num_threads > 0 ){
        std::string num_threads = std::to_string( ((op->axes)[0]).num_threads );
        add(std::string("# pragma omp parallel for num_threads (")+num_threads+")\n");
      }
      size_t num_loops = (op->axes).size();
      for(int i=0; i< num_loops; i++){
        auto v = (op->axes)[i];
        std::string split_flag="";
        if(v.split_flag==1)split_flag="_outer";
        if(v.split_flag==2)split_flag="_inner";
        //
        std::string name_var = std::string("var")+ std::to_string(v.id_var) + split_flag;
        add(std::string("for( int ")+ name_var +" = " +std::to_string(v.begin) + "; ");
        add(name_var + std::string(" < ")+ std::to_string(v.end)+"; ");
        add(name_var + std::string(" += ")+ std::to_string(v.stride)+" ) {\n");
      }
      //
      std::queue<Axis> que;
      for( auto v: op->axes ){
        que.push(v);
      }
      while(!que.empty()){
        Axis axis = que.front();
        que.pop();
        if(axis.split_flag==0)continue;
        bool found =false;
        std::queue<Axis> temp_que;
        Axis temp_axis;
        while(!que.empty()){
          temp_axis = que.front();
          que.pop();
          if(axis.id_var == temp_axis.id_var){
            found=true;
            break;
          }
          temp_que.push(temp_axis);
        }
        if(found==false){
          fmt::print("error : not found the corresponding outer and inner axis!");
          exit(-1);
        }
        while(!temp_que.empty()){
          que.push(temp_que.front());
          temp_que.pop();
        }
        if(axis.split_flag==2){
          std::swap(axis, temp_axis);
        }
        //
        std::string name_var = std::string("var")+ std::to_string(axis.id_var) ;
        add(std::string(" int ") + name_var + " = " + name_var+"_outer * "+std::to_string(temp_axis.end)+" + "+name_var+"_inner ;\n");
      }
      //
      op->body->accept(this);
      for(int i=0; i< num_loops; i++){
        add("}\n");
      }
    }

    void visit(FrontendForStmt * op) override {
      //
      size_t num_loops = (op->begins).size();
      //std::string num_threads = std::to_string(max_num_threads);
      if(op->num_threads>0){
        std::string num_threads = std::to_string(op->num_threads);
        add(std::string("# pragma omp parallel for num_threads (")+num_threads+")\n");
      }
      for(int i=0; i< num_loops; i++){
        add("for( int ");
        ((op->eg)[i])->accept(this);
        add(std::string("= ")+ std::to_string((op->begins)[i])+"; ");
        ((op->eg)[i])->accept(this);
        add(std::string("< ")+ std::to_string((op->ends)[i])+"; ");
        ((op->eg)[i])->accept(this);
        //add("+=1 ) {\n");
        add(std::string("+=")+ std::to_string((op->strides)[i])+" ) {\n");
      }
      op->body->accept(this);
      for(int i=0; i< num_loops; i++){
        add("}\n");
      }

      /*
      add("StructFor( ");
      for(int i=0;i<(op->begin_ends).size();i++){
        ((op->eg)[i])->accept(this);
        add(std::string(": 0 -> ") + std::to_string((op->begin_ends)[i]));
        if(i!=(op->begin_ends).size()-1){
          add(",");
        }
      }
      add(") {\n");
      op->body->accept(this);
      add("}\n");
      */
    }

    void visit(RangeForStmt * op) override {
      add("for( int ");
      op->loop_var->accept(this);
      add("= ");
      op->begin->accept(this);
      add("; ");
      op->loop_var->accept(this);
      add("<= ");
      op->end->accept(this);
      add("; ");
      op->loop_var->accept(this);
      add("+=1 ");
      add(") {\n");
      op->body->accept(this);
      add("}\n");
    }

    void visit(FrontendAllocaStmt * op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        add("err : This DataType has not been implemented!");
      }

      if( same_vec(op->shape(), SHAPE_SCALAR) ){
        if( op->VariableExprPtr == nullptr ){
          add(type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->VariableExprPtr->is_assigned){
          std::string value ;
          if(type == "float")   value = std::to_string(op->VariableExprPtr->value_f32);
          if(type == "double")  value = std::to_string(op->VariableExprPtr->value_f64);
          if(type == "int")     value = std::to_string(op->VariableExprPtr->value_i32);
          add(type+" var"+std::to_string(op->id)+" = "+ value + " ;");
        }
        else{
          add("err : Variable is not assigned and try to been allocated!");
        }
        /*
        if(op->is_variable){
          add(type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->is_const_var){
          add(std::string("const ")+type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->is_immediate){
          add(std::to_string(op->value_f64));
        }
        else{ add("err : not found correct Variable!"); }
        */
      }
      else if(op->alloc_spnode()){
        add(type+" SpNode"+std::to_string(op->id));
        add(std::string("[")+ std::to_string(op->SpNodePtr->temp_size()) +"]");
        for(size_t i=0; i<(op->shape()).size(); i++){
          add(std::string("[")+ std::to_string((op->shape())[i]) +"]");
        }
        add(";");
      }
      else {//alloc TeNode
        //
        add(type+" TeNode"+std::to_string(op->id));
        for(size_t i=0; i<(op->shape()).size(); i++){
          add(std::string("[")+ std::to_string((op->shape())[i]) +"]");
        }
        add(";");
      }
      add("\n");
    }

    void visit(Block * op) override {
      for(auto state: op->statements){
        state->accept(this);
      }
    }

    void visit(FunctionStmt * op) override {
      add(std::string("void func")+ std::to_string(op->id)+"(");
      for(int i=0; i < (op->args).size(); i++){
        auto arg = (op->args)[i];
        std::string type;
        auto dt = arg.first;
        if(dt == DataType::f64){
          type="double";
        }
        else if(dt == DataType::f32){
          type="float";
        }
        else if(dt == DataType::i32){
          type="int";
        }
        else{
          add("err : This DataType has not been implemented!");
        }
        add(type);
        (arg.second)->accept(this);
        if(i!=(op->args).size()-1){
          add(", ");
        }
      }
      add(") {\n");
      op->body->accept(this);
      add("}\n");
    }

};

std::string CodeGenSunway::Code_C="";



class IRPrinter : public IRVisitor {
  public:
    static void run(IRNode *node) {
      auto p = IRPrinter();
      fmt::print("==========\n");
      fmt::print("IR {{\n");
      node->accept(&p);
      fmt::print("}}\n");
      fmt::print("==========\n");
    }

    // print rules: one Expr as one atomic operation
    void visit(BinaryOpExpr *op) override {
      fmt::print(" (");
      //fmt::print("(");
      op->lhs->accept(this);
      //fmt::print(")");

      switch( op->type ){
        DEFINE_CASE_OP_BINARY(+, add)
        DEFINE_CASE_OP_BINARY(-, sub)
        DEFINE_CASE_OP_BINARY(*, mul)
        DEFINE_CASE_OP_BINARY(/, div)
        DEFINE_CASE_OP_BINARY(%, mod)
        DEFINE_CASE_OP_BINARY(&&, bit_and)
        DEFINE_CASE_OP_BINARY(||, bit_or)
        DEFINE_CASE_OP_BINARY (^, bit_xor)
        DEFINE_CASE_OP_BINARY(<, cmp_lt)
        DEFINE_CASE_OP_BINARY(<=, cmp_le)
        DEFINE_CASE_OP_BINARY(>, cmp_gt)
        DEFINE_CASE_OP_BINARY(>=, cmp_ge)
        DEFINE_CASE_OP_BINARY(==, cmp_eq)
        DEFINE_CASE_OP_BINARY(!=, cmp_ne)

        default:
          fmt::print("err : not correct Binary Op!");
          break;
      }

      //fmt::print("(");
      op->rhs->accept(this);
      //fmt::print(")");
      fmt::print(") ");

    }

    void visit(FuncOpExpr *op) override {
      switch( op->type ){
        DEFINE_CASE_OP_FUNC(min)
        DEFINE_CASE_OP_FUNC(max)
        DEFINE_CASE_OP_FUNC(atan2)
        DEFINE_CASE_OP_FUNC(truediv)
        DEFINE_CASE_OP_FUNC(floordiv)

        default:
          fmt::print("err : not correct Func Op!");
          break;
      }

      fmt::print("(");
      op->lhs->accept(this);
      fmt::print(",");
      op->rhs->accept(this);
      fmt::print(") ");
    }

    void visit(UnaryOpExpr *op) override{
      switch( op->type ){

        DEFINE_CASE_OP_UNARY(sqrt)
        DEFINE_CASE_OP_UNARY(floor)
        DEFINE_CASE_OP_UNARY(abs)
        DEFINE_CASE_OP_UNARY(sin)
        DEFINE_CASE_OP_UNARY(asin)
        DEFINE_CASE_OP_UNARY(cos)
        DEFINE_CASE_OP_UNARY(acos)
        DEFINE_CASE_OP_UNARY(tan)
        DEFINE_CASE_OP_UNARY(tanh)
        DEFINE_CASE_OP_UNARY(inv)
        DEFINE_CASE_OP_UNARY(rcp)
        DEFINE_CASE_OP_UNARY(rsqrt)
        DEFINE_CASE_OP_UNARY(exp)
        DEFINE_CASE_OP_UNARY(log)

        default:
          fmt::print("err : not correct Unary Op!");
          break;
      }
      fmt::print("(");
      op->operand->accept(this);
      fmt::print(") ");

    }

    void visit(VariableExpr *op) override {
      //print("{}{} = {} {} {}", op->type_hint(), op->name(),
      //      binary_op_type_name(op->op_type), op->lhs->name(),
      //      op->rhs->name());
      fmt::print(" (");
      if(op->is_variable){
        fmt::print(std::string("var")+std::to_string(op->id));
      }
      else if(op->is_const_var){
        fmt::print(std::string("const_var")+std::to_string(op->id));
      }
      else if(op->is_immediate || op->is_assigned){
        if(op->dt() == DataType::f64)
          fmt::print(std::to_string(op->value_f64));
        else if(op->dt() == DataType::f32)
          fmt::print(std::to_string(op->value_f32));
        else if(op->dt() == DataType::i32)
          fmt::print(std::to_string(op->value_i32));
        else
          fmt::print("err : This DataType has not been implemented!");
      }
      else{ fmt::print("err : not found correct Variable!"); }
      fmt::print(") ");
    }

    void visit(SpNodeIndicesExpr *op) override {
      fmt::print(std::string(" SpNode")+std::to_string(op->SpNodePtr->id));
      if(op->SpNodePtr->time_window_size()>0){
        fmt::print("[ (");
        op->ArgTempPtr->accept(this);
        fmt::print(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ) ] ");
      }
      else if(op->SpNodePtr->time_window_size() < 0){
        fmt::print("[");
        op->ArgTempPtr->accept(this);
        fmt::print("] ");
      }
      for(size_t i=0; i<(op->eg).size(); i++){
        fmt::print("[");
        ((op->eg)[i])->accept(this);
        fmt::print(std::string(" + ")+std::to_string(op->SpNodePtr->halo_size()));
        fmt::print("] ");
        //if(i!=(op->eg).size()-1){
        //  fmt::print(",");
        //}
      }
    }

    void visit(MPIExpr *op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        fmt::print("err : This DataType has not been implemented!");
      }
      //
      std::string type_call = "";
      if(op->type_call == 1)      type_call = "_start";
      else if(op->type_call == 2) type_call = "_end";
      else if(op->type_call == 3) type_call = "";
      //
      fmt::print(std::string(" exchange_halo_")+std::to_string(op->SpNodePtr->ndim())+"D_"+type+type_call+" ( ");
      fmt::print(std::string(" &SpNode")+std::to_string(op->SpNodePtr->id)+"[");
      for(auto v:op->args){
        v->accept(this);
      }
      // in specific time step
      if(op->SpNodePtr->use_time_window()){
        fmt::print(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ]");
      }
      else{
        fmt::print("]");
      }
      //
      for(int i=0;i<(op->SpNodePtr->shape()).size();i++){
        fmt::print("[0]");
      }
      //
      fmt::print(", ");
      fmt::print("0, ");
      for(int i=0;i<(op->SpNodePtr->shape()).size();i++){
        fmt::print(std::to_string((op->SpNodePtr->shape())[i])+", ");
      }
      fmt::print(std::to_string(op->SpNodePtr->halo_size())+", ");
      fmt::print("mpi_rank, ");
      for(int i=0;i<(op->shape_mpi).size();i++){
        fmt::print(std::to_string((op->shape_mpi)[i])+", ");
      }
      fmt::print("mpi_size ) ;\n");
    }

    void visit(OutputExpr *op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        fmt::print("err : This DataType has not been implemented!");
      }
      //
      fmt::print(std::string(" OutputData_")+type+" ( ");
      fmt::print(std::string(" &SpNode")+std::to_string(op->SpNodePtr->id)+"[");
      for(auto v:op->args){
        v->accept(this);
      }
      // in specific time step
      if(op->SpNodePtr->use_time_window()){
        fmt::print(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ]");
      }
      else{
        fmt::print("]");
      }
      //
      for(int i=0;i<(op->shape()).size();i++){
        fmt::print("[0]");
      }
      fmt::print(", ");
      int size=1;
      for(auto v: op->shape()){
        size = size * (v + 2*(op->SpNodePtr->halo_size()));
      }
      fmt::print(std::to_string(size)+", ");
      fmt::print(std::string("\"")+op->name_file+"\" );\n");
    }

    void visit(InputExpr *op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        fmt::print("err : This DataType has not been implemented!");
      }
      //
      fmt::print(std::string(" InputData_")+type+" ( ");
      fmt::print(std::string(" &SpNode")+std::to_string(op->SpNodePtr->id)+"[");
      fmt::print("0]");
      /*
      for(auto v:op->args){
        v->accept(this);
      }
      // in specific time step
      if(op->SpNodePtr->use_time_window()){
        fmt::print(std::string("% ")+std::to_string(op->SpNodePtr->time_window_size())+" ]");
      }
      else{
        fmt::print("]");
      }
      */
      //
      for(int i=0;i<(op->shape()).size();i++){
        fmt::print("[0]");
      }
      fmt::print(", ");
      //
      int size = op->SpNodePtr->temp_size();
      if(op->SpNodePtr->use_time_window()){
        size = op->SpNodePtr->time_window_size();
      }
      for(auto v: op->shape()){
        size = size * (v + 2*(op->SpNodePtr->halo_size()));
      }
      fmt::print(std::to_string(size)+", ");
      fmt::print(std::string("\"")+op->name_file+"\" );\n");
    }

    void visit(TeNodeIndicesExpr *op) override {
      fmt::print(std::string(" TeNode")+std::to_string(op->TeNodePtr->id));
      for(size_t i=0; i<(op->eg).size(); i++){
        fmt::print("[");
        ((op->eg)[i])->accept(this);
        fmt::print("] ");
      }
    }

    void visit(TeNodeExpr * op) override {
      fmt::print(std::string(" TeNode")+std::to_string(op->TeNodePtr->id));
      for(size_t i=0; i<(op->shape()).size(); i++){
        fmt::print(std::string("[")+ std::to_string((op->shape())[i]) +"]");
      }
      std::cout<<" ";
    }

    void visit(PlaceHolderExpr * op) override {
      op->TeNodeIndicesExprPtr->accept(this);
    }
    
    void visit(AssignOpExpr * op) override {
      //fmt::print(" (");
      op->lhs->accept(this);
      fmt::print("=");
      op->rhs->accept(this);
      //fmt::print(") ");
      fmt::print(";\n");
    }

    void visit(CallFuncExpr * op) override {
      fmt::print(std::string(" func") + std::to_string(op->id_func) + "(");
      for(int i=0;i< (op->args).size();i++){
        ((op->args)[i])->accept(this);
        if(i!=(op->args).size()-1){
          fmt::print(",");   
        }
      }
      fmt::print(");\n");
    }

    void visit(CallNodeExpr * op) override {
      fmt::print(std::string(" TeNode")+std::to_string(op->TeNodePtr->id)+" ");
    }

    void visit(ScheduleForStmt * op) override {
      assert((op->axes).size()>0);
      /*
      for(auto v: op->axes){
        std::string split_flag="";
        if(v.split_flag==1)split_flag="_outer";
        if(v.split_flag==2)split_flag="_inner";
        fmt::print( std::string(" int var") + std::to_string(v.id_var)+ split_flag+" ;\n" );
      }
      */
      if( ((op->axes)[0]).num_threads > 0 ){
        std::string num_threads = std::to_string( ((op->axes)[0]).num_threads );
        fmt::print(std::string("# pragma omp parallel for num_threads (")+num_threads+")\n");
      }
      size_t num_loops = (op->axes).size();
      for(int i=0; i< num_loops; i++){
        auto v = (op->axes)[i];
        std::string split_flag="";
        if(v.split_flag==1)split_flag="_outer";
        if(v.split_flag==2)split_flag="_inner";
        //
        std::string name_var = std::string("var")+ std::to_string(v.id_var) + split_flag;
        fmt::print(std::string("for( int ")+ name_var +" = " +std::to_string(v.begin) + "; ");
        fmt::print(name_var + std::string(" < ")+ std::to_string(v.end)+"; ");
        fmt::print(name_var + std::string(" += ")+ std::to_string(v.stride)+" ) {\n");
      }
      //
      std::queue<Axis> que;
      for( auto v: op->axes ){
        que.push(v);
      }
      while(!que.empty()){
        Axis axis = que.front();
        que.pop();
        if(axis.split_flag==0)continue;
        bool found =false;
        std::queue<Axis> temp_que;
        Axis temp_axis;
        while(!que.empty()){
          temp_axis = que.front();
          que.pop();
          if(axis.id_var == temp_axis.id_var){
            found=true;
            break;
          }
          temp_que.push(temp_axis);
        }
        if(found==false){
          fmt::print("error : not found the corresponding outer and inner axis!");
          exit(-1);
        }
        while(!temp_que.empty()){
          que.push(temp_que.front());
          temp_que.pop();
        }
        if(axis.split_flag==2){
          std::swap(axis, temp_axis);
        }
        //
        std::string name_var = std::string("var")+ std::to_string(axis.id_var) ;
        fmt::print(std::string(" int ") + name_var + " = " + name_var+"_outer * "+std::to_string(temp_axis.end)+" + "+name_var+"_inner ;\n");
      }
      //
      op->body->accept(this);
      for(int i=0; i< num_loops; i++){
        fmt::print("}\n");
      }
    }

    void visit(FrontendForStmt * op) override {
      if(op->num_threads>0){
        std::string num_threads = std::to_string(op->num_threads);
        fmt::print(std::string("# pragma omp parallel for num_threads (")+num_threads+")\n");
      }
      fmt::print("StructFor( ");
      for(int i=0;i<(op->begins).size();i++){
        ((op->eg)[i])->accept(this);
        fmt::print(std::string(": ") + std::to_string(op->begins[i]) + " -> " + std::to_string((op->ends)[i]) +" (stride,"+ std::to_string((op->strides)[i])+") ");
        if(i!=(op->begins).size()-1){
          fmt::print(",");
        }
      }
      fmt::print(") {\n");
      op->body->accept(this);
      /*
      for(auto state: op->body->statements){
        state->accept(this);
      }
      */
      fmt::print("}\n");
    }

    void visit(RangeForStmt * op) override {
      fmt::print("RangeFor( ");
      op->loop_var->accept(this);
      fmt::print(": ");
      op->begin->accept(this);
      fmt::print(" -> ");
      op->end->accept(this);
      fmt::print(") {\n");
      op->body->accept(this);
      fmt::print("}\n");
    }

    void visit(FrontendAllocaStmt * op) override {
      std::string type;
      if(op->dt() == DataType::f64){
        type="double";
      }
      else if(op->dt() == DataType::f32){
        type="float";
      }
      else if(op->dt() == DataType::i32){
        type="int";
      }
      else{
        fmt::print("err : This DataType has not been implemented!");
      }

      if( same_vec(op->shape(), SHAPE_SCALAR) ){
        if( op->VariableExprPtr == nullptr ){
          fmt::print(type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->VariableExprPtr->is_assigned){
          std::string value ;
          if(type == "float")   value = std::to_string(op->VariableExprPtr->value_f32);
          if(type == "double")  value = std::to_string(op->VariableExprPtr->value_f64);
          if(type == "int")     value = std::to_string(op->VariableExprPtr->value_i32);
          fmt::print(type+" var"+std::to_string(op->id)+" = "+ value + " ;");
        }
        else{
          fmt::print("err : Variable is not assigned and try to been allocated!");
        }
        /*
        if(op->is_variable){
          fmt::print(type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->is_const_var){
          fmt::print(std::string("const ")+type+" var"+std::to_string(op->id)+std::string(";"));
        }
        else if(op->is_immediate){
          fmt::print(std::to_string(op->value_f64));
        }
        else{ fmt::print("err : not found correct Variable!"); }
        */
      }
      else if(op->alloc_spnode()){
        fmt::print(type+" SpNode"+std::to_string(op->id));
        fmt::print(std::string("[")+ std::to_string(op->SpNodePtr->temp_size()) +"]");
        for(size_t i=0; i<(op->shape()).size(); i++){
          fmt::print(std::string("[")+ std::to_string((op->shape())[i]) +"]");
        }
        fmt::print(";");
      }
      else {//alloc TeNode
        //
        fmt::print(type+" TeNode"+std::to_string(op->id));
        for(size_t i=0; i<(op->shape()).size(); i++){
          fmt::print(std::string("[")+ std::to_string((op->shape())[i]) +"]");
        }
        fmt::print(";");
      }
      fmt::print("\n");
    }

    void visit(Block * op) override {
      for(auto state: op->statements){
        state->accept(this);
      }
    }

    void visit(FunctionStmt * op) override {
      fmt::print(std::string("func")+ std::to_string(op->id)+"(");
      for(int i=0; i < (op->args).size(); i++){
        auto arg = (op->args)[i];
        std::string type;
        auto dt = arg.first;
        if(dt == DataType::f64){
          type="double";
        }
        else if(dt == DataType::f32){
          type="float";
        }
        else if(dt == DataType::i32){
          type="int";
        }
        else{
          fmt::print("err : This DataType has not been implemented!");
        }
        fmt::print(type);
        (arg.second)->accept(this);
        if(i!=(op->args).size()-1){
          fmt::print(", ");
        }
      }
      fmt::print(") {\n");
      op->body->accept(this);
      fmt::print("}\n");
    }

};

void print(IRNode * root){
  IRPrinter::run(root);
}

void PrintIR(){
  IRPrinter::run( context->root() );      
}

class Result {
  public:
    ExprPtr _root;
    ExprGroup _eg;

    Result(){}
    Result(ExprGroup eg, ExprPtr root){
      this->_root=std::move(root);
      this->_eg = eg;
    }
    ExprPtr operator[](int temporal){
      return this->temporal(temporal);
    }
    ExprPtr temporal(int temporal){
      auto arg_temp = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(temporal));
      KernelIRVisitor::run( (this->_root).get(), arg_temp );
      return this->_root;
    }
    ~Result(){}
};
 
class Kernel {
  public:
    ExprPtr _root;
    ExprGroup _eg;
    //the result of Kernel computation
    std::shared_ptr<Expr> TeNodeExprPtr;
    std::shared_ptr<TeNode> TeNodePtr;
    std::shared_ptr<SpNode> SpNodePtr;
    std::shared_ptr<FunctionStmt> FunctionStmtPtr;
    std::vector<Axis> axes;
    CacheRead  cacheread;
    CacheWrite cachewrite;


    Kernel(){}
    Kernel( ExprPtr root ){
      this->_root=std::move(root);
    }
    Kernel(ExprGroup eg, ExprPtr root, std::string is_schedule ){
      assert(is_schedule == schedule);
      //TODO
      //pass transform shape of spatial to temporal
      TypeCheck::run( root.get() );
      VecInt shape = root->shape();
      DataType dt = root->dt();
      this->SpNodePtr = std::make_shared<SpNode>(shape.size(), shape, dt);
      //
      auto TeNodePtr = std::make_shared<TeNode>(0, shape.size(), shape, dt);
      this->TeNodePtr = TeNodePtr;
      std::shared_ptr<Expr> TeNodeExprPtr = std::dynamic_pointer_cast<Expr>(std::make_shared<TeNodeExpr>(TeNodePtr));
      std::shared_ptr<Expr> TeNodeIndicesExprPtr = std::dynamic_pointer_cast<Expr>(std::make_shared<TeNodeIndicesExpr>(TeNodePtr, eg));
      //
      auto kernelexpr = TeNodeIndicesExprPtr << root;
      //
      this->_root=std::move(kernelexpr);
      this->_eg = eg;
      //TypeCheck for two times
      TypeCheck::run( (this->_root).get() );

      this->TeNodeExprPtr = TeNodeExprPtr;
      //axes
      this->axes.clear();
      int ndim = shape.size();
      std::vector<int> ids_eg;
      for(int i=0; i<eg.size(); i++){
        VariableVisitor::run( eg[i].get() );
        ids_eg.push_back( VariableVisitor::id );
      }
      for(int i=0; i<ndim; i++){
        Axis axis (ids_eg[i], i, 0, 0, shape[i], 1, 0, "") ;
        //std::cout<<"id_eg "<<ids_eg[i]<<" "<<i<<" "<<shape[i]<<std::endl;
        axes.push_back(axis);
      }
    }
    //1D
    void tile(int tile_size_x, Axis &xo, Axis &xi ){
      assert( this->axes.size() == 1 );
      std::vector<Axis> _axes = this->axes;
      this->axes.clear();
      const int ndim = _axes.size();
      for(auto v: _axes){
        assert(v.end % tile_size_x==0);
        assert(v.split_flag==0);
        Axis axis_outer(v.id_var, v.order*ndim,   1, 0, v.end/tile_size_x,  1, 0, "");
        Axis axis_inner(v.id_var, v.order*ndim+1, 2, 0, tile_size_x,        1, 0, "");
        this->axes.push_back(axis_outer);
        this->axes.push_back(axis_inner);
      }
      xo=this->axes[0];
      xi=this->axes[1];
    }
    //2D
    void tile(int tile_size_x, int tile_size_y, Axis &xo, Axis &xi, Axis &yo, Axis &yi ){
      printf("%d\n", this->axes.size());
      assert( this->axes.size() == 2 );
      std::vector<Axis> _axes = this->axes;
      this->axes.clear();
      //the original id of variable and for loop, which is not splited.
      //the order of nested for loop with splited or not.
      //0:not splited, 1:splited and outer, 2:splited and inner
      const int ndim = _axes.size();
      std::vector<int> vec_tile_size {tile_size_x, tile_size_y};
      for(int i=0; i< ndim; i++){
        auto v = _axes[i];
        assert(v.end % vec_tile_size[i]==0);
        assert(v.split_flag==0);
        Axis axis_outer(v.id_var, v.order*ndim,   1, 0, v.end/vec_tile_size[i],  1, 0, "");
        Axis axis_inner(v.id_var, v.order*ndim+1, 2, 0, vec_tile_size[i],        1, 0, "");
        this->axes.push_back(axis_outer);
        this->axes.push_back(axis_inner);
      }
      xo=this->axes[0];
      xi=this->axes[1];
      yo=this->axes[2];
      yi=this->axes[3];
    }
    //3D
    void tile(int tile_size_x, int tile_size_y, int tile_size_z, Axis &xo, Axis &xi, Axis &yo, Axis &yi, Axis &zo, Axis &zi  ){
      assert( this->axes.size() == 3 );
      std::vector<Axis> _axes = this->axes;
      this->axes.clear();
      const int ndim = _axes.size();
      std::vector<int> vec_tile_size {tile_size_x, tile_size_y, tile_size_z};
      for(int i=0; i< ndim; i++){
        auto v = _axes[i];
        assert(v.end % vec_tile_size[i]==0);
        assert(v.split_flag==0);
        Axis axis_outer(v.id_var, v.order*ndim,   1, 0, v.end/vec_tile_size[i],  1, 0, "");
        Axis axis_inner(v.id_var, v.order*ndim+1, 2, 0, vec_tile_size[i],        1, 0, "");
        this->axes.push_back(axis_outer);
        this->axes.push_back(axis_inner);
      }
      xo=this->axes[0];
      xi=this->axes[1];
      yo=this->axes[2];
      yi=this->axes[3];
      zo=this->axes[4];
      zi=this->axes[5];
    }

    bool _same_axes( std::vector<Axis> a, std::vector<Axis> b ){
      assert(a.size() == b.size());
      std::sort(a.begin(), a.end());
      std::sort(b.begin(), b.end());
      for(int i=0; i<a.size(); i++){
        if(a[i] != b[i])return false;
      }
      return true;
    }
    //1D
    void reorder( Axis &x0, Axis &x1 ){
      assert( this->axes.size() == 2 );
      std::vector<Axis> reordered_axes {x0, x1};
      if( !this->_same_axes(this->axes, reordered_axes) ){
        fmt::print("the input axis is not same as axes in kernel when reorder!");
        exit(-1);
      }
      this->axes[0] = x0;
      this->axes[1] = x1;
    }
    //2D
    void reorder( Axis &x0, Axis &x1, Axis &x2, Axis &x3 ){
      assert( this->axes.size() == 4 );
      std::vector<Axis> reordered_axes {x0, x1, x2, x3};
      if( !this->_same_axes(this->axes, reordered_axes) ){
        fmt::print("the input axis is not same as axes in kernel when reorder!");
        exit(-1);
      }
      this->axes[0] = x0;
      this->axes[1] = x1;
      this->axes[2] = x2;
      this->axes[3] = x3;
    }
    //3D
    void reorder( Axis &x0, Axis &x1, Axis &x2, Axis &x3, Axis &x4, Axis &x5 ){
      assert( this->axes.size() == 6 );
      std::vector<Axis> reordered_axes {x0, x1, x2, x3, x4, x5};
      if( !this->_same_axes(this->axes, reordered_axes) ){
        fmt::print("the input axis is not same as axes in kernel when reorder!");
        exit(-1);
      }
      this->axes[0] = x0;
      this->axes[1] = x1;
      this->axes[2] = x2;
      this->axes[3] = x3;
      this->axes[4] = x4;
      this->axes[5] = x5;
    }
    void parallel(Axis &x, int num_threads){
      assert( this->axes.size()>0 );
      //we can only parallel the outer most for loop with multi-threads.
      assert( this->axes[0] == x );
      assert( num_threads >=0 );
      (this->axes[0]).set_num_threads(num_threads);
    }
    
    void cache_read (const FrontendTensor &x, CacheRead & _cacheread, const std::string &position) {
      //the allocated buffer in ldm is global by default.
      assert(position == "global");
      _cacheread.set(x.SpNodePtr);
      this->cacheread = _cacheread;
      this->SpNodePtr = x.SpNodePtr;
    }
    void cache_write( CacheWrite & _cachewrite, const std::string &position) {
      //the allocated buffer in ldm is global by default.
      assert(position == "global");
      _cachewrite.set(this->TeNodePtr);
      this->cachewrite = _cachewrite;
    }
    void compute_at( const CacheRead & _cacheread, const Axis & x ){
      // inorder to compute on CPE of Sunway, we only tile the outer-most loop for continuous data acess.
      // assert(this->cacheread == _cacheread && this->axes.size()>0 && this->axes[0] == x);
      this->cacheread.set(x);
    }
    void compute_at( const CacheWrite & _cachewrite, const Axis & x ){
      // inorder to compute on CPE of Sunway, we only tile the outer-most loop for continuous data acess.
      // assert(this->cachewrite == _cachewrite && this->axes.size()>0 && this->axes[0] == x);
      this->cachewrite.set(x);
    }
    void build( std::string target ){
      if( target == "sunway" ) {
        //scheduleforstmt
        std::shared_ptr<Block> body_scheduleforstmt = std::make_shared<Block>();
        body_scheduleforstmt->insert(std::dynamic_pointer_cast<IRNode>(this->_root));
        //
        assert(this->cacheread.use_cache_read() && this->cachewrite.use_cache_write());
        auto scheduleforstmt = std::make_shared<ScheduleForStmt>(this->_eg, this->axes, this->cacheread, this->cachewrite, body_scheduleforstmt);

        //
        VecInt shape = this->_root->shape();
        DataType dt = this->_root->dt();
        //
        auto arg_temp = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
        std::vector<std::pair<DataType, std::shared_ptr<Expr>> > args;
        args.push_back( make_pair(DataType::i32, arg_temp) );
        KernelIRVisitor::run( (this->_root).get(), arg_temp );
        args.push_back( make_pair(dt, this->TeNodeExprPtr) );

        //for slave code
        auto arg_struct_slave = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(dt));
        std::vector<std::pair<DataType, std::shared_ptr<Expr>> > args_slave;
        args_slave.push_back( make_pair(dt, arg_struct_slave) );
        std::shared_ptr<Block> body_func_slave_stmt = std::make_shared<Block>();
        body_func_slave_stmt->insert(std::dynamic_pointer_cast<IRNode>(scheduleforstmt));
        auto FunctionStmtSlavePtr= std::make_shared<FunctionStmt>(args_slave, body_func_slave_stmt);
        //
        auto extern_func = std::make_pair(FunctionStmtSlavePtr->id, dt);
        call_extern_func_slave.push_back(extern_func);
        //
        //slave_sunway_ast_builder
        slave_sunway_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(FunctionStmtSlavePtr));

        //
        //AthreadStmt
        std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
        std::shared_ptr<AthreadStmt> athreadstmtptr = std::make_shared<AthreadStmt>(FunctionStmtSlavePtr->id, arg_temp, this->SpNodePtr, this->TeNodePtr);
        body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(athreadstmtptr));
        this->FunctionStmtPtr = std::make_shared<FunctionStmt>(args, body_funcstmt);
        //current_ast_builder
        current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(this->FunctionStmtPtr));
      }
      else{
        //scheduleforstmt
        std::shared_ptr<Block> body_scheduleforstmt = std::make_shared<Block>();
        body_scheduleforstmt->insert(std::dynamic_pointer_cast<IRNode>(this->_root));
        //
        CacheRead _cacheread;
        CacheWrite _cachewrite;
        auto scheduleforstmt = std::make_shared<ScheduleForStmt>(this->_eg, this->axes, _cacheread, _cachewrite, body_scheduleforstmt);

        //
        VecInt shape = this->_root->shape();
        DataType dt = this->_root->dt();
        //
        auto arg_temp = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
        std::vector<std::pair<DataType, std::shared_ptr<Expr>> > args;
        args.push_back( make_pair(DataType::i32, arg_temp) );
        KernelIRVisitor::run( (this->_root).get(), arg_temp );
        args.push_back( make_pair(dt, this->TeNodeExprPtr) );
        //
        std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
        body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(scheduleforstmt));
        this->FunctionStmtPtr = std::make_shared<FunctionStmt>(args, body_funcstmt);
        current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(this->FunctionStmtPtr));
      }
    }

    Kernel(ExprGroup eg, ExprPtr root, int num_threads=0 ){
      assert(num_threads >= 0);
      //TODO
      //pass transform shape of spatial to temporal
      TypeCheck::run( root.get() );
      VecInt shape = root->shape();
      DataType dt = root->dt();
      //this->SpNodePtr = std::make_shared<SpNode>(shape.size(), shape, dt);
      this->SpNodePtr = nullptr;
      //
      auto TeNodePtr = std::make_shared<TeNode>(0, shape.size(), shape, dt);
      this->TeNodePtr = TeNodePtr;
      std::shared_ptr<Expr> TeNodeExprPtr = std::dynamic_pointer_cast<Expr>(std::make_shared<TeNodeExpr>(TeNodePtr));
      std::shared_ptr<Expr> TeNodeIndicesExprPtr = std::dynamic_pointer_cast<Expr>(std::make_shared<TeNodeIndicesExpr>(TeNodePtr, eg));
      //
      auto kernelexpr = TeNodeIndicesExprPtr << root;
      //
      this->_root=std::move(kernelexpr);
      this->_eg = eg;
      //TypeCheck for two times
      TypeCheck::run( (this->_root).get() );

      // to get TeNodeExprPtr for cache_read and cache_write
      this->TeNodeExprPtr = TeNodeExprPtr;

      std::shared_ptr<Block> body_frontendforstmt = std::make_shared<Block>();
      body_frontendforstmt->insert(std::dynamic_pointer_cast<IRNode>(this->_root));
      //
      std::vector<int> begins;
      for(size_t i=0;i<shape.size();i++){
        begins.push_back(0);
      }
      auto frontendforstmt = std::make_shared<FrontendForStmt>(eg, begins, shape, num_threads, body_frontendforstmt);

      // args should be a temporary writed buffer and the shape and dt are infered
      // there are two args, t (which time step) and buffer (writing results).
      auto arg_temp = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      std::vector<std::pair<DataType, std::shared_ptr<Expr>> > args;
      args.push_back( make_pair(DataType::i32, arg_temp) );
      KernelIRVisitor::run( (this->_root).get(), arg_temp );
      args.push_back( make_pair(dt, TeNodeExprPtr) );
      //
      std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
      body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(frontendforstmt));
      this->FunctionStmtPtr = std::make_shared<FunctionStmt>(args, body_funcstmt);
      current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(this->FunctionStmtPtr));
    }
    ~Kernel(){}

    bool is_spatial(){
      return !((this->_eg).empty());
    }
    void compile(){}

    void print_ir(){
      run(_root.get());
    }
    void run(IRNode * root){
      IRPrinter::run( root );      
    }

    ExprPtr operator[](int temporal){
      return this->temporal(temporal);
    }

    ExprPtr temporal(int temporal){
      VecInt shape = this->_root->shape();
      DataType dt = this->_root->dt();
      auto id_func = this->FunctionStmtPtr->id;
      std::shared_ptr<Expr> _TempExprPtr;

      if(found.find(std::make_pair(id_func, temporal))==found.end() ){
        auto TeNodePtr = std::make_shared<TeNode>(0, shape.size(), shape, dt);
        auto FrontendAllocaStmtPtr = std::make_shared<FrontendAllocaStmt>(shape, dt, TeNodePtr->id); 
        auto TeNodeIndicesExprPtr = std::make_shared<TeNodeIndicesExpr>(TeNodePtr, this->_eg);
        //
        auto arg_temp = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(temporal));
        std::vector<std::shared_ptr<Expr> > args;
        args.push_back(arg_temp);
        auto CallNodeExprPtr = std::make_shared<CallNodeExpr>(TeNodePtr, shape, dt);
        //args.push_back(std::dynamic_pointer_cast<Expr>(TeNodeIndicesExprPtr));
        args.push_back(std::dynamic_pointer_cast<Expr>(CallNodeExprPtr));
        auto CallFuncExprPtr = std::make_shared<CallFuncExpr>(args, id_func, shape, dt);
        //
        
        std::shared_ptr<PlaceHolderExpr> PlaceHolderExprPtr = std::make_shared<PlaceHolderExpr>(FrontendAllocaStmtPtr, CallFuncExprPtr, TeNodeIndicesExprPtr);
        _TempExprPtr = std::dynamic_pointer_cast<Expr>(PlaceHolderExprPtr);
        //current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(FrontendAllocaStmtPtr));
        //current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(CallFuncExprPtr));
        found[std::make_pair(id_func, temporal)]=_TempExprPtr;
        
      }
      else {
        _TempExprPtr=found[std::make_pair(id_func, temporal)];
      }
      return _TempExprPtr;
    }

    /*
    ExprPtr _temporal(int temporal){
      VecInt shape = this->_root->shape();
      DataType dt = this->_root->dt();
      //auto args = this->FunctionStmtPtr->args;
      auto arg_temp = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(temporal));
      std::vector<std::shared_ptr<Expr> > args;
      args.push_back(arg_temp);
      auto id_func = this->FunctionStmtPtr->id;
      auto callfuncexpr = std::dynamic_pointer_cast<Expr>(std::make_shared<CallFuncExpr>(args, id_func, shape, dt));
      std::shared_ptr<Expr> tenodeexpr;
      if(found.find(std::make_pair(id_func, temporal))==found.end() ){
        //
        auto TeNodePtr = std::make_shared<TeNode>(0, shape.size(), shape, dt);
        current_ast_builder().insert(std::move(std::dynamic_pointer_cast<IRNode>(std::make_shared<FrontendAllocaStmt>(shape , dt, TeNodePtr->id))));
        //
        tenodeexpr = std::dynamic_pointer_cast<Expr>(std::make_shared<TeNodeExpr>(temporal, shape, dt));
        auto allocstmt = tenodeexpr << callfuncexpr;
        current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(allocstmt));
        found[std::make_pair(id_func, temporal)]=tenodeexpr;
      }
      else {
        tenodeexpr=found[std::make_pair(id_func, temporal)];
      }
      //
      //return callfuncexpr;
      return tenodeexpr;
    }
    */
};


class Stencil{
  public:
    ExprPtr _root;
    std::shared_ptr<TeNodeExpr> TeNodeExprPtr;
    std::shared_ptr<FunctionStmt> FunctionStmtPtr;
    std::vector<int> vec_run_id;
    bool use_mpi;
    std::vector<int> shape_mpi;
    int func_id_timing;
    int start;
    int end;

    static int t;

    Stencil(){}
    //this is a old one and does not fit for asyn
    /*
    Stencil(std::vector<int> shape_mpi, const FrontendTensor &x, ExprGroup eg, ExprPtr root){
      this->_root=std::move(root);
      TypeCheck::run( (this->_root).get() );
      GraphIRVisitor::run( (this->_root).get() ); 
      //
      VecInt shape = this->_root->shape();
      DataType dt = this->_root->dt();
      std::cout<<"shape size in graph is "<<shape.size()<<std::endl;
      this->TeNodeExprPtr = std::make_shared<TeNodeExpr>(this->t, shape, dt);
      // 
      auto loop_var = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      //GraphIRVisitor::set( (this->_root).get(), loop_var ); 
      GraphIRVisitorPlus::run( (this->_root).get(), loop_var );
      // Alloca and CallFunc in RangeForStmt
      std::shared_ptr<Block> body_rangeforstmt = std::make_shared<Block>();
      std::cout<<"PlaceHolderExprPtrVecs.size(): "<<GraphIRVisitor::PlaceHolderExprPtrVecs.size()<<std::endl;
      for(auto ptr: GraphIRVisitor::PlaceHolderExprPtrVecs){
        body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(ptr->FrontendAllocaStmtPtr));
        ptr->CallFuncExprPtr->args[0]=loop_var+(ptr->CallFuncExprPtr->args[0]);
        body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(ptr->CallFuncExprPtr));
      }

      //TODO
      //insert IfStmt into body_frontendforstmt
      //body_frontendforstmt
      std::shared_ptr<Block> body_frontendforstmt = std::make_shared<Block>();
      body_frontendforstmt->insert(std::dynamic_pointer_cast<IRNode>(this->_root));
      //begins and ends of intra and outer
      std::vector<int> begins_intra;
      std::vector<int> begins_outer;
      std::vector<int> ends_intra;
      std::vector<int> ends_outer = shape;
      for(size_t i=0;i<shape.size();i++){
        begins_intra.push_back((x.SpNodePtr)->halo_size());
        begins_outer.push_back(0);
        ends_intra.push_back(shape[i]-((x.SpNodePtr)->halo_size()));
      }
      //frontendforstmt of intra and outer
      auto frontendforstmt_intra = std::make_shared<FrontendForStmt>(eg, begins_intra, ends_intra, body_frontendforstmt);
      auto frontendforstmt_outer = std::make_shared<FrontendForStmt>(eg, begins_outer, ends_outer, body_frontendforstmt);

      //exchange_halo_start
      std::vector<std::shared_ptr<Expr> > args_mpiexpr;
      args_mpiexpr.push_back( (loop_var - 1) );
      auto mpiexprptr_start= std::make_shared<MPIExpr>(args_mpiexpr, x.SpNodePtr, shape_mpi, 1, shape, dt);
      body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(mpiexprptr_start));
      //intra area
      body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(frontendforstmt_intra));
      //exchange_halo_end
      auto mpiexprptr_end= std::make_shared<MPIExpr>(args_mpiexpr, x.SpNodePtr, shape_mpi, 2, shape, dt);
      body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(mpiexprptr_end));
      //outer area
      body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(frontendforstmt_outer));

      // rangefor
      auto arg_start = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      auto arg_end = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      auto rangeforstmtptr = std::make_shared<RangeForStmt>(loop_var, arg_start, arg_end, body_rangeforstmt, 0, 0);
      //
      // start func
      std::vector<std::pair<DataType, std::shared_ptr<Expr>> > args;
      args.push_back( make_pair(DataType::i32, arg_start) );
      args.push_back( make_pair(DataType::i32, arg_end) );
      std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
      body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(rangeforstmtptr));
      auto FunctionStmtPtr = std::make_shared<FunctionStmt>(args, body_funcstmt);
      // end func
      current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(FunctionStmtPtr));
      this->FunctionStmtPtr = FunctionStmtPtr;
    }
    */

    Stencil(std::vector<int> shape_mpi, const FrontendTensor &x, ExprGroup eg, ExprPtr root, int num_threads=0 , int is_assigned=1){
      assert(num_threads >= 0);
      this->use_mpi = true;
      this->shape_mpi = shape_mpi;
      //
      this->_root=std::move(root);
      TypeCheck::run( (this->_root).get() );
      GraphIRVisitor::run( (this->_root).get() ); 
      //
      VecInt shape = this->_root->shape();
      DataType dt = this->_root->dt();
      std::cout<<"shape size in graph is "<<shape.size()<<std::endl;
      this->TeNodeExprPtr = std::make_shared<TeNodeExpr>(this->t, shape, dt);
      // 
      auto loop_var = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      //GraphIRVisitor::set( (this->_root).get(), loop_var ); 
      GraphIRVisitorPlus::run( (this->_root).get(), loop_var );
      // Alloca and CallFunc in RangeForStmt
      std::shared_ptr<Block> body_rangeforstmt = std::make_shared<Block>();
      std::cout<<"PlaceHolderExprPtrVecs.size(): "<<GraphIRVisitor::PlaceHolderExprPtrVecs.size()<<std::endl;
      for(auto ptr: GraphIRVisitor::PlaceHolderExprPtrVecs){
        //allocate TeNode global
        current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(ptr->FrontendAllocaStmtPtr));
        //origin allocate TeNode in function
        //body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(ptr->FrontendAllocaStmtPtr));
        ptr->CallFuncExprPtr->args[0]=loop_var+(ptr->CallFuncExprPtr->args[0]);
        body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(ptr->CallFuncExprPtr));
      }

      //TODO
      //body_frontendforstmt
      std::shared_ptr<Block> body_frontendforstmt = std::make_shared<Block>();
      body_frontendforstmt->insert(std::dynamic_pointer_cast<IRNode>(this->_root));
      //begins and ends
      std::vector<int> begins;
      std::vector<int> ends = shape;
      for(size_t i=0;i<shape.size();i++){
        begins.push_back(0);
      }
      //frontendforstmt
      auto frontendforstmt = std::make_shared<FrontendForStmt>(eg, begins, ends, num_threads, body_frontendforstmt);
      //
      if(is_assigned==1){
        body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(frontendforstmt));
      }

      //exchange_halo_syn
      std::vector<std::shared_ptr<Expr> > args_mpiexpr;
      args_mpiexpr.push_back( loop_var );
      auto mpiexprptr_syn= std::make_shared<MPIExpr>(args_mpiexpr, x.SpNodePtr, shape_mpi, 3, shape, dt);
      body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(mpiexprptr_syn));

      // rangefor
      auto arg_start = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      auto arg_end = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      auto rangeforstmtptr = std::make_shared<RangeForStmt>(loop_var, arg_start, arg_end, body_rangeforstmt, 0, 0);
      //
      // start func
      std::vector<std::pair<DataType, std::shared_ptr<Expr>> > args;
      args.push_back( make_pair(DataType::i32, arg_start) );
      args.push_back( make_pair(DataType::i32, arg_end) );
      std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
      body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(rangeforstmtptr));
      auto FunctionStmtPtr = std::make_shared<FunctionStmt>(args, body_funcstmt);
      // end func
      current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(FunctionStmtPtr));
      this->FunctionStmtPtr = FunctionStmtPtr;
    }

    Stencil(ExprGroup eg, ExprPtr root, int num_threads=0 , int is_assigned=1){
      assert(num_threads >= 0);
      this->use_mpi = false;
      this->shape_mpi = std::vector<int> {};
      //
      this->_root=std::move(root);
      TypeCheck::run( (this->_root).get() );
      GraphIRVisitor::run( (this->_root).get() ); 
      //
      VecInt shape = this->_root->shape();
      DataType dt = this->_root->dt();
      std::cout<<"shape size in graph is "<<shape.size()<<std::endl;
      this->TeNodeExprPtr = std::make_shared<TeNodeExpr>(this->t, shape, dt);
      // 
      auto loop_var = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      //GraphIRVisitor::set( (this->_root).get(), loop_var ); 
      GraphIRVisitorPlus::run( (this->_root).get(), loop_var );
      //
      std::shared_ptr<Block> body_rangeforstmt = std::make_shared<Block>();
      std::cout<<"PlaceHolderExprPtrVecs.size(): "<<GraphIRVisitor::PlaceHolderExprPtrVecs.size()<<std::endl;
      for(auto ptr: GraphIRVisitor::PlaceHolderExprPtrVecs){
        //allocate TeNode global
        current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(ptr->FrontendAllocaStmtPtr));
        //origin allocate TeNode in function
        //body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(ptr->FrontendAllocaStmtPtr));
        ptr->CallFuncExprPtr->args[0]=loop_var+(ptr->CallFuncExprPtr->args[0]);
        body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(ptr->CallFuncExprPtr));
      }
      //
      std::shared_ptr<Block> body_frontendforstmt = std::make_shared<Block>();
      body_frontendforstmt->insert(std::dynamic_pointer_cast<IRNode>(this->_root));
      //
      std::vector<int> begins;
      for(size_t i=0;i<shape.size();i++){
        begins.push_back(0);
      }
      auto frontendforstmt = std::make_shared<FrontendForStmt>(eg, begins, shape, num_threads, body_frontendforstmt);
      //
      if(is_assigned==1){
        body_rangeforstmt->insert(std::dynamic_pointer_cast<IRNode>(frontendforstmt));
      }
      //
      auto arg_start = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      auto arg_end = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(DataType::i32));
      auto rangeforstmtptr = std::make_shared<RangeForStmt>(loop_var, arg_start, arg_end, body_rangeforstmt, 0, 0);
      //
      // start func
      std::vector<std::pair<DataType, std::shared_ptr<Expr>> > args;
      args.push_back( make_pair(DataType::i32, arg_start) );
      args.push_back( make_pair(DataType::i32, arg_end) );
      std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
      body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(rangeforstmtptr));
      auto FunctionStmtPtr = std::make_shared<FunctionStmt>(args, body_funcstmt);
      // end func
      current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(FunctionStmtPtr));
      this->FunctionStmtPtr = FunctionStmtPtr;
    }

    Stencil(ExprPtr root){
      this->use_mpi = false;
      this->shape_mpi = std::vector<int> {};
      //
      this->_root=std::move(root);
      TypeCheck::run( (this->_root).get() );
      GraphIRVisitor::run( (this->_root).get() ); 
      //
      VecInt shape = this->_root->shape();
      DataType dt = this->_root->dt();
      this->TeNodeExprPtr = std::make_shared<TeNodeExpr>(this->t, shape, dt);
      //
      current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(this->_root));
    }
    ~Stencil(){}

    void print_ir(){
      IRPrinter::run( (this->_root).get() );      
    }

    void run(int start , int end){
      this->start = start;
      this->end = end;
      auto arg_start = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(start));
      auto arg_end = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(end));
      std::vector<std::shared_ptr<Expr> > args;
      args.push_back(arg_start);
      args.push_back(arg_end);
      auto CallFuncExprPtr = std::make_shared<CallFuncExpr>(args, this->FunctionStmtPtr->id, this->TeNodeExprPtr->shape(), this->TeNodeExprPtr->dt());
      //
      std::vector<std::pair<DataType, std::shared_ptr<Expr> > > args_func;
      std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
      body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(CallFuncExprPtr));
      auto FunctionStmtPtr = std::make_shared<FunctionStmt>(args_func, body_funcstmt);
      this->vec_run_id.push_back( FunctionStmtPtr->id );
      this->func_id_timing = FunctionStmtPtr->id;
      //
      //current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(CallFuncExprPtr));
      current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(FunctionStmtPtr));
    }
    // mpi
    void input(std::vector<int> shape_mpi, const FrontendTensor &x, std::string name_file){
      //auto arg_time_step = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(time_step));
      std::vector<std::shared_ptr<Expr> > args;
      //args.push_back(arg_time_step);
      auto InputExprPtr = std::make_shared<InputExpr>(args, x.SpNodePtr, name_file, shape_mpi, this->TeNodeExprPtr->shape(), this->TeNodeExprPtr->dt());
      //body_funcstmt
      std::vector<std::pair<DataType, std::shared_ptr<Expr> > > args_func;
      std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
      body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(InputExprPtr));
      //exchange_halo_syn
      int temp_size = (x.SpNodePtr)->temp_size();
      if((x.SpNodePtr)->use_time_window()){
        temp_size = (x.SpNodePtr)->time_window_size();
      }
      for(int i=0; i<temp_size; i++){
        std::vector<std::shared_ptr<Expr> > args_mpiexpr;
        auto arg_time_step = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(i));
        args_mpiexpr.push_back( arg_time_step );
        auto mpiexprptr_syn= std::make_shared<MPIExpr>(args_mpiexpr, x.SpNodePtr, shape_mpi, 3, (x.SpNodePtr)->shape(), (x.SpNodePtr)->dt());
        body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(mpiexprptr_syn));
      }
      //
      auto FunctionStmtPtr = std::make_shared<FunctionStmt>(args_func, body_funcstmt);
      this->vec_run_id.push_back( FunctionStmtPtr->id );
      //
      current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(FunctionStmtPtr));
    }
    void input(const FrontendTensor &x, std::string name_file){
      //auto arg_time_step = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(time_step));
      std::vector<std::shared_ptr<Expr> > args;
      //args.push_back(arg_time_step);
      auto InputExprPtr = std::make_shared<InputExpr>(args, x.SpNodePtr, name_file, this->TeNodeExprPtr->shape(), this->TeNodeExprPtr->dt());
      //
      std::vector<std::pair<DataType, std::shared_ptr<Expr> > > args_func;
      std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
      body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(InputExprPtr));
      auto FunctionStmtPtr = std::make_shared<FunctionStmt>(args_func, body_funcstmt);
      this->vec_run_id.push_back( FunctionStmtPtr->id );
      //
      current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(FunctionStmtPtr));
    }
    void output(const FrontendTensor &x, const int time_step, std::string name_file){
      auto arg_time_step = std::dynamic_pointer_cast<Expr>(std::make_shared<VariableExpr>(time_step));
      std::vector<std::shared_ptr<Expr> > args;
      args.push_back(arg_time_step);
      auto OutputExprPtr = std::make_shared<OutputExpr>(args, x.SpNodePtr, name_file, this->TeNodeExprPtr->shape(), this->TeNodeExprPtr->dt());
      //
      std::vector<std::pair<DataType, std::shared_ptr<Expr> > > args_func;
      std::shared_ptr<Block> body_funcstmt = std::make_shared<Block>();
      body_funcstmt->insert(std::dynamic_pointer_cast<IRNode>(OutputExprPtr));
      auto FunctionStmtPtr = std::make_shared<FunctionStmt>(args_func, body_funcstmt);
      this->vec_run_id.push_back( FunctionStmtPtr->id );
      //
      current_ast_builder().insert(std::dynamic_pointer_cast<IRNode>(FunctionStmtPtr));
    }
    void compile(){
      std::string Code_C = CodeGenCPU::run( context->root() , this->vec_run_id );
      fmt::print(Code_C);
    }
    int _num_processes(){
      int num_processes=1;
      if(this->use_mpi){
        for(auto v: this->shape_mpi){
          num_processes *= v;
        }
      }
      else{
        fmt::print("there is no shape of mpi!");
        exit(-1);
      }
      return num_processes;
    }
    void _target(std::string name_file, std::string target){
      //
      std::string cc = "gcc ";
      std::string omp_flag = "-fopenmp -O3 ";
      std::string run = " ";
      std::string run_flag = " ";
      //
      std::string name_compile_file = std::string("compile_")+name_file+".sh";
      std::string name_run_file = std::string("run_")+name_file+".sh";
      //
      if(target == "x86"){
      }
      else if(target == "macos"){
        omp_flag = std::string("-Xpreprocessor ") + omp_flag + "-lomp ";
      }
      else if(target == "feiteng" || target == "sunway" || target == "sunway_slave"){
        fmt::print(std::string("Please use mpi on ")+target+" !");
        exit(-1);
      }
      else{
        fmt::print(std::string("The supported targets are : x86, feiteng, sunway and macos! Please check spelling!"));
        exit(-1);
      }
      //
      std::string command_compile = cc + omp_flag + "-o " + name_file + " " + name_file + ".c";
      std::string command_run = run + run_flag + "./" + name_file;
      write2file(name_compile_file, command_compile);
      write2file(name_run_file, command_run);
    }

    void _target_mpi(std::string name_file, std::string target){
      int num_processes = this->_num_processes();
      //
      std::string cc = "mpicc ";
      std::string omp_flag = "-fopenmp -O3 ";
      std::string run = "mpirun ";
      std::string run_flag = std::string("-np ")+std::to_string(num_processes)+" ";
      //
      std::string name_compile_file = std::string("compile_")+name_file+".sh";
      std::string name_run_file = std::string("run_")+name_file+".sh";
      //
      if(target == "x86"){
        run_flag = std::string("--bind-to none ")+run_flag;
        //
      	std::string command_compile = cc + omp_flag + "-o " + name_file + " " + name_file + ".c";
      	std::string command_run = run + run_flag + "./" + name_file+" > "+name_file+".x86.out ";
      	write2file(name_compile_file, command_compile);
      	write2file(name_run_file, command_run);
        return ;
      }
      else if(target == "feiteng"){
        std::string command_compile = cc + omp_flag + "-o " + name_file + " " + name_file + ".c";
        std::string command_run = std::string("#!/bin/bash\nexport OMP_NUM_THREADS=32\n")+"yhrun -N "+std::to_string(num_processes)+" -n "+ std::to_string(num_processes)+" -c 32 -p th_mt "+"./" + name_file;
        std::string command_batch = std::string("yhbatch -N ")+std::to_string(num_processes)+" -n "+ std::to_string(num_processes)+" -p th_mt "+"./" + name_run_file;
        write2file(name_compile_file, command_compile);
        write2file(name_run_file, command_run);
        std::string name_batch_file = std::string("batch_")+name_file+".sh";
        write2file(name_batch_file, command_batch);
        return ;
      }
      else if(target == "macos"){
        omp_flag = std::string("-Xpreprocessor ") + omp_flag + "-lomp ";
      }
      else if(target == "sunway"){
        omp_flag = " ";
        std::string sw5cc = "sw5cc -host -I/usr/sw-mpp/mpi2/include -std=c99 -c ";
        std::string command_compile_host = sw5cc + name_file + ".c";
        std::string command_compile_hybrid = cc + "-hybrid " + "-o " + name_file + " " + name_file + ".o";
        std::string job_queue = std::string("q_sw_share ");
        std::string command_run = std::string("bsub -q ")+job_queue+"-b -cgsp 64 -n "+ std::to_string(num_processes) + " -host_stack 512 -share_size 7450 " + "-o ./" + name_file + ".out " + "./" + name_file;
        write2file(name_compile_file, command_compile_host+"\n"+command_compile_hybrid);
        write2file(name_run_file, command_run);
        return ;
      }
      else if(target == "sunway_slave"){
        omp_flag = " ";
        std::string sw5cc = "sw5cc -host -I/usr/sw-mpp/mpi2/include -std=c99 -c ";
        std::string command_compile_host = sw5cc + name_file + ".c";
        std::string command_compile_slave = std::string("sw5cc -slave -c ") + name_file + "_slave.c";
        std::string command_compile_hybrid = cc + "-hybrid " + "-o " + name_file + " " + name_file + ".o "+name_file+"_slave.o";
        std::string job_queue = std::string("q_sw_share ");
        //std::string command_run = std::string("bsub -q ")+job_queue+"-I -b -cgsp 64 -n "+ std::to_string(num_processes) + " -host_stack 128 -share_size 2048 " + "./" + name_file;
        std::string command_run = std::string("bsub -q ")+job_queue+"-I -b -cgsp 64 -n "+ std::to_string(num_processes) + " -host_stack 512 -share_size 7450 " + "-o ./" + name_file + ".out " + "./" + name_file;
        write2file(name_compile_file, command_compile_host+"\n"+command_compile_slave+"\n"+command_compile_hybrid);
        write2file(name_run_file, command_run);
        return ;
      }
      else{
        fmt::print(std::string("The supported targets are : x86, feiteng, sunway and macos! Please check spelling!"));
        exit(-1);
      }
      //
      std::string command_compile = cc + omp_flag + "-o " + name_file + " " + name_file + ".c";
      std::string command_run = run + run_flag + "./" + name_file;
      write2file(name_compile_file, command_compile);
      write2file(name_run_file, command_run);
    }
    void compile_to_source_code(std::string name_file){
      std::string Code_C = CodeGenCPU::run( context->root() , this->vec_run_id );
      name_file += ".c";
      std::ofstream	OsWrite(name_file, std::ofstream::out);
	    OsWrite<<Code_C;
	    OsWrite.close();
    }
    void compile_to_source_code(std::string name_file, std::string target = "", bool use_timing = false){
      //
      int timing_id =-1;
      if(use_timing) timing_id = this->func_id_timing;
      std::cout<<"timing id "<<timing_id<<std::endl;
      if(target == "sunway"){
        //host
        std::string Code_C = CodeGenSunway::run( context->root() , this->vec_run_id , this->use_mpi, timing_id, this->start, this->end);
        std::string name_file_host = name_file + ".c";
        std::ofstream	OsWrite(name_file_host, std::ofstream::out);
	      OsWrite<<Code_C;
	      OsWrite.close();
        //slave
        if(call_extern_func_slave.size()>0){
          std::string Code_C_Slave = CodeGenSunwaySlave::run( context_slave_sunway->root() );
          std::string name_file_slave = name_file + "_slave.c";
          write2file(name_file_slave, Code_C_Slave);
        }
      }
      else{
        std::string Code_C = CodeGenCPU::run( context->root() , this->vec_run_id , this->use_mpi, timing_id, this->start, this->end);
        std::ofstream	OsWrite(name_file + ".c", std::ofstream::out);
	      OsWrite<<Code_C;
	      OsWrite.close();
      }
      if(target != ""){
        if(target =="sunway" && call_extern_func_slave.size()>0)
          target+="_slave";
        if(this->use_mpi){
          this->_target_mpi(name_file, target);
        }
        else{
          this->_target(name_file, target);
        }
      }
    }
    void compile_to_source_code_mpi(std::string name_file, std::string target = "", bool use_timing = false){
      //
      int timing_id =-1;
      if(use_timing) timing_id = this->func_id_timing;
      std::cout<<"timing id "<<timing_id<<std::endl;
      if(target == "sunway"){
        //host
        std::string Code_C = CodeGenSunway::run( context->root() , this->vec_run_id , this->use_mpi, timing_id, this->start, this->end);
        std::string name_file_host = name_file + ".c";
        std::ofstream	OsWrite(name_file_host, std::ofstream::out);
	      OsWrite<<Code_C;
	      OsWrite.close();
        //slave
        if(call_extern_func_slave.size()>0){
          std::string Code_C_Slave = CodeGenSunwaySlave::run( context_slave_sunway->root() );
          std::string name_file_slave = name_file + "_slave.c";
          write2file(name_file_slave, Code_C_Slave);
        }
      }
      else{
        std::string Code_C = CodeGenCPU::run( context->root() , this->vec_run_id , this->use_mpi, timing_id, this->start, this->end);
        std::ofstream	OsWrite(name_file + ".c", std::ofstream::out);
	      OsWrite<<Code_C;
	      OsWrite.close();
      }
      if(target != ""){
        if(target =="sunway" && call_extern_func_slave.size()>0)
          target+="_slave";
        this->_target_mpi(name_file, target);
      }
    }

};

int Stencil::t = 0;




/*
class For {
 public:
  For(Expr i, Expr s, Expr e, const std::function<void()> &func) {
    auto stmt_unique = std::make_unique<FrontendForStmt>(i, s, e);
    auto stmt = stmt_unique.get();
    current_ast_builder().insert(std::move(stmt_unique));
    auto _ = current_ast_builder().create_scope(stmt->body);
    func();
  }
};

*/
