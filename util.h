
#include <iostream>
#include <cstring>

namespace fmt{
  void print(const std::string &str){
    std::cout<<str;
  }

  template<typename Type>
  void print(Type * array, int size){
    std::cout<<"\n[";
    for(int i=0;i< size-1;i++){
      std::cout<<array[i]<<", ";
    }
    std::cout<<array[size-1]<<"";
    std::cout<<"]\n";
  }
}

typedef std::vector<int> VecInt;

#define SHAPE_SCALAR (VecInt {-1})

bool same_vec(VecInt a, VecInt b){
  if(a.size() != b.size())return false;
  bool same=true;
  for(size_t i=0; i<a.size(); i++){
    if(a[i] != b[i]){
      same=false;
      break;
    }
  }
  return same;
}

bool is_scalar(VecInt &a){
  return same_vec(a, SHAPE_SCALAR);
}

bool same_datatype(DataType a, DataType b){
  return a==b;
}

bool is_float(DataType a){
  return a==DataType::f16 || a==DataType::f32 || a==DataType::f64;
}

bool is_int(DataType a){
  return a==DataType::i1 || a==DataType::i8 || a==DataType::i16 || a==DataType::i32 || a==DataType::i64;
}

bool is_uint(DataType a){
  return a==DataType::u8 || a==DataType::u16 || a==DataType::u32 || a==DataType::u64;
}

DataType convert_up(DataType a, DataType b){
/*
  if(a>b)
    std::cout<<"DataType a>b "<<std::endl;
  else
    std::cout<<"DataType a<=b "<<std::endl;
    */
  return a>b?a:b;
}

void write2file(std::string name_file, std::string content){
  std::ofstream	OsWrite(name_file, std::ofstream::out);
	OsWrite<<content;
	OsWrite.close();
}

template<typename Type>
bool in_vec(const std::vector<Type> &vec, Type &a){
  bool found = false;
  for(auto v: vec){
    if( v==a )found =true;
  }
  return found;
}


