#ifndef COMMON_UTILITY_HPP
#define COMMON_UTILITY_HPP
#include <sstream>
#include <string>
template<typename T> string toString(const T& t){
    ostringstream oss;  //创建一个格式化输出流
    oss<<t;             //把值传递如流中
    return oss.str();   
}





#endif