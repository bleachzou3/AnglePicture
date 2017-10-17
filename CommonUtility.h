#ifndef COMMON_UTILITY_HPP
#define COMMON_UTILITY_HPP
#include <sstream>
#include <string>
template<typename T> string toString(const T& t){
    ostringstream oss;  //����һ����ʽ�������
    oss<<t;             //��ֵ����������
    return oss.str();   
}





#endif