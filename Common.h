#ifndef COMMON_H_
#define COMMON_H_
#include <log4cpp/Category.hh>
#include <log4cpp/PropertyConfigurator.hh>
log4cpp::Category& rootLog  = log4cpp::Category::getRoot();
log4cpp::Category& subLog = log4cpp::Category::getInstance(std::string("sub1"));

#endif