#ifndef STATUS_LOG_H
#define STATUS_LOG_H

#include "timer.h"
#include <iostream>
#include <string>

class StatusLog{
	StatusLog()=delete;
	StatusLog(const StatusLog&)=delete;
	StatusLog&operator=(const StatusLog&)=delete;
	long long time;
public:
	StatusLog(std::string msg){
		std::cout << msg << " ... "<<std::flush;
		time = -get_micro_time();
	}
	~StatusLog(){
		time += get_micro_time();
		std::cout << "done ["<<time/1000 << "ms]"<< std::endl;
	}
};

#endif
