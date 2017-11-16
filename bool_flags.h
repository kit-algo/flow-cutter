#ifndef BOOL_FLAGS_H
#define BOOL_FLAGS_H

#include <vector>
#include <cassert>
#include <algorithm>

class BoolFlags{
public:
	BoolFlags(){}
	explicit BoolFlags(int flag_count)
		:flag(flag_count, false){
	}

	int flag_count()const{
		return flag.size();
	}

	void set(int f){
		assert(flag_count () != 0 && "You forgot to pass the flag_count to the constructor");
		assert(0 <= f && f < flag_count() && "flag is out of bounds");
		flag[f] = true;
	}

	void unset(int f){
		assert(flag_count () != 0 && "You forgot to pass the flag_count to the constructor");
		assert(0 <= f && f < flag_count() && "flag is out of bounds");
		flag[f] = false;
	}

	bool is_set(int f)const{
		assert(flag_count () != 0 && "You forgot to pass the flag_count to the constructor");
		assert(0 <= f && f < flag_count() && "flag is out of bounds");
		return flag[f];
	}

	void unset_all(){
		std::fill(flag.begin(), flag.end(), false);
	}
private:
	std::vector<bool>flag;
};
#endif

