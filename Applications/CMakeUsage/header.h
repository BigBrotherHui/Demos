#include <iostream>
#include <cmakeusage_export.h>
#include "config.h"
class CMAKEUSAGE_EXPORT Test
{
public:
	void what()
	{
		std::cout << "hello,world" << std::endl;
		std::cout << var1 << std::endl;
	}
};