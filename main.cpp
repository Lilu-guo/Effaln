#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include "tokenize.h"
#include "ds.h"
using namespace std;
extern "C"
{
	int effaln(int argc, const char **argv);
}
int main(int argc, const char **argv)
{
	if (argc > 2 && strcmp(argv[1], "-A") == 0)
	{
		const char *file = argv[2];
		ifstream in;
		in.open(file);
		char buf[4096];
		int lastret = -1;
		while (in.getline(buf, 4095))
		{
			EList<string> args;
			args.push_back(string(argv[0]));
			tokenize(buf, " \t", args);
			const char **myargs = (const char **)malloc(sizeof(char *) * args.size());
			for (size_t i = 0; i < args.size(); i++)
			{
				myargs[i] = args[i].c_str();
			}
			if (args.size() == 1)
				continue;
			lastret = effaln((int)args.size(), myargs);
			free(myargs);
		}
		if (lastret == -1)
		{
			cerr << "Warning: No arg strings parsed from " << file << endl;
			return 0;
		}
		return lastret;
	}
	else
	{
		return effaln(argc, argv);
	}
}
