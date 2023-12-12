#include <iostream>
#include "presets.h"
#include "opts.h"
using namespace std;
void PresetsV0::apply(
	const std::string &preset,
	std::string &policy,
	EList<std::pair<int, std::string>> &opts)
{
	if (preset == "very-fast")
	{
		policy += ";SEED=0";
		policy += ";SEEDLEN=22";
		policy += ";DPS=5";
		policy += ";ROUNDS=1";
		policy += ";IVAL=S,0,2.50";
	}
	else if (preset == "fast")
	{
		policy += ";SEED=0";
		policy += ";SEEDLEN=22";
		policy += ";DPS=10";
		policy += ";ROUNDS=2";
		policy += ";IVAL=S,0,2.50";
	}
	else if (preset == "sensitive")
	{
		policy += ";SEED=0";
		policy += ";SEEDLEN=22";
		policy += ";DPS=15";
		policy += ";ROUNDS=2";
		policy += ";IVAL=S,1,1.15";
	}
	else if (preset == "very-sensitive")
	{
		policy += ";SEED=0";
		policy += ";SEEDLEN=20";
		policy += ";DPS=20";
		policy += ";ROUNDS=3";
		policy += ";IVAL=S,1,0.50";
	}
	else if (preset == "very-fast-local")
	{
		policy += ";SEED=0";
		policy += ";SEEDLEN=25";
		policy += ";DPS=5";
		policy += ";ROUNDS=1";
		policy += ";IVAL=S,1,2.00";
	}
	else if (preset == "fast-local")
	{
		policy += ";SEED=0";
		policy += ";SEEDLEN=22";
		policy += ";DPS=10";
		policy += ";ROUNDS=2";
		policy += ";IVAL=S,1,1.75";
	}
	else if (preset == "sensitive-local")
	{
		policy += ";SEED=0";
		policy += ";SEEDLEN=20";
		policy += ";DPS=15";
		policy += ";ROUNDS=2";
		policy += ";IVAL=S,1,0.75";
	}
	else if (preset == "very-sensitive-local")
	{
		policy += ";SEED=0";
		policy += ";SEEDLEN=20";
		policy += ";DPS=20";
		policy += ";ROUNDS=3";
		policy += ";IVAL=S,1,0.50";
	}
	else
	{
		cerr << "Unknown preset: " << preset.c_str() << endl;
	}
}
