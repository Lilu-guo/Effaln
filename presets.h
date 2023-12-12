#ifndef PRESETS_H_
#define PRESETS_H_
#include <string>
#include <utility>
#include "ds.h"
class Presets
{
public:
	Presets() {}
	virtual ~Presets() {}
	virtual void apply(
		const std::string &preset,
		std::string &policy,
		EList<std::pair<int, std::string>> &opts) = 0;
	virtual const char *name() = 0;
};
class PresetsV0 : public Presets
{
public:
	PresetsV0() : Presets() {}
	virtual ~PresetsV0() {}
	virtual void apply(
		const std::string &preset,
		std::string &policy,
		EList<std::pair<int, std::string>> &opts);
	virtual const char *name() { return "V0"; }
};
#endif
