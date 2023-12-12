#ifndef SIMPLE_FUNC_H_
#define SIMPLE_FUNC_H_
#include <math.h>
#include <cassert>
#include <limits>
#include "tokenize.h"
#define SIMPLE_FUNC_CONST 1
#define SIMPLE_FUNC_LINEAR 2
#define SIMPLE_FUNC_SQRT 3
#define SIMPLE_FUNC_LOG 4
class SimpleFunc
{
public:
	SimpleFunc() : type_(0), I_(0.0), X_(0.0), C_(0.0), L_(0.0) {}
	SimpleFunc(int type, double I, double X, double C, double L)
	{
		init(type, I, X, C, L);
	}
	void init(int type, double I, double X, double C, double L)
	{
		type_ = type;
		I_ = I;
		X_ = X;
		C_ = C;
		L_ = L;
	}
	void init(int type, double C, double L)
	{
		type_ = type;
		C_ = C;
		L_ = L;
		I_ = -std::numeric_limits<double>::max();
		X_ = std::numeric_limits<double>::max();
	}
	void setType(int type) { type_ = type; }
	void setMin(double mn) { I_ = mn; }
	void setMax(double mx) { X_ = mx; }
	void setConst(double co) { C_ = co; }
	void setCoeff(double ce) { L_ = ce; }
	int getType() const { return type_; }
	double getMin() const { return I_; }
	double getMax() const { return X_; }
	double getConst() const { return C_; }
	double getCoeff() const { return L_; }
	void mult(double x)
	{
		if (I_ < std::numeric_limits<double>::max())
		{
			I_ *= x;
			X_ *= x;
			C_ *= x;
			L_ *= x;
		}
	}
	bool initialized() const { return type_ != 0; }
	void reset() { type_ = 0; }
	bool alwaysPositive() const
	{
		return f<double>(1.0) > 0 && (SIMPLE_FUNC_CONST || L_ >= 0.0);
	}
	template <typename T>
	T f(double x) const
	{
		assert(type_ >= SIMPLE_FUNC_CONST && type_ <= SIMPLE_FUNC_LOG);
		double X;
		if (type_ == SIMPLE_FUNC_CONST)
		{
			X = 0.0;
		}
		else if (type_ == SIMPLE_FUNC_LINEAR)
		{
			X = x;
		}
		else if (type_ == SIMPLE_FUNC_SQRT)
		{
			X = sqrt(x);
		}
		else if (type_ == SIMPLE_FUNC_LOG)
		{
			X = log(x);
		}
		else
		{
			throw 1;
		}
		double ret = std::max(I_, std::min(X_, C_ + L_ * X));
		if (ret == std::numeric_limits<double>::max())
		{
			return std::numeric_limits<T>::max();
		}
		else if (ret == std::numeric_limits<double>::min())
		{
			return std::numeric_limits<T>::min();
		}
		else
		{
			return (T)ret;
		}
	}
	static int parseType(const std::string &otype);
	static SimpleFunc parse(
		const std::string &s,
		double defaultConst = 0.0,
		double defaultLinear = 0.0,
		double defaultMin = 0.0,
		double defaultMax = std::numeric_limits<double>::max());
protected:
	int type_;
	double I_, X_, C_, L_;
};
#endif
