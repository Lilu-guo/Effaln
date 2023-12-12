#ifndef PROCESSOR_SUPPORT_H_
#define PROCESSOR_SUPPORT_H_
#if defined(__INTEL_COMPILER)
#define USING_INTEL_COMPILER
#elif defined(__GNUC__)
#define USING_GCC_COMPILER
#include <cpuid.h>
#elif defined(_MSC_VER)
#define USING MSC_COMPILER
#endif
struct regs_t
{
    unsigned int EAX, EBX, ECX, EDX;
};
#define BIT(n) ((1 << n))
class ProcessorSupport
{
#ifdef POPCNT_CAPABILITY
public:
    ProcessorSupport() {}
    bool POPCNTenabled()
    {
        regs_t regs;
        try
        {
#if (defined(USING_INTEL_COMPILER) || defined(USING_MSC_COMPILER))
            __cpuid((int *)&regs, 0);
            __cpuid((int *)&regs, 0x1);
#elif defined(USING_GCC_COMPILER)
                __get_cpuid(0x1, &regs.EAX, &regs.EBX, &regs.ECX, &regs.EDX);
#else
            std::cerr << "ERROR: please define __cpuid() for this build.\n";
            assert(0);
#endif
            if (!((regs.ECX & BIT(20)) && (regs.ECX & BIT(23))))
                return false;
        }
        catch (int e)
        {
            return false;
        }
        return true;
    }
#endif
};
#endif
