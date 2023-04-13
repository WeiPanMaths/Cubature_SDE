#ifndef _TRAPFPE_H_
#define _TRAPFPE_H_
#ifdef GCC
#include <fenv.h>
#include <sys/cdefs.h>
__BEGIN_DECLS
void trapfpe(void){
	feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW)}
__END_DECLS

#else

#include <float.h>
void trapfpe(void)
{
	unsigned int cw;
	/* could use _controlfp */
	cw = _control87(0, 0) & MCW_EM;
	cw &= ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW);
	_control87(cw, MCW_EM);

}

#endif
#endif /* _TRAPFPE_H_ */
