#ifndef PHI_CPP
#define PHI_CPP

#include "phi.h"
#include "spt_func.h"

double phi(double x)
{	
	return 1. / sqrt(1. - sqr_(x) );

}
#endif
