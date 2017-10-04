#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include "sequence.h"
#include "stdinc.h"
#include "afnio.h"
#include "cmsa.h"
#include "cma_gmb.h"
#include "dms_typ.h"

extern int     fa2cma_main(int argc, char *argv[],FILE *ofp);

int	main(int argc, char *argv[]){ return fa2cma_main(argc, argv,0); }


