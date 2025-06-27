#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _GluTrans_reg();
extern void _Kir4_reg();
extern void _kdrglia_reg();
extern void _kpump_reg();
extern void _potassiumd_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," GluTrans.mod");
fprintf(stderr," Kir4.mod");
fprintf(stderr," kdrglia.mod");
fprintf(stderr," kpump.mod");
fprintf(stderr," potassiumd.mod");
fprintf(stderr, "\n");
    }
_GluTrans_reg();
_Kir4_reg();
_kdrglia_reg();
_kpump_reg();
_potassiumd_reg();
}
