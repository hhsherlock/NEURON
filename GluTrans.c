/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__GluTrans
#define _nrn_initial _nrn_initial__GluTrans
#define nrn_cur _nrn_cur__GluTrans
#define _nrn_current _nrn_current__GluTrans
#define nrn_jacob _nrn_jacob__GluTrans
#define nrn_state _nrn_state__GluTrans
#define _net_receive _net_receive__GluTrans 
#define kstates kstates__GluTrans 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define Gluout _p[0]
#define density _p[1]
#define itrans _p[2]
#define itransLog _p[3]
#define C1 _p[4]
#define C2 _p[5]
#define C3 _p[6]
#define C4 _p[7]
#define C5 _p[8]
#define C6 _p[9]
#define surf _p[10]
#define volin _p[11]
#define volout _p[12]
#define ko _p[13]
#define ki _p[14]
#define DC1 _p[15]
#define DC2 _p[16]
#define DC3 _p[17]
#define DC4 _p[18]
#define DC5 _p[19]
#define DC6 _p[20]
#define v _p[21]
#define _g _p[22]
#define _ion_ko	*_ppvar[0]._pval
#define _ion_ki	*_ppvar[1]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_u(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_GluTrans", _hoc_setdata,
 "u_GluTrans", _hoc_u,
 0, 0
};
#define u u_GluTrans
 extern double u( _threadargsprotocomma_ double , double );
 /* declare global and static user variables */
#define Gluin Gluin_GluTrans
 double Gluin = 0.3;
#define Kout Kout_GluTrans
 double Kout = 3;
#define Kin Kin_GluTrans
 double Kin = 120;
#define Naout Naout_GluTrans
 double Naout = 150;
#define Nain Nain_GluTrans
 double Nain = 15;
#define charge charge_GluTrans
 double charge = 1.6e-019;
#define k61 k61_GluTrans
 double k61 = 0.0002;
#define k16 k16_GluTrans
 double k16 = 0.0016;
#define k65 k65_GluTrans
 double k65 = 0.1;
#define k56 k56_GluTrans
 double k56 = 1;
#define k54 k54_GluTrans
 double k54 = 10;
#define k45 k45_GluTrans
 double k45 = 4;
#define k43 k43_GluTrans
 double k43 = 0.6;
#define k34 k34_GluTrans
 double k34 = 0.2;
#define k32 k32_GluTrans
 double k32 = 0.5;
#define k23 k23_GluTrans
 double k23 = 0.015;
#define k21 k21_GluTrans
 double k21 = 0.1;
#define k12 k12_GluTrans
 double k12 = 20;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "k12_GluTrans", "l",
 "k21_GluTrans", "/ms",
 "k23_GluTrans", "l",
 "k32_GluTrans", "/ms",
 "k34_GluTrans", "/ms",
 "k43_GluTrans", "/ms",
 "k45_GluTrans", "/ms",
 "k54_GluTrans", "l",
 "k56_GluTrans", "/ms",
 "k65_GluTrans", "l",
 "k16_GluTrans", "l",
 "k61_GluTrans", "l",
 "Nain_GluTrans", "mM/l",
 "Naout_GluTrans", "mM/l",
 "Kin_GluTrans", "mM/l",
 "Kout_GluTrans", "mM/l",
 "Gluin_GluTrans", "mM/l",
 "charge_GluTrans", "coulombs",
 "Gluout_GluTrans", "mM/l",
 "C1_GluTrans", "/cm2",
 "C2_GluTrans", "/cm2",
 "C3_GluTrans", "/cm2",
 "C4_GluTrans", "/cm2",
 "C5_GluTrans", "/cm2",
 "C6_GluTrans", "/cm2",
 "itrans_GluTrans", "mA/cm2",
 0,0
};
 static double C60 = 0;
 static double C50 = 0;
 static double C40 = 0;
 static double C30 = 0;
 static double C20 = 0;
 static double C10 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "k12_GluTrans", &k12_GluTrans,
 "k21_GluTrans", &k21_GluTrans,
 "k23_GluTrans", &k23_GluTrans,
 "k32_GluTrans", &k32_GluTrans,
 "k34_GluTrans", &k34_GluTrans,
 "k43_GluTrans", &k43_GluTrans,
 "k45_GluTrans", &k45_GluTrans,
 "k54_GluTrans", &k54_GluTrans,
 "k56_GluTrans", &k56_GluTrans,
 "k65_GluTrans", &k65_GluTrans,
 "k16_GluTrans", &k16_GluTrans,
 "k61_GluTrans", &k61_GluTrans,
 "Nain_GluTrans", &Nain_GluTrans,
 "Naout_GluTrans", &Naout_GluTrans,
 "Kin_GluTrans", &Kin_GluTrans,
 "Kout_GluTrans", &Kout_GluTrans,
 "Gluin_GluTrans", &Gluin_GluTrans,
 "charge_GluTrans", &charge_GluTrans,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GluTrans",
 "Gluout_GluTrans",
 "density_GluTrans",
 0,
 "itrans_GluTrans",
 "itransLog_GluTrans",
 0,
 "C1_GluTrans",
 "C2_GluTrans",
 "C3_GluTrans",
 "C4_GluTrans",
 "C5_GluTrans",
 "C6_GluTrans",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 23, _prop);
 	/*initialize range parameters*/
 	Gluout = 2e-005;
 	density = 1e+012;
 	_prop->param = _p;
 	_prop->param_size = 23;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[2]; /* ko */
 	_ppvar[1]._pval = &prop_ion->param[1]; /* ki */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _GluTrans_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 3);
  _extcall_thread = (Datum*)ecalloc(2, sizeof(Datum));
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 23, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GluTrans C:/Users/savtc/Dropbox/Task2Astro/GluTrans.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double F = 96485.3;
 static double PI = 3.14159;
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 extern double *_nrn_thread_getelm();
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[6], _dlist1[6]; static double *_temp1;
 static int kstates();
 
static int kstates (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<6;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 /* COMPARTMENT volin {
     }
   */
 /* COMPARTMENT volout {
     }
   */
 /* ~ C1 <-> C2 ( Gluout * k12 * u ( _threadargscomma_ v , - 0.1 ) , k21 )*/
 f_flux =  Gluout * k12 * u ( _threadargscomma_ v , - 0.1 ) * C1 ;
 b_flux =  k21 * C2 ;
 _RHS1( 5) -= (f_flux - b_flux);
 _RHS1( 4) += (f_flux - b_flux);
 
 _term =  Gluout * k12 * u ( _threadargscomma_ v , - 0.1 ) ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 4 ,5)  -= _term;
 _term =  k21 ;
 _MATELM1( 5 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ C2 <-> C3 ( Naout * k23 * u ( _threadargscomma_ v , 0.5 ) , k32 )*/
 f_flux =  Naout * k23 * u ( _threadargscomma_ v , 0.5 ) * C2 ;
 b_flux =  k32 * C3 ;
 _RHS1( 4) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  Naout * k23 * u ( _threadargscomma_ v , 0.5 ) ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 3 ,4)  -= _term;
 _term =  k32 ;
 _MATELM1( 4 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ C3 <-> C4 ( k34 * u ( _threadargscomma_ v , 0.4 ) , k43 )*/
 f_flux =  k34 * u ( _threadargscomma_ v , 0.4 ) * C3 ;
 b_flux =  k43 * C4 ;
 _RHS1( 3) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  k34 * u ( _threadargscomma_ v , 0.4 ) ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  k43 ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ C4 <-> C5 ( k45 , k54 * Gluin )*/
 f_flux =  k45 * C4 ;
 b_flux =  k54 * Gluin * C5 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  k45 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  k54 * Gluin ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ C5 <-> C6 ( k56 * u ( _threadargscomma_ v , 0.6 ) , k65 * Nain )*/
 f_flux =  k56 * u ( _threadargscomma_ v , 0.6 ) * C5 ;
 b_flux =  k65 * Nain * C6 ;
 _RHS1( 1) -= (f_flux - b_flux);
 
 _term =  k56 * u ( _threadargscomma_ v , 0.6 ) ;
 _MATELM1( 1 ,1)  += _term;
 _term =  k65 * Nain ;
 _MATELM1( 1 ,0)  -= _term;
 /*REACTION*/
  /* ~ C6 <-> C1 ( ki * k61 , k16 * u ( _threadargscomma_ v , 0.6 ) * ko )*/
 f_flux =  ki * k61 * C6 ;
 b_flux =  k16 * u ( _threadargscomma_ v , 0.6 ) * ko * C1 ;
 _RHS1( 5) += (f_flux - b_flux);
 
 _term =  ki * k61 ;
 _MATELM1( 5 ,0)  -= _term;
 _term =  k16 * u ( _threadargscomma_ v , 0.6 ) * ko ;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
   /* C1 + C2 + C3 + C4 + C5 + C6 = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= C6 ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= C5 ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= C4 ;
 _MATELM1(0, 3) = 1;
 _RHS1(0) -= C3 ;
 _MATELM1(0, 4) = 1;
 _RHS1(0) -= C2 ;
 _MATELM1(0, 5) = 1;
 _RHS1(0) -= C1 ;
 /*CONSERVATION*/
   } return _reset;
 }
 
double u ( _threadargsprotocomma_ double _lx , double _lth ) {
   double _lu;
 _lu = exp ( _lth * _lx / ( 2.0 * ( 26.7 ) ) ) ;
   
return _lu;
 }
 
static void _hoc_u(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  u ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<6;_i++) _p[_dlist1[_i]] = 0.0;}
 /* COMPARTMENT volin {
   }
 */
 /* COMPARTMENT volout {
   }
 */
 /* ~ C1 <-> C2 ( Gluout * k12 * u ( _threadargscomma_ v , - 0.1 ) , k21 )*/
 f_flux =  Gluout * k12 * u ( _threadargscomma_ v , - 0.1 ) * C1 ;
 b_flux =  k21 * C2 ;
 DC1 -= (f_flux - b_flux);
 DC2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C2 <-> C3 ( Naout * k23 * u ( _threadargscomma_ v , 0.5 ) , k32 )*/
 f_flux =  Naout * k23 * u ( _threadargscomma_ v , 0.5 ) * C2 ;
 b_flux =  k32 * C3 ;
 DC2 -= (f_flux - b_flux);
 DC3 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C3 <-> C4 ( k34 * u ( _threadargscomma_ v , 0.4 ) , k43 )*/
 f_flux =  k34 * u ( _threadargscomma_ v , 0.4 ) * C3 ;
 b_flux =  k43 * C4 ;
 DC3 -= (f_flux - b_flux);
 DC4 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C4 <-> C5 ( k45 , k54 * Gluin )*/
 f_flux =  k45 * C4 ;
 b_flux =  k54 * Gluin * C5 ;
 DC4 -= (f_flux - b_flux);
 DC5 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C5 <-> C6 ( k56 * u ( _threadargscomma_ v , 0.6 ) , k65 * Nain )*/
 f_flux =  k56 * u ( _threadargscomma_ v , 0.6 ) * C5 ;
 b_flux =  k65 * Nain * C6 ;
 DC5 -= (f_flux - b_flux);
 DC6 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ C6 <-> C1 ( ki * k61 , k16 * u ( _threadargscomma_ v , 0.6 ) * ko )*/
 f_flux =  ki * k61 * C6 ;
 b_flux =  k16 * u ( _threadargscomma_ v , 0.6 ) * ko * C1 ;
 DC6 -= (f_flux - b_flux);
 DC1 += (f_flux - b_flux);
 
 /*REACTION*/
   /* C1 + C2 + C3 + C4 + C5 + C6 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<6;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 /* COMPARTMENT volin {
 }
 */
 /* COMPARTMENT volout {
 }
 */
 /* ~ C1 <-> C2 ( Gluout * k12 * u ( _threadargscomma_ v , - 0.1 ) , k21 )*/
 _term =  Gluout * k12 * u ( _threadargscomma_ v , - 0.1 ) ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 4 ,5)  -= _term;
 _term =  k21 ;
 _MATELM1( 5 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ C2 <-> C3 ( Naout * k23 * u ( _threadargscomma_ v , 0.5 ) , k32 )*/
 _term =  Naout * k23 * u ( _threadargscomma_ v , 0.5 ) ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 3 ,4)  -= _term;
 _term =  k32 ;
 _MATELM1( 4 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ C3 <-> C4 ( k34 * u ( _threadargscomma_ v , 0.4 ) , k43 )*/
 _term =  k34 * u ( _threadargscomma_ v , 0.4 ) ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  k43 ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ C4 <-> C5 ( k45 , k54 * Gluin )*/
 _term =  k45 ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  k54 * Gluin ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ C5 <-> C6 ( k56 * u ( _threadargscomma_ v , 0.6 ) , k65 * Nain )*/
 _term =  k56 * u ( _threadargscomma_ v , 0.6 ) ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  k65 * Nain ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ C6 <-> C1 ( ki * k61 , k16 * u ( _threadargscomma_ v , 0.6 ) * ko )*/
 _term =  ki * k61 ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 5 ,0)  -= _term;
 _term =  k16 * u ( _threadargscomma_ v , 0.6 ) * ko ;
 _MATELM1( 0 ,5)  -= _term;
 _MATELM1( 5 ,5)  += _term;
 /*REACTION*/
   /* C1 + C2 + C3 + C4 + C5 + C6 = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 6;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ko = _ion_ko;
  ki = _ion_ki;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 6; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 6, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ko = _ion_ko;
  ki = _ion_ki;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 2);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  C6 = C60;
  C5 = C50;
  C4 = C40;
  C3 = C30;
  C2 = C20;
  C1 = C10;
 {
   C1 = 0.9074 ;
   C2 = 0.0199 ;
   C3 = 0.0435 ;
   C4 = 0.0103 ;
   C5 = 0.0142 ;
   C6 = 0.0047 ;
   volin = 1.0 ;
   volout = 1.0 ;
   surf = 1.0 ;
   ko = 3.0 ;
   ki = 110.0 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ko = _ion_ko;
  ki = _ion_ki;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   itrans = - charge * density * ( 1e+006 ) * ( 0.6 * ( C1 * k16 * ko * u ( _threadargscomma_ v , 0.6 ) - C6 * k61 * ki ) - 0.1 * ( C1 * k12 * Gluout * u ( _threadargscomma_ v , - 0.1 ) - C2 * k21 ) + 0.5 * ( C2 * k23 * Naout * u ( _threadargscomma_ v , 0.5 ) - C3 * k32 ) + 0.4 * ( C3 * k34 - C4 * k43 ) + 0.6 * ( C5 * k56 * u ( _threadargscomma_ v , 0.6 ) - C6 * k65 * Nain ) ) ;
   }
 _current += itrans;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ko = _ion_ko;
  ki = _ion_ki;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ko = _ion_ko;
  ki = _ion_ki;
 {  sparse_thread(&_thread[_spth1]._pvoid, 6, _slist1, _dlist1, _p, &t, dt, kstates, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 6; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 }}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(C6) - _p;  _dlist1[0] = &(DC6) - _p;
 _slist1[1] = &(C5) - _p;  _dlist1[1] = &(DC5) - _p;
 _slist1[2] = &(C4) - _p;  _dlist1[2] = &(DC4) - _p;
 _slist1[3] = &(C3) - _p;  _dlist1[3] = &(DC3) - _p;
 _slist1[4] = &(C2) - _p;  _dlist1[4] = &(DC2) - _p;
 _slist1[5] = &(C1) - _p;  _dlist1[5] = &(DC1) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "GluTrans.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "\n"
  "The model of Glutamate  transporter.\n"
  "is based on two papers, \n"
  "\n"
  "from the paper \n"
  "\n"
  "1. Zhang Z1, Tao Z, Gameiro A, Barcelona S, Braams S, Rauen T, Grewer C. \n"
  "Transport direction determines the kinetics of substrate transport by the glutamate transporter EAAC1.\n"
  "Proc Natl Acad Sci U S A. 2007 Nov 13;104(46):18025-30. Epub 2007 Nov 8.\n"
  "\n"
  "we determine the basic kinetic scheme for glutamate transporters, \n"
  "\n"
  "from the  paper\n"
  " \n"
  "2. Bergles, D.E. & Jahr, C.E. \n"
  "Synaptic activation of glutamate transporters in hippocampal astrocytes. Neuron 19, 1297-1308 (1997).\n"
  "\n"
  "we corrected the numerical values of the kinetic constants corresponding to the dynamics of glutamate transporters in astrocytes\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX  GluTrans\n"
  "	USEION k READ ko, ki\n"
  "    RANGE part, C1, C2, C3, C4, C5, C6\n"
  "    GLOBAL k12, k21, k23, k32, k34, k43, k45, k54, k56, k65, k16, k61\n"
  "    GLOBAL Nain, Naout, Kin, Kout, Gluin, charge \n"
  "    RANGE  itrans, Gluout, density, itransLog\n"
  "    NONSPECIFIC_CURRENT itrans\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (l) = (liter)\n"
  "    (nA) = (nanoamp)\n"
  "    (mV) = (millivolt)\n"
  "    (mA) = (milliamp)\n"
  "    (pS) = (picosiemens)\n"
  "    (umho) = (micromho)\n"
  "    (mM) = (milli/liter)\n"
  "    (uM) = (micro/liter)\n"
  "    F = (faraday) (coulombs)\n"
  "        PI      = (pi)       (1)\n"
  "}\n"
  "\n"
  "PARAMETER {	\n"
  "    : Rates\n"
  "\n"
  "    k12 = 20           (l /mM /ms)\n"
  "    k21 = 0.1          (/ms)\n"
  "    k23 = 0.015       (l /mM /ms)\n"
  "    k32 = 0.5          (/ms)\n"
  "    k34 = 0.2          (/ms)\n"
  "    k43 = 0.6          (/ms)\n"
  "    k45 = 4            (/ms)\n"
  "    k54 = 10           (l /mM /ms)\n"
  "    k56 = 1            (/ms) \n"
  "    k65 = 0.1          (l /mM /ms) \n"
  "    k16 = 0.0016          (l /mM /ms)\n"
  "    k61 =  2e-4        (l /mM /ms)\n"
  "\n"
  "    Nain = 15        (mM/l)\n"
  "    Naout = 150   (mM/l)\n"
  "    Kin = 120       (mM/l)\n"
  "    Kout = 3        (mM/l)\n"
  "    Gluin = 0.3      (mM/l)\n"
  "    Gluout = 20e-6	(mM/l)\n"
  "\n"
  "    density =1e12  : (/cm2) : 10000 per um2\n"
  "    charge = 1.6e-19 (coulombs)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    v	   (mV)		:  voltage\n"
  "    itrans (mA/cm2)            : \n"
  "    surf   (cm2)\n"
  "    volin  (liter)\n"
  "    volout (liter)\n"
  "    itransLog\n"
  "	ko      (mM)\n"
  "	ki      (mM)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    : Transporter  states (all fractions)\n"
  "            : \n"
  "    C1	(/cm2)	:  \n"
  "    C2	(/cm2)	:  \n"
  "    C3	(/cm2)	: \n"
  "    C4	(/cm2)	: \n"
  "    C5	(/cm2)	: \n"
  "    C6  (/cm2)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    C1= 0.9074    \n"
  "    C2= 0.0199    \n"
  "    C3= 0.0435    \n"
  "    C4= 0.0103    \n"
  "    C5= 0.0142    \n"
  "    C6= 0.0047\n"
  "    volin = 1\n"
  "    volout = 1\n"
  "    surf = 1\n"
  "	ko = 3      \n"
  "	ki = 110 \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE kstates METHOD sparse\n"
  "    \n"
  "	itrans=-charge*density*(1e+006)*(0.6*(C1*k16*ko*u(v,0.6)-C6*k61*ki) -0.1*(C1*k12*Gluout*u(v,-0.1)-C2*k21)+0.5*(C2*k23*Naout*u(v,0.5)-C3*k32)+0.4*( C3*k34-C4*k43)+0.6*(C5*k56*u(v,0.6)-C6*k65*Nain) )\n"
  "	\n"
  "    : itrans=-charge*density*(1e+006)*(0.6*(C1*k16*Kout*u(v,0.6)-C6*k61*Kin) -0.1*(C1*k12*Gluout*u(v,-0.1)-C2*k21)+0.5*(C2*k23*Naout*u(v,0.5)-C3*k32)+0.4*( C3*k34-C4*k43)+0.6*(C5*k56*u(v,0.6)-C6*k65*Nain) )\n"
  "    : itransLog=log(-itrans*(1e+006))\n"
  "\n"
  "    :itrans=-charge*density*(1e+006)*(0.6*(C1*k16*Kout*u(v,0.6)-C6*k61*Kin) +0.4*( C3*k34-C4*k43)+0.6*(C5*k56*u(v,0.6)-C6*k65*Nain) )	  \n"
  "}\n"
  "\n"
  "KINETIC kstates {\n"
  "            COMPARTMENT volin { Nain Kin Gluin}\n"
  "            COMPARTMENT volout { Naout Kout Gluout}\n"
  "            : COMPARTMENT surf { C1 C2 C3 C4 C5 C6}\n"
  "        : surf=1 : !!!!!!!\n"
  "        ~ C1   <-> C2      (Gluout*k12*u(v,-0.1), k21)\n"
  "        ~ C2  <-> C3       (Naout*k23*u(v,0.5),k32)\n"
  "        ~ C3 <-> C4	       (k34*u(v,0.4),k43)\n"
  "        ~ C4 <-> C5 	   (k45,k54*Gluin)\n"
  "        ~ C5 <-> C6	       (k56*u(v,0.6),k65*Nain)\n"
  "        ~ C6  <-> C1       (ki*k61, k16*u(v,0.6)*ko)\n"
  "        \n"
  "    CONSERVE C1+C2+C3+C4+C5+C6= 1\n"
  "}\n"
  "\n"
  "FUNCTION u(x(mV), th) {\n"
  "    u = exp(th*x/(2*(26.7 (mV))))\n"
  "}\n"
  ;
#endif
