/* typedef struct {
    double *restrict s0, *restrict ss, *restrict sd, *restrict mleaf;
} HRUState;

typedef struct {
    double *restrict sg, *restrict sr;
    HRUState hru[2];
} States;

typedef struct {
    const double *ne, *height, *hypsperc;
} Hypsometry;

typedef struct {
    const double * fhru;
    const double * hveg;
    const double * laimax;
} HRUSpatial;

typedef struct {
    const double * k_rout;
    const double * kssat;
    const double * prefr;
    const double * s0max;
    const double * slope;
    const double * ssmax;
    const double * k_gw;
    const double * kr_sd;
    const double * kr_0s;
    const double * k0sat;
    const double * sdmax;
    const double * kdsat;
    const double * ne;
} Spatial;

typedef struct {
    const double slope_coeff;
    const double pair;
    const double kr_coeff;
} Parameters;

typedef struct {
    const double * tat;
    const double * rgt;
    const double * pt;
    const double * avpt;
    const double * u2t;
    const double * radcskyt;
} Forcing;

typedef struct {
    const double alb_dry;
    const double alb_wet;
    const double cgsmax;
    const double er_frac_ref;
    const double fsoilemax;
    const double lairef;
    const double rd;
    const double s_sls;
    const double sla;
    const double tgrow;
    const double tsenc;
    const double ud0;
    const double us0;
    const double vc;
    const double w0lime;
    const double w0ref_alb;
    const double wdlimu;
    const double wslimu;
} HRUParameters;

typedef struct {
    double *restrict s0;
    double *restrict ss;
    double *restrict sd;
    double *restrict mleaf;
} HRUOutputs;

typedef struct {
    double *restrict e0;
    double *restrict etot;
    double *restrict dd;
    double *restrict s0;
    double *restrict ss;
    double *restrict sd;
    double *restrict qtot;
    double *restrict sr;
    double *restrict sg;
    HRUOutputs hru[2];
} Outputs;
*/

/* NOTE: Wendy Sharples Dec 2018: ALL the outputs need to be passed into the C program in order to use LIS io. 
NB: LIS passes in 1 cell at a time so needed to refactor based on this; cells is always 1 in LIS
*/

/* OLD declaration: void awral(Forcing inputs, Outputs outputs, States initial_states, States final_states, 
           Parameters params, Spatial spatial, Hypsometry hypso, HRUParameters *hruparams, HRUSpatial *hruspatial, int timesteps, int cells) {
*/

void awral_driver_600_(double *tat, double *rgt, double *pt, double *avpt, double *u2t, double *radcskyt, double *e0, double *etot, double *dd, double *s0_avg, double *ss_avg, double *sd_avg, double *qtot, double *sr, double *sg, double *s0, double * ss, double * sd, double * mleaf, double *slope_coeff, double *pair, double *kr_coeff, double *k_rout, double *kssat, double *prefr, double *s0max, double *slope, double *ssmax, double *k_gw, double *kr_sd, double *kr_0s, double *k0sat, double *sdmax, double *kdsat, double *ne, double * height, double * hypsperc, double * alb_dry, double * alb_wet, double * cgsmax, double * er_frac_ref, double * fsoilemax, double * lairef, double * rd, double * s_sls, double * sla, double * tgrow, double * tsenc, double * ud0, double * us0, double * vc, double * w0lime, double * w0ref_alb, double * wdlimu, double * wslimu, double * fhru, double * hveg, double * laimax, int *timesteps);


