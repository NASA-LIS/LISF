#include "awral.h"
#include <tgmath.h>
#include <stdio.h>
#include <stdlib.h>

/* Overall comment: wrapping means that some changes to the code were necessary */
const int HYPS_LEN = 20; // +++ Could be supplied as input value
const int NUM_CELLS = 1; // +++ LIS can only take one tile at a time

void copycells(const double *in, double *out, int len) {
    int i;
    for (i = 0; i < len; i ++){
        out[i] = in[i];
    }
}

void build_hypso(double * sgtemp, double * height, double ne, int cells) {
    int c, l;
    for (c = 0; c < cells; c ++) {
        for (l = 0; l < 20; l++) {
            int idx = c*HYPS_LEN + l;
            sgtemp[idx] = height[idx] * ne * 1000.0;
        }
    }
}

double hyps_fsat(double * sgtemp, double * hypsfsat, double sg, double offset) {

    double sgmin = sgtemp[0] + offset;
    double sg_eff = sgmin + sg;
    if (sg + offset <= 0.) {
        return(0);
    } else if (sg_eff >= sgtemp[HYPS_LEN-1]) {
        return(1.0);
    } else {
        int id0 = -1;
        int i;
        for (i = 0; i < HYPS_LEN; i ++) {
            if (sg_eff >= sgtemp[i]) {
                id0++;
            } else {
                break;
            }
        }
        int id1 = id0+1;

        double srange = sgtemp[id1]-sgtemp[id0];
        double hdiff = hypsfsat[id1]-hypsfsat[id0];

        return (hypsfsat[id0]+hdiff*((sg_eff-sgtemp[id0])/srange));
    }

}

#define HRU_SUM(HRUVAL) hru_sum((double *)HRUVAL,fhru_hru,hidx);

double hru_sum(double *hru_data,double fhru[2],int hidx[2]) {
    return (fhru[0]*hru_data[hidx[0]]+fhru[1]*hru_data[hidx[1]]);
}

void calc_soil_flows(double *s, double *i, double *e, double *drain, double *iflow, const double smax, const double ksat, const double rh, const double km) {
    if ((*s + *i) <= *e) {
        *e = *s + *i;
        *s = 0.0;
        *drain = 0.0;
        *iflow = 0.0;
    } else {
        double a = ksat / (smax*smax);
        const double b = 1.0;
        double c = -(*s + *i - *e);
        if ((smax - *s + ksat) <= (*i - *e)) {
            *drain = (1 - rh) * ksat;
            *iflow = (rh * ksat) + (*s - smax - ksat + *i - *e);
            *s = smax;
        } else {
            *s = (-b + (pow(b * b - 4. * a * c,0.5))) / (2. * a);
            if (*s < smax) {
                double pfull = (*s/smax);
                *drain = (1 - rh) * ksat * pfull * pfull;
                *iflow = rh * ksat * pfull * pfull;
            } else {
                *drain = (1 - rh) * ksat;
                *iflow = (rh * ksat) + (-c - smax -ksat);
                *s = smax;
            }
        }
    }
}

/* New declaration for LIS - NOTE that LIS passes in one cell at at time- need cell number for building hypsometric curves */
void awral_driver_600_(double *tat, double *rgt, double *pt, double *avpt, double *u2t, double *radcskyt, double *e0, double *etot, double *dd, double *s0_avg, double *ss_avg, double *sd_avg, double *qtot, double *sr, double *sg, double *s0, double * ss, double * sd, double * mleaf, double *slope_coeff, double *pair, double *kr_coeff, double *k_rout, double *kssat, double *prefr, double *s0max, double *slope, double *ssmax, double *k_gw, double *kr_sd, double *kr_0s, double *k0sat, double *sdmax, double *kdsat, double *ne, double * height, double * hypsperc, double * alb_dry, double * alb_wet, double * cgsmax, double * er_frac_ref, double * fsoilemax, double * lairef, double * rd, double * s_sls, double * sla, double * tgrow, double * tsenc, double * ud0, double * us0, double * vc, double * w0lime, double * w0ref_alb, double * wdlimu, double * wslimu, double * fhru, double * hveg, double * laimax, int *timesteps) {    

    // HRU States
    // 1D arrays
    double *s0_ = malloc(2*sizeof(double));
    double *ss_ = malloc(2*sizeof(double));
    double *sd_ = malloc(2*sizeof(double));
    double *mleaf_ = malloc(2*sizeof(double));

    // HRU Values exported to non-HRU regions

    double *dd_ = malloc(2*sizeof(double));
    double *eg_ = malloc(2*sizeof(double));
    double *y_ = malloc(2*sizeof(double));
    double *qif_ = malloc(2*sizeof(double));
    double *qr_ = malloc(2*sizeof(double));

    int hru, c, ts;
    #pragma ivdep
    for (hru=0; hru<2; hru++){
        s0_[hru] = s0[hru];
        ss_[hru] = ss[hru];
        sd_[hru] = sd[hru];
        mleaf_[hru] = mleaf[hru];
    }

    const double stefbolz  = 5.67e-8;

    double *sgtemp_ = malloc(HYPS_LEN*sizeof(double));

    build_hypso(sgtemp_,height,*ne,NUM_CELLS);

    double *eff_rd = malloc(2*sizeof(double));//[2][cells];

    for (hru=0; hru<2; hru++){
        eff_rd[hru] = rd[hru] * 1000.0 * *ne;
    }

    double fsat_;//[cells];
    double fegt_;//[cells];

    // Load scalar constants:
    double slope_coeff_const = *slope_coeff;
    double pair_const = *pair;
    double kr_coeff_const = *kr_coeff;
    
    for (ts=0; ts<*timesteps; ts++) {

        // Hypso
        for (c=0; c<NUM_CELLS; c++) {
            fsat_ = hyps_fsat(&sgtemp_[c*HYPS_LEN],hypsperc,*sg,0.0);
            int idx = NUM_CELLS*ts + c;

            //Zero fill fhru-averaged outputs (NB if this was more than 1 cell - need to pass in as arrays)
            *e0 = 0.0;
            *etot = 0.0;
            *dd = 0.0;
            *s0_avg = 0.0;
            *ss_avg = 0.0;
            *sd_avg = 0.0;
        }

        // HRU loop
        for (hru=0; hru<2; hru++) {
            for (c=0; c<NUM_CELLS; c++) {
                fegt_ = hyps_fsat(&sgtemp_[c*HYPS_LEN],hypsperc,*sg,eff_rd[hru]);
            }

            double alb_dry_hru = alb_dry[hru];
            double alb_wet_hru = alb_wet[hru];
            double cgsmax_hru = cgsmax[hru];
            double er_frac_ref_hru = er_frac_ref[hru];
            double fsoilemax_hru = fsoilemax[hru];
            double lairef_hru = lairef[hru];
            double rd_hru = rd[hru];
            double s_sls_hru = s_sls[hru];
            double sla_hru = sla[hru];
            double tgrow_hru = tgrow[hru];
            double tsenc_hru = tsenc[hru];
            double ud0_hru = ud0[hru];
            double us0_hru = us0[hru];
            double vc_hru = vc[hru];
            double w0lime_hru = w0lime[hru];
            double w0ref_alb_hru = w0ref_alb[hru];
            double wdlimu_hru = wdlimu[hru];
            double wslimu_hru = wslimu[hru];

            #pragma ivdep
            #pragma vector always
            for (c=0; c<NUM_CELLS; c++) {
                int idx = NUM_CELLS*ts + c;
                int hru_cidx = hru*NUM_CELLS+c;

                double s0_c = s0_[hru];
                double ss_c = ss_[hru];
                double sd_c = sd_[hru];
                double mleaf_c = mleaf_[hru];

                //Load spatial data
                double k_rout_c = *k_rout;
                double slope_c = *slope;
                double kssat_c = *kssat;
                double k_gw_c = *k_gw;
                double sdmax_c = *sdmax;
                double prefr_c = *prefr;
                double kr_0s_c = *kr_0s;
                double s0max_c = *s0max;
                double kdsat_c = *kdsat;
                double kr_sd_c = *kr_sd;
                double k0sat_c = *k0sat;
                double ssmax_c = *ssmax;
                double ne_c = *ne;

                //Load HRU spatial data
                double fhru_c = fhru[hru];
                double hveg_c = hveg[hru];
                double laimax_c = laimax[hru];
                
                //Load forcing data
                double rgt_idx = *rgt;
                double tat_idx = *tat;
                double avpt_idx = *avpt;
                double u2t_idx = *u2t;
                double pt_idx = *pt;
                double radcskyt_idx = *radcskyt;

                double ta = tat_idx;
                double pg = pt_idx;
                double rg = rgt_idx;
                double pe = avpt_idx;
                double u2 = u2t_idx;
                double radclearsky = radcskyt_idx;
                
                double lai = sla_hru * mleaf_c;
                double fveg = 1.0 - exp(-lai / lairef_hru);
                double fsoil = 1.0 - fveg;

                double w0 = s0_c / s0max_c;
                double ws = ss_c / ssmax_c;
                double wd = sd_c / sdmax_c;

                double fsat = fsat_;
                double fegt = fegt_;

                double pes = 610.8 * exp(17.27 * ta / (237.3 + ta));
                double frh = pe / pes;

                double keps = 1.4e-3 * ((ta / 187.0) * (ta / 187.0) + ta / 107.0 + 1.0) * (6.36 * pair_const + pe) / pes;

                double delta = 4217.457 / ((240.97 + ta)*(240.97 + ta)) * pes;
                double gamma = 0.000646 * pair_const * (1.0 + 0.000946 * ta);
                double lambda = 2.501 - (0.002361 * ta);

                double rgeff = rg;

                double alb_veg = 0.452 * vc_hru;
                double alb_soil = alb_wet_hru + (alb_dry_hru - alb_wet_hru) * exp(-w0 / w0ref_alb_hru);
                double alb = fveg * alb_veg + fsoil * alb_soil;
                double rsn = (1.0 - alb) * rgeff;

                double tkelv = ta + 273.15;
                double rlin = stefbolz * pow(tkelv,4.0) * (1.0 - (1.0 - 0.65 * pow((pe / tkelv),0.14)) * (1.35 * rgeff / radclearsky - 0.35));
                double rlout = 1.0 * stefbolz * pow(tkelv,4);
                double rln = (rlin - rlout) * 0.0864;
                double rneff = rsn + rln;

                double fh = log(813.0 / hveg_c - 5.45);
                double ku2 = 0.305 / (fh * (fh + 2.3));
                double ga = ku2*u2;

                //  Potential evaporation
                //New for v4 (different to original in Van Dijk (2010)) 
                double e0_c = (delta * rneff + 6.43 / 1000.0 * gamma * (pes - pe) * (1.0 + 0.536 * u2)) / (delta + gamma) / lambda;

                e0_c = fmax(e0_c,0.);
                double usmax = us0_hru*fmin(ws/wslimu_hru,1.0);
                double udmax = ud0_hru*fmin(wd/wdlimu_hru,1.0);
                double u0 = fmax(usmax,udmax);

                double gsmax = cgsmax_hru * vc_hru;
                double gs = fveg * gsmax;
                double ft = 1.0 / (1.0 + (keps / (1.0 + keps)) * ga / gs);
                double etmax = ft*e0_c;

                double et = fmin(u0,etmax);

                double us, ud;

                if (u0 > 0.0) {
                    us = fmin( (usmax/(usmax+udmax)) * et,ss_c-0.01);
                    ud = fmin( (udmax/(usmax+udmax)) * et,sd_c-0.01);
                } else {
                    us = 0.0;
                    ud = 0.0; 
                }

                et = us + ud;

                // Remaining ET
                double et_rem = e0_c-et;

                double fsoile = fsoilemax_hru * fmin(w0/w0lime_hru,1.);
                double es = (1.0 - fsat)*fsoile*(et_rem);
                double eg_c = fsat*fsoilemax_hru*(et_rem);

                double sveg = s_sls_hru * lai;
                double fer = er_frac_ref_hru * fveg;
                double pwet = -log(1.0 - fer/fveg) * sveg / fer;

                // +++
                // Need to limit this and update remaining ET
                double ei;

                if (pg<pwet) {
                    ei = fveg * pg;
                } else {
                    ei = fveg * pwet + fer * (pg-pwet);
                }

                double pn = fmax(0.0,pg-ei);
                double rhof = (1.0 - fsat)*(pn - (prefr_c*tanh(pn/prefr_c)));
                double rsof = fsat * pn;
                double qr_c = rhof + rsof;
                double i = pn - qr_c;

                //+++ Purely spatial, should be in spatial_mapping
                double slope_eff = slope_coeff_const * slope_c;

                double rh_0s = tanh(slope_eff*w0)*tanh(kr_coeff_const*(kr_0s_c-1.0)*w0);
                double rh_sd = tanh(slope_eff*ws)*tanh(kr_coeff_const*(kr_sd_c-1.0)*ws);
                double rh_dd = 0.;

                double km_s0 = pow((k0sat_c * kssat_c),0.5);
                double km_ss = pow((kssat_c * kdsat_c),0.5);
                double km_sd = kdsat_c;

                double d0, if0, ds, ifs, dd_c, ifd;

                calc_soil_flows(&s0_c,&i,&es,&d0,&if0,s0max_c,k0sat_c, rh_0s, km_s0);
                calc_soil_flows(&ss_c,&d0,&us,&ds,&ifs,ssmax_c,kssat_c, rh_sd, km_ss);
                calc_soil_flows(&sd_c,&ds,&ud,&dd_c,&ifd,sdmax_c,kdsat_c, rh_dd, km_sd);

                et = us + ud;

                double wz = fmax(0.01, sd_c/sdmax_c);

                double y_c;

                // +++
                // Need to limit this and update remaining ET
                if ((fegt-fsat) > 0.0) {
                    y_c = (fegt-fsat) * fsoilemax_hru * et_rem;
                } else {
                    y_c = 0.0;
                }

                double fvmax = 1. - exp(-fmax(laimax_c,0.0027777778450399637)/lairef_hru);

                double fveq = (1.0 / fmax((e0_c/u0)-1.0,1e-3)) * (keps/(1.0+keps))*(ga/gsmax);
                fveq = fmin(fveq,fvmax);

                double dmleaf = -log(1.0 - fveq) * lairef_hru / sla_hru - mleaf_c;
                double mleafnet = 0.0;

                if (dmleaf > 0.) {
                    mleafnet = dmleaf/tgrow_hru;
                } else if (dmleaf < 0.) {
                    mleafnet = dmleaf/tsenc_hru;
                }
                
                mleaf_c += mleafnet;

                s0_[hru] = s0_c;
                ss_[hru] = ss_c;
                sd_[hru] = sd_c;
                mleaf_[hru] = mleaf_c;

                dd_[hru] = dd_c;
                eg_[hru] = eg_c;
                y_[hru] = y_c;
                qif_[hru] = if0 + ifs + ifd;
                qr_[hru] = qr_c;

                double ee = es + eg_c + ei + y_c;
                double etot_c = ee + et;

                //Hardcoded HRU specific outputs
                s0[hru] = s0_c;
                ss[hru] = ss_c;
                sd[hru] = sd_c;
                mleaf[hru] = mleaf_c;

                //Hardcoded combined (FHRU-weighted) outputs
                *e0 += e0_c*fhru_c;
                *etot += etot_c*fhru_c;
                *dd += dd_c*fhru_c;
                *s0_avg += s0_c*fhru_c;
                *ss_avg += ss_c*fhru_c;
                *sd_avg += sd_c*fhru_c;
            } // end cell loop
        } // end HRU loop

        // Post HRU
        #pragma ivdep
        #pragma vector always
        for (c=0; c<NUM_CELLS; c++) {
            double k_rout_c = *k_rout;
            double slope_c = *slope;
            double kssat_c = *kssat;
            double k_gw_c = *k_gw;
            double sdmax_c = *sdmax;
            double prefr_c = *prefr;
            double kr_0s_c = *kr_0s;
            double s0max_c = *s0max;
            double kdsat_c = *kdsat;
            double kr_sd_c = *kr_sd;
            double k0sat_c = *k0sat;
            double ssmax_c = *ssmax;
            double ne_c = *ne;

            int hidx[2];
            hidx[0] = c;
            hidx[1] = NUM_CELLS + c;

            double fhru_hru[2];
            fhru_hru[0] = fhru[0];
            fhru_hru[1] = fhru[1];

            double dd_c = HRU_SUM(dd_);
            double eg_c = HRU_SUM(eg_);
            double y_c = HRU_SUM(y_);
            double qif_c = HRU_SUM(qif_);
            double qr_c = HRU_SUM(qr_);
            
	    // set from initial states
            double sg_c = *sg;
            double sr_c = *sr;

            sg_c = sg_c + dd_c;
            double qg = (1.0 - exp(-k_gw_c)) * fmax(0.0,sg_c);
            sg_c = sg_c - qg - eg_c - y_c;

            // calculate groundwater levels as AHD

            double gw_level = height[c*HYPS_LEN] + sg_c/(ne_c * 1000.0);

            sr_c = fmax(0.0,(sr_c + qr_c + qg + qif_c));

            double kd = (1.0 - exp(-k_rout_c));

            double qtot_c = fmin(sr_c,(sr_c * kd));
            sr_c = sr_c - qtot_c;

            int idx = NUM_CELLS*ts + c;

            //Cell level hardcoded outputs
            *qtot = qtot_c;
            // Set final states
            *sr  = sr_c;
            *sg = sg_c;
        }
    }

    free(s0_);
    free(ss_);
    free(sd_);
    free(mleaf_);
    free(dd_);
    free(eg_);
    free(y_);
    free(qif_);
    free(qr_);

    free(sgtemp_);
    free(eff_rd);

}


