''' 
Python main script for hindcast processing
'''
import os
import sys
import argparse
from ghis2s.main.hcast_module import S2SHcast
from ghis2s.main import walltime
from ghis2s.shared import logging_utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
    parser.add_argument('-y', '--year', required=False, type=int, help='forecast year')
    parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
    parser.add_argument('-p', '--preprocess', required=False, default=None, type=str,
                        help='preprocess options: regrid, clim, stats')
    parser.add_argument('-j', '--submit_job', action='store_true',
                        help='Submit SLURM jobs (default: False)?')
    parser.add_argument('-r', '--report', action='store_true',
                        help='Print report')
    parser.add_argument('-l', '--logging', action='store_true',
                        help='Write centralized log file')

    args = parser.parse_args()
    if  args.preprocess is not None:
        s2s = S2SHcast(month=args.month, config_file=args.config_file, pre_process=True)
    else:
        s2s = S2SHcast(month=args.month, config_file=args.config_file, year=args.year)

    # Print job report
    if  args.report:
        JOB_SCHEDULE = os.path.join(s2s.scrdir, 'SLURM_JOB_SCHEDULE')
        if os.path.exists(JOB_SCHEDULE):
            walltime.slurm(s2s.scrdir)
        else:
            if args.step is None:
                CYLC_WORKFLOW = f'cylc_e2e_{s2s.yyyy}{s2s.mm}'
            else:
                CYLC_WORKFLOW = f'cylc_{args.step.lower()}_{s2s.yyyy}{s2s.mm}'
            walltime.cylc(CYLC_WORKFLOW)
        sys.exit()

    # Write LOG file
    if args.logging:
        logging_utils.write_centralized_logging(s2s.scrdir)
        sys.exit()

    if  args.preprocess is not None:
        if args.preprocess.lower() == 'regrid':
            s2s.bcsd(clim=True)
        elif args.preprocess.lower() == 'clim':
            s2s.clim_metforce()
            s2s.clim_nmme()
            s2s.clim_nafpa()
            s2s.write_cylc_snippet()
        elif args.preprocess.lower() == 'stats':
            s2s.s2smetric()
        else:
            print(f"Invalid preprocessor {args.preprocess}")
            sys.exit()
    else:
        s2s.bcsd()
        s2s.lis_fcst()
        s2s.s2spost()
        s2s.write_cylc_snippet()

    # save S2S schedule dictionary
    logging_utils.save_schedule(s2s.scrdir, s2s.schedule)

    # Submit SLURM jobs
    # -----------------
    if  args.submit_job:
        s2s.submit_jobs()
