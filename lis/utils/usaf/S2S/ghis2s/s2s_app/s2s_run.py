''' 
GHI-S2S run script
'''
import os
import sys
import argparse
from ghis2s.main.experiment_setup import S2Srun
from ghis2s.main import walltime
from ghis2s.shared import logging_utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
    parser.add_argument('-y', '--year', required=True, type=int, help='forecast year')
    parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
    parser.add_argument('-s', '--step', required=False, default=None, type=str,
                        help='S2S step: LISDA, LDTICS, BCSD, FCST, POST, METRICS or PLOTS')
    parser.add_argument('-o', '--one_step', action='store_true',
                        help='Is only one step (default: False)?')
    parser.add_argument('-j', '--submit_job', action='store_true',
                        help='Submit SLURM jobs (default: False)?')
    parser.add_argument('-r', '--report', action='store_true',
                        help='Print report')
    parser.add_argument('-l', '--logging', action='store_true',
                        help='Write centralized log file')
    parser.add_argument('-d', '--delete_forecast', action='store_true',
                        help='delete forecast')

    args = parser.parse_args()

    s2s = S2Srun(year=args.year, month=args.month, config_file=args.config_file)

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

    # Delete forecast
    if args.delete_forecast:
        s2s.delete_forecast()
        sys.exit()

    if args.step is not None:
        if args.step == 'LISDA':
            s2s.lis_darun()
        elif args.step == 'LDTICS':
            s2s.ldt_ics()
            if not args.one_step:
                s2s.bcsd()
                s2s.lis_fcst()
                s2s.s2spost()
                s2s.s2smetric()
                s2s.s2splots()
        elif args.step == 'BCSD':
            s2s.nmme_file_checker()
            s2s.clim_files_checker()
            s2s.cfsv2_file_checker()
            s2s.bcsd()
            if not args.one_step:
                s2s.lis_fcst()
                s2s.s2spost()
                s2s.s2smetric()
                s2s.s2splots()
        elif args.step == 'FCST':
            s2s.lis_fcst()
            if not args.one_step:
                s2s.s2spost()
                s2s.s2smetric()
                s2s.s2splots()
        elif args.step == 'POST':
            s2s.s2spost()
            if not args.one_step:
                s2s.s2smetric()
                s2s.s2splots()
        elif args.step == 'METRICS':
            s2s.s2smetric()
            if not args.one_step:
                s2s.s2splots()
        elif args.step == 'PLOTS':
            s2s.s2splots()

        # Write CYLC workflow runtime snippet
        # -----------------------------------
        s2s.write_cylc_snippet()
    else:
        s2s.main()

    # save S2S schedule dictionary
    logging_utils.save_schedule(s2s.scrdir, s2s.schedule)

    # Submit SLURM jobs
    # -----------------
    if  args.submit_job:
        s2s.write_log_monitoring_script()
        s2s.submit_jobs()
