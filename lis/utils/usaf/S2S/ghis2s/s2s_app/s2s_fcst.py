import os
import sys
import yaml
import argparse

#######################################################################################
#                 End-to_end Forecast                         
#  python s2s_run.py -y year -m month -c config_file 
#######################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', required=True, type=str, help='config file')
parser.add_argument('-y', '--year', required=True, type=int, help='forecast year')
parser.add_argument('-m', '--month', required=True, type=int, help='forecast month')
parser.add_argument('-r', '--report', required=False, type=str, help='Report')
parser.add_argument('-s', '--step', required=False, type=str, help='S2S step: LISDA, LDTICS, BCSD, FCST, POST, METRICS or PLOTS')
parser.add_argument('-o', '--one_step', required=False, type=str, help='Is only one step?')
args = parser.parse_args()

with open(args.config_file, 'r', encoding="utf-8") as file:
    config = yaml.safe_load(file)

# Import ghis2s
# -------------
sys.path.append(config['SETUP']['LISFDIR'] + 'lis/utils/usaf/S2S/')
from ghis2s.shared.s2s_run import S2Srun

s2s = S2Srun(year=args.year, month=args.month, config_file=args.config_file)
s2s.main()
