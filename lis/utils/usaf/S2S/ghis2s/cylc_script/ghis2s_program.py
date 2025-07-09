#!/usr/bin/env python3

"""Integrates ghis2s S2S forecasting system into the GHI workflow."""

from copy import copy
from pathlib import Path
import logging
import os
import shutil
import subprocess
import sys
#from sharpy import RoseProgram

ENV_DEFINITION = {
    "CONFIG_FILE": (str, None),
    "FORECAST_YEAR": (int, None),
    "FORECAST_MONTH": (int, None),
    "USER_EMAIL": (str, None),
    "S2S_STEP": (str, "E2E"),
    "ONE_STEP": (bool, False),
    "SUBMIT_JOB": (bool, False),
    "E2ESDIR": (str, "/discover/nobackup/projects/ghilis/smahanam/Cylc_test/"),
    "OUTPUT_ROOT": (str, "/discover/nobackup/projects/ghilis/smahanam/GHI-repos/ghi-config/pkg/"),
    "MODEL": (str, "S2S"),
    "VARIANT": (str, ""),
    "PYTHONPATH": (str, None),
}

#class Ghis2sProgram(RoseProgram):
class Ghis2sProgram():
    """Main program class for ghis2s integration""" 
    def __init__(self, env_definition=None):
        #super().__init__(copy(env_definition or ENV_DEFINITION))
        # Temporary environment setup
        self.env = {}
        env_def = env_definition or ENV_DEFINITION
        for key, (type_func, default) in env_def.items():
            env_value = os.environ.get(key, default)
            if env_value is not None:
                if type_func == bool and isinstance(env_value, str):
                    self.env[key] = env_value.lower() in ('true', '1', 'yes', 'on')
                else:
                    self.env[key] = type_func(env_value)
            else:
                self.env[key] = default
        self._validate_required_env()
        self._setup_paths()
        
    def _validate_required_env(self):
        """Validate required environment variables are set."""
        required_vars = ["CONFIG_FILE", "FORECAST_YEAR", "FORECAST_MONTH", "USER_EMAIL"]
        missing_vars = []
        
        for var in required_vars:
            if self.env[var] is None:
                missing_vars.append(var)
                
        if missing_vars:
            logging.error("Missing required environment variables: %s", ", ".join(missing_vars))
            sys.exit(1)
            
    def _setup_paths(self):
        """Setup necessary paths for the workflow."""
        self._e2es_dir = Path(self.env["E2ESDIR"])
        self._cylc_home = Path(self.env["OUTPUT_ROOT"])
        self._log_dir = self._e2es_dir / "scratch" / f"{self.env['FORECAST_YEAR']:04d}{self.env['FORECAST_MONTH']:02d}"
        
    @property
    def model_variant(self):
        """Convenience property for joining the model and variant name."""
        return "_".join(filter(None, [self.env["MODEL"], self.env["VARIANT"]]))
        
    @property
    def forecast_date(self):
        """The forecast date in YYYY-MM format."""
        return f"{self.env['FORECAST_YEAR']:04d}-{self.env['FORECAST_MONTH']:02d}"
        
    @property
    def workflow_name(self):
        """The Cylc workflow name."""
        return f"S2S-{self.env['FORECAST_YEAR']:04d}{self.env['FORECAST_MONTH']:02d}"

    def _import_ghis2s(self):
        """Import ghis2s module with proper PYTHONPATH setup."""
        # Update PYTHONPATH if specified in environment
        if self.env["PYTHONPATH"]:
            current_pythonpath = os.environ.get("PYTHONPATH", "")
            if current_pythonpath:
                os.environ["PYTHONPATH"] = f"{self.env['PYTHONPATH']}:{current_pythonpath}"
            else:
                os.environ["PYTHONPATH"] = self.env["PYTHONPATH"]
                
            # Add to sys.path for immediate availability
            import sys
            for path in self.env["PYTHONPATH"].split(":"):
                if path and path not in sys.path:
                    sys.path.insert(0, path)
    
        try:
            from ghis2s.s2s_app.s2s_run import S2Srun
            return S2Srun
        except ImportError as err:
            logging.error("Failed to import ghis2s: %s", err)
            if self.env["PYTHONPATH"]:
                logging.error("PYTHONPATH was set to: %s", self.env["PYTHONPATH"])
            else:
                logging.error("PYTHONPATH not set. Please set PYTHONPATH environment variable to include ghis2s installation directory")
            sys.exit(1)
            
    def _run_s2s_workflow(self):
        """Execute the S2S workflow steps."""
        S2Srun = self._import_ghis2s()
        
        logging.info("Creating working directories and job files for S2S forecast initialized on %s-01", 
                    self.forecast_date)
        
        # Change to E2ES directory
        original_dir = os.getcwd()
        os.chdir(self._e2es_dir)
        
        try:
            # Initialize S2S run
            s2s = S2Srun(
                year=self.env["FORECAST_YEAR"],
                month=self.env["FORECAST_MONTH"],
                config_file=self.env["CONFIG_FILE"]
            )
            
            # Execute specified step or full workflow
            step = self.env["S2S_STEP"]
            one_step = self.env["ONE_STEP"]
            
            if step is not None:
                self._execute_s2s_step(s2s, step, one_step)
            else:
                s2s.main()
                
            # Write CYLC workflow runtime snippet
            s2s.write_cylc_snippet()
            
            # Submit SLURM jobs if requested
            if self.env["SUBMIT_JOB"]:
                s2s.submit_jobs()
                return
                
        finally:
            os.chdir(original_dir)
            
    def _execute_s2s_step(self, s2s, step, one_step):
        """Execute specific S2S workflow step."""
        step_map = {
            'E2E': lambda: s2s.main(),
            'LISDA': lambda: s2s.lis_darun(),
            'LDTICS': lambda: self._execute_ldtics_chain(s2s, one_step),
            'BCSD': lambda: self._execute_bcsd_chain(s2s, one_step),
            'FCST': lambda: self._execute_fcst_chain(s2s, one_step),
            'POST': lambda: self._execute_post_chain(s2s, one_step),
            'METRICS': lambda: self._execute_metrics_chain(s2s, one_step),
            'PLOTS': lambda: s2s.s2splots()
        }
        
        if step in step_map:
            step_map[step]()
        else:
            logging.error("Unknown S2S step: %s", step)
            sys.exit(1)
            
    def _execute_ldtics_chain(self, s2s, one_step):
        """Execute LDTICS and subsequent steps if not one_step."""
        s2s.ldt_ics()
        if not one_step:
            s2s.bcsd()
            s2s.lis_fcst()
            s2s.s2spost()
            s2s.s2smetric()
            s2s.s2splots()
            
    def _execute_bcsd_chain(self, s2s, one_step):
        """Execute BCSD and subsequent steps if not one_step."""
        s2s.bcsd()
        if not one_step:
            s2s.lis_fcst()
            s2s.s2spost()
            s2s.s2smetric()
            s2s.s2splots()
            
    def _execute_fcst_chain(self, s2s, one_step):
        """Execute FCST and subsequent steps if not one_step."""
        s2s.lis_fcst()
        if not one_step:
            s2s.s2spost()
            s2s.s2smetric()
            s2s.s2splots()
            
    def _execute_post_chain(self, s2s, one_step):
        """Execute POST and subsequent steps if not one_step."""
        s2s.s2spost()
        if not one_step:
            s2s.s2smetric()
            s2s.s2splots()
            
    def _execute_metrics_chain(self, s2s, one_step):
        """Execute METRICS and subsequent steps if not one_step."""
        s2s.s2smetric()
        if not one_step:
            s2s.s2splots()
            
    def _setup_cylc_workflow(self):
        """Setup the Cylc workflow directory and configuration."""
        # Create workflow directory
        workflow_dir = self._cylc_home / self.workflow_name
        workflow_dir.mkdir(exist_ok=True)
        os.chdir(workflow_dir)
        
        # Copy flow.cylc from generated CYLC_workflow.rc
        cylc_config_source = self._log_dir / "CYLC_workflow.rc"
        cylc_config_target = "flow.cylc" 
        
        if not cylc_config_source.exists():
            logging.error("CYLC workflow config not found: %s", cylc_config_source)
            sys.exit(1)
            
        shutil.copy(cylc_config_source, cylc_config_target)
        
        # Make shell scripts executable
        result = subprocess.run(
            f'chmod 755 {self._log_dir}/*/*sh', 
            shell=True, 
            capture_output=True, 
            text=True
        )
        if result.returncode != 0:
            logging.warning("Failed to chmod shell scripts: %s", result.stderr)
            
    def _update_cylc_config(self):
        """Update the Cylc configuration with user-specific values."""
        config_file = Path("flow.cylc")
        
        with config_file.open('r') as file:
            filedata = file.read()
            
        # Replace placeholders
        filedata = filedata.replace('USEREMAIL', self.env["USER_EMAIL"])
        filedata = filedata.replace('YYYY-MM', self.forecast_date)
        filedata = filedata.replace('--output ', '--output=')
        filedata = filedata.replace('--error ', '--error=')
        
        with config_file.open('w') as file:
            file.write(filedata)
            
    def _install_cylc_workflow(self):
        """Install the Cylc workflow."""
        command = f"cylc install --symlink-dirs=run={self._log_dir}"
        result = subprocess.run(command, shell=True)
        
        if result.returncode != 0:
            logging.error("Failed to install Cylc workflow")
            sys.exit(1)
            
    def _print_cylc_commands(self):
        """Print useful Cylc commands for the user."""
        print(f'Useful CYLC commands from {os.getcwd()}')
        print("   ")
        print(f'Run {self.workflow_name}: cylc play {self.workflow_name}')
        print(f'Monitor {self.workflow_name}: cylc tui {self.workflow_name}')
        print(f'Show status {self.workflow_name}: cylc show {self.workflow_name}')
        print(f'Stop {self.workflow_name}: cylc stop --now {self.workflow_name}')
        print(f'Cat log: cylc cat-log {self.workflow_name}')
        
    def run(self):
        """Run main program"""
        logging.info("Starting ghis2s S2S workflow integration")
        
        # Run S2S workflow
        self._run_s2s_workflow()
        
        # Skip Cylc setup if submitting jobs directly
        if self.env["SUBMIT_JOB"]:
            logging.info("Jobs submitted directly, skipping Cylc workflow setup")
            return
            
        # Setup Cylc workflow
        self._setup_cylc_workflow()
        self._update_cylc_config()
        self._install_cylc_workflow()
        
        # Print helpful commands
        self._print_cylc_commands()
        
        logging.info("ghis2s S2S workflow integration completed successfully")
        
    def main(self):
        """Temporary main method - TODO: Remove when sharpy available"""
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.run()
        
if __name__ == "__main__":
    Ghis2sProgram().main()
