
"""
Logging utilities for GHIS2S with 16WS compatibility.
"""

import os
import glob
import subprocess
import json
from pathlib import Path
from datetime import datetime
import re

# 16WS compatible formats
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
LOG_FORMAT = "%(asctime)s [%(levelname)s]: %(message)s"

class TaskLogger:
    """Logger for individual S2S tasks with 16WS compatibility"""

    def __init__(self, task_name, log_dir, description):
        """
        Initialize task logger
        
        Parameters:
        -----------
        task_name : str
            Name of the task (with or without .j extension)
        """
        # Clean task name (remove .j if present)
        self.task_name = task_name.replace('.j', '') if task_name.endswith('.j') else task_name
        pid = os.getpid()
        self.log_file = f"{log_dir}/logs/{self.task_name}_{pid}.log"

        # Create file and write header
        with open(self.log_file, 'w', encoding="utf-8") as f:
            f.write("[HEADER]\n")
            f.write(f"Script: {self.task_name}\n")
            if "\n" in description:
                desc_lines = description.split("\n")
                f.write(f"Description: {desc_lines[0]}\n")
                for line in desc_lines[1:]:
                    f.write(f"             {line}\n")
            else:
                f.write(f"Description: {description}\n")
            f.write(f"Process ID: {pid}\n")
            f.write("\n")

    def _write_log(self, level, message, subtask=None):
        """Write a log entry with 16WS log formatting"""
        timestamp = datetime.now().strftime(DATE_FORMAT)
        subtask_tag = f"[{subtask}] " if subtask else ""

        with open(self.log_file, 'a', encoding="utf-8") as f:
            f.write(f"[{timestamp}] [{level}] {subtask_tag}{message}\n")

    def info(self, message, subtask=None):
        """Log info level message"""
        self._write_log("INFO", message, subtask)

    def warning(self, message, subtask=None):
        """Log warning level message"""
        self._write_log("WARNING", message, subtask)

    def error(self, message, subtask=None):
        """Log error level message"""
        self._write_log("ERROR", message, subtask)

    def debug(self, message, subtask=None):
        """Log debug level message"""
        self._write_log("DEBUG", message, subtask)

class GHIS2SLogger:
    """
    Centralized logger for GHIS2S tasks
    Reads individual task logs and combines them into a centralized log
    """

    def __init__(self, scratch_path):
        """
        Initialize centralized logger

        Parameters:
        -----------
        scratch_path : str or Path
            Path to scratch directory containing task logs and schedule file
        """
        self.scratch_path = Path(scratch_path)
        self.schedule_path = self.scratch_path / 'ghis2s_schedule.json'
        self.central_log_path = self.scratch_path / 'ghis2s_main.log'

    def load_schedule(self):
        """Load schedule dictionary from disk"""
        try:
            if self.schedule_path.exists():
                with open(self.schedule_path, 'r', encoding="utf-8") as f:
                    return json.load(f)
            return {}
        except Exception as e:
            print(f"Error loading schedule from {self.schedule_path}: {e}")
            return {}

    def group_log_by_subtask(self, log_content):
        """
        Group log entries by subtask while maintaining chronological order within groups
        
        Returns organized content with:
        1. Header section (unchanged)
        2. General entries (no subtask)
        3. Grouped subtask entries
        4. Final general entries
        """

        mem_patterns = ["out-of-memory", "oom-kill", "Memory limit exceeded", "cgroup.*memory"]
        time_patterns = ["TIME LIMIT", "TIMEOUT", "CANCELLED.*TIME"]
        def check_error_patterns(file_pattern, patterns):
            ''' checks if SLURM standard error message contains the pattern '''
            for pattern in patterns:
                result = subprocess.run(f'grep -i "{pattern}" {file_pattern}',
                                        shell=True, capture_output=True,
                                        text=True, check=False)
                if result.returncode == 0:
                    return True
            return False

        lines = log_content.split('\n')

        header_lines = []
        general_entries = []
        subtask_groups = {}

        current_section = 'header'
        mem_error = False
        time_error = False
        current_task = ''

        for line in lines:
            if current_section == 'header' and re.match(r'\[20\d{2}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\]'
                                                        , line):
                current_section = 'body'

            if current_section == 'header':
                if mem_error:
                    timestamp = datetime.now().strftime(DATE_FORMAT)
                    general_entries.append(
                        f"[{timestamp}] [ERROR] {current_task}run killed "
                        f"by the cgroup out-of-memory handler")
                    mem_error = False
                if time_error:
                    timestamp = datetime.now().strftime(DATE_FORMAT)
                    general_entries.append(
                        f"[{timestamp}] [ERROR] {current_task}run CANCELLED DUE TO TIME LIMIT")
                    time_error = False
                header_lines.append(line)
                if line.startswith("Script:"):
                    current_task = line[7:].strip().replace('run', '')
                    mem_error = check_error_patterns(
                        f'{self.scratch_path}/*/logs/{current_task}*.err', mem_patterns)
                    time_error = check_error_patterns(
                        f'{self.scratch_path}/*/logs/{current_task}*.err', time_patterns)
            else:
                # Check for subtask pattern: [timestamp] [LEVEL] [subtask] message
                subtask_match = re.search(
                    r'\[20\d{2}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\] '
                    r'\[(INFO|WARNING|ERROR|DEBUG)\] '
                    r'\[([^\]]+)\] (.+)', line)
                if subtask_match:
                    timestamp_level = line[:line.find('] [', line.find('] [') + 1) + 1]
                    subtask = subtask_match.group(2)
                    message = subtask_match.group(3)

                    if subtask not in subtask_groups:
                        subtask_groups[subtask] = []
                    subtask_groups[subtask].append(f"{timestamp_level} {message}")
                else:
                    if line.strip():
                        general_entries.append(line)

        # Reconstruct the log after grouping to subtasks
        result_lines = []

        # 1. Add header
        result_lines.extend(header_lines)
        result_lines.append('')

        # 2. Add initial general entries (like "Starting parallel processing")
        if general_entries:
            result_lines.extend(general_entries[:1])
            result_lines.append('')

        # 3. Add grouped subtask entries
        if subtask_groups:
            for subtask in sorted(subtask_groups.keys()):
                result_lines.append(f"--- {subtask} Processing ---")
                result_lines.extend(subtask_groups[subtask])
                result_lines.append('')

        # 4. Add remaining general entries
        if len(general_entries) > 1:
            result_lines.extend(general_entries[1:])

        return '\n'.join(result_lines)

    def update_centralized_log(self):
        """
        Create/overwrite centralized log by copying all task logs entirely with optional grouping
        """
        log_files = self.find_log_files()

        if not log_files:
            print("No log files found to process")
            return 0

        with open(self.central_log_path, 'w', encoding="utf-8") as central:
            # Write main header with creation time
            timestamp = datetime.now().strftime(DATE_FORMAT)
            central.write(f"=== GHIS2S Centralized Log - Created {timestamp} ===\n\n")

            # Copy each log file entirely in schedule order
            for log_file_path in log_files:
                log_file = Path(log_file_path)

                try:
                    # Add separator for each log file
                    central.write(f"\n{'='*80}\n")
                    central.write(f"LOG FILE: {log_file.name}\n")
                    central.write(f"{'='*80}\n\n")

                    # Read and group log file content
                    with open(log_file, 'r', encoding="utf-8") as task_log:
                        content = task_log.read()

                    # Apply subtask grouping
                    grouped_content = self.group_log_by_subtask(content)
                    central.write(grouped_content)

                    # Add spacing between log files
                    central.write("\n\n")

                except Exception as e:
                    central.write(f"ERROR: Could not read log file {log_file}: {e}\n\n")

        print(f"Centralized log created with {len(log_files)} task logs")
        return len(log_files)

    def find_log_files(self):
        """
        Find all task log files based on schedule order

        Returns:
        --------
        list : ordered list of log file paths based on schedule order
        """
        log_files = []
        schedule = self.load_schedule()

        if schedule:
            # Process in schedule order (maintains order as written)
            for task_key, task_info in schedule.items():
                task_name = task_key.replace('.j', '') if task_key.endswith('.j') else task_key
                subdir = task_info.get('subdir', '')

                task_dir = self.scratch_path / subdir
                log_pattern = str(task_dir / f"logs/{task_name}*.log")
                matching_logs = glob.glob(log_pattern)

                # Add found log files in the schedule order
                for log_file in matching_logs:
                    if Path(log_file).exists() and Path(log_file).stat().st_size > 0:
                        log_files.append(log_file)

        return log_files

def save_schedule(scratch_path, schedule_dict):
    """
    Save schedule dictionary to in SCRDIR
    
    Parameters:
    -----------
    scratch_path : str
        Path to scratch directory
    schedule_dict : dict
        Schedule dictionary to save
    """
    try:
        schedule_path = Path(scratch_path) / 'ghis2s_schedule.json'
        with open(schedule_path, 'w', encoding="utf-8") as f:
            json.dump(schedule_dict, f, indent=2)
    except Exception as e:
        print(f"Error saving schedule: {e}")

def write_centralized_logging(scratch_path):
    """
    Update centralized log file with entries from all task logs
    
    Parameters:
    -----------
    scratch_path : str or Path
        Path to scratch directory with task logs
    """
    logger = GHIS2SLogger(scratch_path)
    entry_count = logger.update_centralized_log()
    print(f"Updated centralized log with {entry_count} new entries")
    return f"Updated centralized log with {entry_count} new entries"
