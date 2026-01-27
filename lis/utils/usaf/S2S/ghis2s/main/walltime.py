#!/usr/bin/env python3
''' The module contains walltime computing functions '''
import os
import re
import glob
import subprocess
from datetime import datetime

def cylc(workflow):
    ''' compute walltimes of individual tasks using flow.cylc and job.status files '''
    flow_cylc_path = f"~/cylc-run/{workflow}/runN/flow.cylc"
    flow_cylc_path = os.path.expanduser(flow_cylc_path)
    task_names = []
    try:
        with open(flow_cylc_path, 'r', encoding="utf-8") as f:
            for line in f:
                match = re.match(r'^\s*\[\[(.+)\]\]', line.strip())
                if match:
                    task_name = match.group(1)
                    if task_name.strip():
                        task_names.append(task_name)
    except FileNotFoundError:
        print(f"Error: Could not find flow.cylc at {flow_cylc_path}")
        return

    print("#" * 71)
    print("                          STATUS OF CYLC JOBS")
    print("#" * 71)
    print()
    print("            TASK NAME                    WALLTIME (HH:MM:SS)")
    print()

    earliest_start = None
    latest_end = None
    excluded_tasks = {"log_monitor", "final_log_collect", "stop_log_monitor"}

    for task in task_names:
        if task in excluded_tasks:
            continue

        job_pattern = os.path.expanduser(f"~/cylc-run/{workflow}/runN/log/job/*/{task}/NN/")
        job_dirs = glob.glob(job_pattern)

        for jobdir in job_dirs:
            job_status_path = os.path.join(jobdir, "job.status")
            if os.path.isfile(job_status_path):
                init_time = None
                exit_time = None

                try:
                    with open(job_status_path, 'r', encoding="utf-8") as f:
                        for line in f:
                            if line.startswith("CYLC_JOB_INIT_TIME="):
                                init_time = line.split('=', 1)[1].strip()
                            elif line.startswith("CYLC_JOB_EXIT_TIME="):
                                exit_time = line.split('=', 1)[1].strip()
                except:
                    continue

                if init_time and exit_time:
                    try:
                        init_dt = datetime.fromisoformat(init_time.replace('Z', '+00:00'))
                        exit_dt = datetime.fromisoformat(exit_time.replace('Z', '+00:00'))
                        init_epoch = int(init_dt.timestamp())
                        exit_epoch = int(exit_dt.timestamp())

                        if earliest_start is None or init_epoch < earliest_start:
                            earliest_start = init_epoch
                        if latest_end is None or exit_epoch > latest_end:
                            latest_end = exit_epoch

                        walltime_seconds = exit_epoch - init_epoch

                        hours = walltime_seconds // 3600
                        minutes = (walltime_seconds % 3600) // 60
                        seconds = walltime_seconds % 60

                        # Format and print task info
                        task_display = task
                        print(f"{task_display:<30}                    {hours:2d}h {minutes:2d}m {seconds:2d}s")

                    except (ValueError, TypeError):
                        continue

    print()

    # Calculate actual elapsed time (latest end - earliest start)
    if earliest_start is not None and latest_end is not None:
        total_elapsed_seconds = latest_end - earliest_start
        total_hours = total_elapsed_seconds // 3600
        total_days = total_hours // 24
        remaining_hours = total_hours % 24
        remaining_minutes = (total_elapsed_seconds % 3600) // 60
        remaining_secs = total_elapsed_seconds % 60

        print(f"ELAPSED TIME : {total_days:2d}d {remaining_hours:2d}h {remaining_minutes:2d}m {remaining_secs:2d}s")
    else:
        print("ELAPSED TIME : Unable to calculate")

def slurm(scrdir):
    ''' SLURM walltimes '''

    print("#" * 71)
    print("                          STATUS OF SLURM JOBS")
    print("#" * 71)
    print(" ")
    print("            JOB FILE                 WALLTIME (HH:MM:SS)")
    print(" ")

    def datetime_to_seconds(dt_str):
        """Convert datetime string to seconds since epoch"""
        try:
            for fmt in ['%Y-%m-%dT%H:%M:%S', '%Y-%m-%d %H:%M:%S', '%Y-%m-%dT%H:%M:%S.%f']:
                try:
                    dt = datetime.strptime(dt_str, fmt)
                    return int(dt.timestamp())
                except ValueError:
                    continue
            # If none of the formats work, try using date command as fallback
            result = subprocess.run(['date', '-d', dt_str, '+%s'],
                                  capture_output=True, text=True, check=True)
            return int(result.stdout.strip())
        except (subprocess.CalledProcessError, ValueError):
            return None

    def seconds_to_datetime(seconds):
        """Convert seconds since epoch to datetime string"""
        try:
            result = subprocess.run(['date', '-d', f'@{seconds}', '+%Y-%m-%dT%H:%M:%S'],
                                  capture_output=True, text=True, check=True)
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            return None

    def compute_elapse(job_id):
        """Compute elapsed time for a job"""
        try:
            # Run sacct command
            result = subprocess.run(['sacct', '-j', str(job_id), '--format=Start,End', '-P'],
                                  capture_output=True, text=True, check=True)

            lines = result.stdout.strip().split('\n')[1:]  # Skip header
            if not lines:
                return None

            start_times = []
            end_times = []

            for line in lines:
                parts = line.split('|')
                if len(parts) >= 2:
                    start_times.append(parts[0])
                    end_times.append(parts[1])

            if not start_times or not end_times:
                return None

            min_start = None
            max_end = None

            for start in start_times:
                start_sec = datetime_to_seconds(start)
                if start_sec is not None:
                    if min_start is None or start_sec < min_start:
                        min_start = start_sec

            for end in end_times:
                end_sec = datetime_to_seconds(end)
                if end_sec is not None:
                    if max_end is None or end_sec > max_end:
                        max_end = end_sec

            if min_start is None or max_end is None:
                return None

            elapsed_seconds = max_end - min_start
            min_start_datetime = seconds_to_datetime(min_start)
            max_end_datetime = seconds_to_datetime(max_end)

            return f"{min_start_datetime}|{max_end_datetime}|{elapsed_seconds}"

        except subprocess.CalledProcessError:
            return None

    # Read job IDs and job files from SLURM_JOB_SCHEDULE
    schedule_file = os.path.join(scrdir, 'SLURM_JOB_SCHEDULE')

    jobids = []
    jobfiles = []
    try:
        with open(schedule_file, 'r', encoding="utf-8") as f:
            for line in f:
                if '.j' in line:
                    # Equivalent to tr -s ' ' | cut -d' ' -f1,2
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        jobids.append(parts[0])
                        jobfiles.append(parts[1])
    except FileNotFoundError:
        print(f"Error: Could not find {schedule_file}")
        return

    if not jobids:
        print("No jobs found in schedule file")
        return

    tlen = len(jobids)
    if jobfiles and jobfiles[-1] != 'set_permission.j':
        pass  # Keep tlen as is
    else:
        tlen -= 1

    start_time = None
    end_job = None

    for cjobs, (jobid, jobfile) in enumerate(zip(jobids, jobfiles)):
        if jobfile != 'set_permission.j':
            try:
                # Get job times using sacct
                result = subprocess.run(['sacct', '-j', jobid, '--format=start,end,elapsed'],
                                      capture_output=True, text=True, check=True)

                lines = result.stdout.strip().split('\n')
                if len(lines) > 1:
                    times_line = lines[-1]  # tail -1 equivalent
                    times_parts = times_line.split()

                    if len(times_parts) >= 2:
                        start_job = times_parts[0]
                        end_job = times_parts[1]

                        if (end_job not in ['Unknown', 'None'] and
                            start_job not in ['Unknown', 'None']):

                            result_str = compute_elapse(jobid)
                            if result_str:
                                start_job, end_job, elapsed_seconds_str = result_str.split('|')
                                elapsed_seconds = int(elapsed_seconds_str)

                                hours = elapsed_seconds // 3600
                                minutes = (elapsed_seconds % 3600) // 60
                                seconds = elapsed_seconds % 60

                                ehms = f"{hours}h {minutes}m {seconds}s"
                                job_number = f"{cjobs + 1}/{tlen}"
                                print(f"{job_number:>7s} {jobfile:<36s} {ehms}")

                                if cjobs == 0:
                                    start_time = start_job

                            else:
                                break
                        else:
                            break
                    else:
                        break
                else:
                    break

            except subprocess.CalledProcessError:
                break

    # Calculate total elapsed time
    if start_time and end_job:
        min_start = datetime_to_seconds(start_time)
        max_end = datetime_to_seconds(end_job)

        if min_start and max_end:
            elapsed_seconds = max_end - min_start

            days = elapsed_seconds // 86400
            hours = (elapsed_seconds % 86400) // 3600
            minutes = (elapsed_seconds % 3600) // 60
            seconds = elapsed_seconds % 60

            ehms = f"{days}d {hours}h {minutes}m {seconds}s"
            print(" ")
            print(f"ELAPSED TIME : {ehms}")
