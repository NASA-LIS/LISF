#!/bin/bash
WORKFLOW=$1

# Extract task names from flow.cylc
TASK_NAMES=$(grep -E '^\s*\[\[.*\]\]' ~/cylc-run/$WORKFLOW/runN/flow.cylc | sed 's/.*\[\[\(.*\)\]\].*/\1/' | grep -v '^[[:space:]]*$')

echo "#######################################################################"
echo "                          STATUS OF CYLC JOBS"
echo "#######################################################################"
echo ""
echo "            TASK NAME                    WALLTIME (HH:MM:SS)"
echo ""

# Variables to track overall start and end times
earliest_start=""
latest_end=""

for task in $TASK_NAMES; do
    for jobdir in ~/cylc-run/$WORKFLOW/runN/log/job/*/$task/NN/; do
        if [ -f "$jobdir/job.status" ]; then
            cycle=$(basename $(dirname $(dirname $jobdir)))
            submit_num=$(basename $jobdir)

            if [ "$task" != "log_monitor" ] && [ "$task" != "final_log_collect" ] && [ "$task" != "stop_log_monitor" ]; then
                # Extract timing info
                init_time=$(grep "CYLC_JOB_INIT_TIME=" "$jobdir/job.status" 2>/dev/null | cut -d'=' -f2)
                exit_time=$(grep "CYLC_JOB_EXIT_TIME=" "$jobdir/job.status" 2>/dev/null | cut -d'=' -f2)

                if [ -n "$init_time" ] && [ -n "$exit_time" ]; then
                    # Convert ISO timestamps to epoch seconds
                    init_epoch=$(date -d "$init_time" +%s 2>/dev/null)
                    exit_epoch=$(date -d "$exit_time" +%s 2>/dev/null)

                    if [ -n "$init_epoch" ] && [ -n "$exit_epoch" ]; then
                        # Track earliest start and latest end
                        if [ -z "$earliest_start" ] || [ "$init_epoch" -lt "$earliest_start" ]; then
                            earliest_start=$init_epoch
                        fi
                        if [ -z "$latest_end" ] || [ "$exit_epoch" -gt "$latest_end" ]; then
                            latest_end=$exit_epoch
                        fi
                        
                        walltime_seconds=$((exit_epoch - init_epoch))

                        # Convert to HH:MM:SS format
                        hours=$((walltime_seconds / 3600))
                        minutes=$(((walltime_seconds % 3600) / 60))
                        seconds=$((walltime_seconds % 60))

                        # Format task name (you might want to include cycle info)
                        task_display="${task}"

                        printf "%-30s                    %2dh %2dm %2ds\n" "$task_display" "$hours" "$minutes" "$seconds"
                    fi
                fi
            fi
        fi
    done
done

echo ""

# Calculate actual elapsed time (latest end - earliest start)
if [ -n "$earliest_start" ] && [ -n "$latest_end" ]; then
    total_elapsed_seconds=$((latest_end - earliest_start))
    
    total_hours=$((total_elapsed_seconds / 3600))
    total_days=$((total_hours / 24))
    remaining_hours=$((total_hours % 24))
    remaining_minutes=$(((total_elapsed_seconds % 3600) / 60))
    remaining_secs=$((total_elapsed_seconds % 60))
    
    printf "ELAPSED TIME : %2dd %2dh %2dm %2ds\n" "$total_days" "$remaining_hours" "$remaining_minutes" "$remaining_secs"
else
    echo "ELAPSED TIME : Unable to calculate"
fi

