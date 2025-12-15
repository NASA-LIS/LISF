#!/bin/bash
# Run datetime configuration - Format: YYYYMMDDHHMM
RUN_DATETIME="$1"  # Read from first command line argument

# Configuration for single datetime job submission
TEMPLATE_FILE="SNIP_LDT_template.sh"
LOG_DIR="./log"
SUBMISSION_LOG="$LOG_DIR/job_submissions.log"

# Single datetime configuration - Format: YYYYMMDDHHMM
# RUN_DATETIME="202501200600"  # Example: Jan 20, 2025 at 06

# Create log directory if it doesn't exist
mkdir -p $LOG_DIR

# Function to parse datetime string
parse_datetime() {
    local datetime=$1
    
    # Extract components from YYYYMMDDHHMMformat
    local year=${datetime:0:4}
    local month=${datetime:4:2}
    local day=${datetime:6:2}
    local hour=${datetime:8:2}
    
    # Validate datetime format
    if [[ ${#datetime} -ne 12 ]]; then
        echo "ERROR: Invalid datetime format '$datetime'. Must be YYYYMMDDHHMM (12 digits)" | tee -a $SUBMISSION_LOG
        exit 1
    fi
    
    # Validate hour
    if [[ ! "$hour" =~ ^(00|06|12|18)$ ]]; then
        echo "ERROR: Invalid hour '$hour'. Must be 00, 06, 12, or 18" | tee -a $SUBMISSION_LOG
        exit 1
    fi
    
    # Return validated datetime
    echo "$datetime"
}

# Function to submit a job for the specified datetime
submit_job() {
    local datetime=$1
    
    # Create job-specific script in log directory
    local job_script="$LOG_DIR/job_${datetime}.sh"
    
    # Copy template and replace placeholder
    cp $TEMPLATE_FILE $job_script
    
    # Replace the single datetime placeholder
    sed -i "s/PLACEHOLDER_DATETIME/$datetime/g" $job_script
    
    # Update job name with datetime
    sed -i "s/#SBATCH --job-name=.*/#SBATCH --job-name=$datetime/g" $job_script
    
    # Update log file locations to use LOG_DIR
    sed -i "s|#SBATCH --output=\.\/log\/|#SBATCH --output=$LOG_DIR\/|g" $job_script
    sed -i "s|#SBATCH --error=\.\/log\/|#SBATCH --error=$LOG_DIR\/|g" $job_script
    
    # Submit the job
    echo "Submitting job for datetime: $datetime" | tee -a $SUBMISSION_LOG
    job_id=$(sbatch $job_script | awk '{print $4}')
    echo "Job submitted with ID: $job_id (Job name: $datetime)" | tee -a $SUBMISSION_LOG
    
    return 0
}

# Main execution
echo "Starting single datetime job submission at $(date)" | tee -a $SUBMISSION_LOG
echo "Template file: $TEMPLATE_FILE" | tee -a $SUBMISSION_LOG
echo "Target datetime: $RUN_DATETIME" | tee -a $SUBMISSION_LOG
echo "----------------------------------------" | tee -a $SUBMISSION_LOG

# Validate the datetime
echo "Validating datetime: $RUN_DATETIME" | tee -a $SUBMISSION_LOG
validated_datetime=$(parse_datetime $RUN_DATETIME)

echo "Validated datetime: $validated_datetime" | tee -a $SUBMISSION_LOG

# Submit the job
submit_job $validated_datetime

echo "" | tee -a $SUBMISSION_LOG
echo "Job submission completed at $(date)" | tee -a $SUBMISSION_LOG
echo "Check '$SUBMISSION_LOG' for submission details" | tee -a $SUBMISSION_LOG
echo "Job script created: $LOG_DIR/job_${validated_datetime}.sh" | tee -a $SUBMISSION_LOG

# Display queue status
echo "" | tee -a $SUBMISSION_LOG
echo "Current queue status:" | tee -a $SUBMISSION_LOG
squeue -u $(whoami) | tee -a $SUBMISSION_LOG
