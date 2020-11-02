BASEDIR=$(dirname $0)
cd $BASEDIR
date +"%d-%m-%Y, %H:%M" >> submission_output.log
python3 job_submission.py >> submission_output.log
