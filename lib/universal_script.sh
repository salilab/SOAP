#!/bin/bash                         
#                                  
#$ -S /bin/bash                     
#$ -o /netapp/sali/gqdong/undefinedrundir/output
#$ -e /netapp/sali/gqdong/undefinedrundir/output
#$ -cwd                            
#$ -r y                            
#$ -j y                            
memundefined            
#$ -l arch=linux-x64
#$ -p 0
#$ -l h_rt=timeundefined               
#$ -t tundefined        

# Anything under here can be a csh script

#unfinishednumber

date
hostname


taskid=$SGE_TASK_ID
input=$taskid


#TMPDIR="/scrapp/gqdong/$JOB_ID/$taskid"
#mkdir -p $TMPDIR
#cd $TMPDIR

# Make a temporary directory on the scratch disk,
# specific to the user and SGE job.
cd /netapp/sali/gqdong/undefinedrundir/
export PYTHONPATH=$PYTHONPATH:SOAPPATH

# Copy input files to $TMPDIR here...

startsec=`date +%s`
# do the actual job (here some processing on 1abc and 2xyz)
/salilab/diva1/home/modeller/modpy.sh python2.6 /netapp/sali/gqdong/undefinedrundir/runme.py $input

endsec=`date +%s`

duration=$(($endsec - $startsec))

# Copy back output files from $TMPDIR here...

rm -rf "/scrapp/gqdong/$JOB_ID/$taskid"
rm -rf "/scratch/gqdong/$JOB_ID/$taskid"
rm -rf "/scrapp2/gqdong/$JOB_ID/$taskid"


touch /netapp/sali/gqdong/undefinedrundir/output/$taskid.$duration.`hostname`
date
exit 0


