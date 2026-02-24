## 1. change the path to SCRIPT_DIR according to your directory tree
SCRIPT_DIR:={YOUR PATH}/NCP-MDAnalysis/15N_relaxation_rates/scripts/

## 2. specify the path to the directory where each trajectory's replicas were processed.
# The averaging is performed over individual tails H4-1 and H4-2 for set of trajectory's replicas.
# You can adapt this template to accommodate any number of replicas by simply adding a new variable named PATH_TO_PROCESSING_DIR_REPLICA_4 and etc.
PATH_TO_PROCESSING_DIR_REPLICA_1 = {YOUR PATH TO REPLICA-1}/H4_tails
PATH_TO_PROCESSING_DIR_REPLICA_2 = {YOUR PATH TO REPLICA-2}/H4_tails
PATH_TO_PROCESSING_DIR_REPLICA_3 =  {YOUR PATH TO REPLICA-3}/H4_tails
# PATH_TO_PROCESSING_DIR is a comma-separated sequence of directories for all replicas.
PATH_TO_PROCESSING_DIRS = ${PATH_TO_PROCESSING_DIR_REPLICA_1},\
						  ${PATH_TO_PROCESSING_DIR_REPLICA_2},\
						  ${PATH_TO_PROCESSING_DIR_REPLICA_3}

## 3. specify fit parameters
FIT_LIMIT_NS=1000 # ns
# you may specify logariphmic resampling of the correlation function (LAG_INDEX="log") and set the corresponding number of points from 0 ns to FIT_LIMIT_NS (N_LAG_POINTS).
# this step may prevent overfitting of the correlation function at large timescales
# If  you want to fit without logariphmic spacing, you should delete the these parameters (LAG_SPACING and N_LAG_POINTS)
LAG_SPACING="log"
N_LAG_POINTS=1000


## 4. specify experimental parameters needed for calculation
NMR_FREQ=1200e6 # Hz
TUMBLING_TIME=163.4 # ns