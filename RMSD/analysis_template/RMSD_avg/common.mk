## 1. change the path to SCRIPT_DIR according to your directory tree
 SCRIPT_DIR:={YOUR PATH}/NCP-MDAnalysis/RMSD/scripts/

## 2. specify the path to the directory where each trajectory's replicas were processed.
# The averaging is performed over individual replicas.
# You can adapt this template to accommodate any number of replicas by simply adding a new variable named PATH_TO_PROCESSING_DIR_REPLICA_4 and etc.
PATH_TO_PROCESSING_DIR_REPLICA_1 = {YOUR PATH TO REPLICA-1}/H4_tails
PATH_TO_PROCESSING_DIR_REPLICA_2 = {YOUR PATH TO REPLICA-2}/H4_tails
PATH_TO_PROCESSING_DIR_REPLICA_3 =  {YOUR PATH TO REPLICA-3}/H4_tails
# PATH_TO_PROCESSING_DIR is a comma-separated sequence of directories for all replicas.
PATH_TO_PROCESSING_DIRS = ${PATH_TO_PROCESSING_DIR_REPLICA_1},\
						  ${PATH_TO_PROCESSING_DIR_REPLICA_2},\
						  ${PATH_TO_PROCESSING_DIR_REPLICA_3}
