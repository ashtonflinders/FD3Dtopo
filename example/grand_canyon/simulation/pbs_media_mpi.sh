#!/bin/bash
# modified from PSSClabs scripts

### Number of nodes n nodes using m cpus Processor Per Node
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00
### Queue name
##PBS -q high
#PBS -q default
### Job name
#PBS -N g.canyon.media

### Merge stderr with stdout
#PBS -j oe
### Mail to user
#PBS -m eb
### Declare job-non-rerunable
#PBS -r n

### generate real time output in user home directory
##PBS -k oe

PBS_PWD="`pwd`";
cd "${PBS_O_WORKDIR}";
#THIS_HOST="`hostname`";
THIS_HOST=$(hostname --short)
DATE=$(date)

#MPICH_ROOT="/opt/mpich/p4-intel";
#MPICH_BUILD="p4.shm-intel10";
#MPICH_ROOT="/opt/mpich1/${MPICH_BUILD}";
#MPI_BUILD=""
#MPI_ROOT=/opt/openmpi/${MPI_BUILD}
MPI_BUILD="3.1-intel"
MPI_ROOT=/opt/mpich3/${MPI_BUILD}

# This job's working directory
echo Job ID: $PBS_JOBID
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
echo Running on host ${THIS_HOST}
echo Time is ${DATE}
echo Directory is $(pwd)
echo This job runs on the following processors:
echo $(cat $PBS_NODEFILE)

# How many processors per node?
PPN=$(grep -c ${THIS_HOST} $PBS_NODEFILE | tr -d '[:blank;]')
echo This job has allocated ${PPN} processors per node.

#Define number of processors
NPROCS=$(wc -l < $PBS_NODEFILE)
echo This job has allocated $NPROCS processors.

# How many nodes?
NHOSTS=$(( ${NPROCS} / ${PPN} ))
echo This job has allocated ${NHOSTS} nodes.


#MPIRUN_BIN="${MPICH_ROOT}/bin/mpirun";
#MPIEXEC_BIN="/opt/mpiexec/bin/mpiexec";
#CLEANIPCS_BIN="${MPICH_ROOT}/sbin/cleanipcs";
##FNM_BIN="./bin/seis3d_media";
#FNM_BIN="./bin/seis3d_media_mpi";
##FNM_BIN="./bin.WithQs/seis3d_media_mpi";
MPIEXEC_BIN="${MPI_ROOT}/bin/mpiexec -n ${NPROCS} -f ${PBS_NODEFILE}"

FNM_BIN="./bin/seis3d_media_mpi";

#NPROCS="`wc -l < ${PBS_NODEFILE} | tr -d '[:blank:]'`";

#P4_GLOBMEMSIZE=$(( 1024 * 1024 * 1024 * 7 ))

hr()
{
	perl -e 'print "\n" . "-"x70 . "\n\n"';
}

#print_job_info()
#{
#	hr;
#	printf "Torque Job ID: %s\n" "${PBS_JOBID}";
#	printf "\nRunning on host %s @ %s\n" "${THIS_HOST}" "`date`";
#	printf "\nStarting directory was %s\n" "${PBS_PWD}";
#	printf "Working directory is %s\n" "${PBS_O_WORKDIR}";
#	printf "The PWD is %s\n" "`pwd`";
#	printf "\nThis job runs on the following processors:\n\n\t";
#	printf "%s " `cat ${PBS_NODEFILE} | sort`;
#	printf "\n\n";
#	printf "This job has allocated %s nodes/processors.\n" "${NPROCS}";
#	hr;
#}

#clean_ipcs()
#{
#	for NODE in `cat ${PBS_NODEFILE} | sort -u`; do
#		ssh ${NODE} ${CLEANIPCS_BIN};
#	done;
#}

#run_mpirun()
#{
#	clean_ipcs;
#	MPI_CMD="${MPIRUN_BIN} -nolocal -np ${NPROCS} -machinefile ${PBS_NODEFILE} ${FNM_BIN}";
#    printf "begin simulation, please go to bed ...\n";
#	printf "%s\n\n" "${MPI_CMD}";
#	time ${MPI_CMD};
	#sleep 10;
#}

run_mpiexec()
{
        MPIEXEC_CMD="${MPIEXEC_BIN} ${FNM_BIN}"
        echo -e "\n${MPIEXEC_CMD}\n"
	time ${MPIEXEC_CMD};
}

main()
{
	run_mpiexec;
	hr;
}

time main;

