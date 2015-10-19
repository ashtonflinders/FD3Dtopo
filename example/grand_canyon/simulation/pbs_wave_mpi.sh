#!/bin/bash
# modified from PSSClabs scripts

### Number of nodes n nodes using m cpus Processor Per Node
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00

### Queue name
#PBS -q default
### Job name
#PBS -N g.canyon.wave

### Merge stderr with stdout
#PBS -j oe
### Mail to user
#PBS -m eb
### Declare job-non-rerunable
#PBS -r n

### Generate real-time output in user home directory
## PBS -k oe

PBS_PWD="`pwd`";
cd "${PBS_O_WORKDIR}";
#THIS_HOST="`hostname`";
THIS_HOST=$(hostname --short)
DATE=$(date)

#MPICH_ROOT="/opt/mpich/p4-intel";
#MPICH_BUILD="p4.shm-intel10";
#MPICH_ROOT="/opt/mpich1/${MPICH_BUILD}";
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

# Intel library path
#export LD_LIBRARY_PATH="/opt/intel/cce/10/lib:/opt/intel/fce/10/lib"

#MPIRUN_BIN="${MPICH_ROOT}/bin/mpirun";
#MPIEXEC_BIN="/opt/mpiexec/bin/mpiexec";
MPIEXEC_BIN="${MPI_ROOT}/bin/mpiexec -n ${NPROCS} -f ${PBS_NODEFILE}"

FNM_BIN="./bin/seis3d_wave_mpi";

hr()
{
    perl -e 'print "\n" . "-"x70 . "\n\n"';
}

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

