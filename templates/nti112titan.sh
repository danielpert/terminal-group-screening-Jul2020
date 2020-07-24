{# Templated in accordance with: https://www.olcf.ornl.gov/for-users/system-user-guides/titan/running-jobs/ #}
{% set mpiexec = "aprun" %}
{% extends "torque.sh" %}
{% block tasks %}
{% set threshold = 0 if force else 0.9 %}
{% set cpu_tasks = operations|calc_tasks('np', parallel, force) %}
{% set gpu_tasks = operations|calc_tasks('ngpu', parallel, force) %}
{% if wraprun %}
{% set nn_cpu = operations|length %}
{% else %}
{% set nn_cpu = cpu_tasks|calc_num_nodes(16) %}
{% endif %}
{% set nn_gpu = gpu_tasks|calc_num_nodes(1) %}
{% set nn = nn|default((nn_cpu, nn_gpu)|max, true) %}
#PBS -l nodes={{ nn|check_utilization(gpu_tasks, 1, threshold, 'GPU') }}
{% endblock %}
{% block header %}
#!/bin/bash -l
{{ super() -}}
{% set account = account|default(environment|get_account_name, true) %}
{% if account %}
#PBS -A {{ account }}
{% endif %}
#PBS -j oe
#PBS -l gres=atlas1%atlas2
#PBS -m abe
#PBS -M gilmerjb.job.scheduler@gmail.com
#PBS -A NTI112
#PBS -V
{% if debug_q %}
#PBS -q debug
{% endif %}
{% endblock %}
{% block project_header %}
set -e
set -u
cd {{ project.config.project_dir }}
{% if lammps %}
module swap PrgEnv-pgi PrgEnv-gnu
module load fftw
module load lammps
{% endif %}
{% if gromacs %}
export CRAY_CUDA_MPS=1
module load gromacs/5.1.0
{% endif %}
{% if wraprun %}
module load python wraprun
{% endif %}
{% if init %}
module load python_anaconda3/3.6.4-anaconda3-5.1.0
export LD_LIBRARY_PATH="${CONDALIBS}:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="/ccs/proj/nti112/incite35/lib:$LD_LIBRARY_PATH"
module load glib/2.36.4
export PATH="/ccs/proj/nti112/incite35/bin:$PATH"
{% endif %}
{% endblock %}
{% block body %}
{% set cmd_suffix = cmd_suffix|default('') ~ (' &' if parallel else '') %}
{% if wraprun %}
wraprun \
{% endif %}
{% for operation in operations %}
{% if operation.directives.nranks and not mpi_prefix %}
{% set mpi_prefix = "%s -n %d "|format(mpiexec|default("mpiexec"), operation.directives.nranks) %}
{% endif %}
{% if operation.directives.omp_num_threads %}
export OMP_NUM_THREADS={{ operation.directives.omp_num_threads }}
{% endif %}
{% if wraprun and not concat_trj %}
{{ ' -n '}}{{ operation.directives.nranks}} {{ '--w-cd ' }}{{ operation.cmd }} {{ ' : \\' if not loop.last else '' }}
{% elif wraprun and concat_trj %}
{{ operation.cmd }} {{ ' : \\' if not loop.last else '' }}
{% else %}
{{ cmd_prefix }}{{ operation.cmd }}{{ cmd_suffix }}
{% endif %}
{% endfor %}
{% endblock %}
{% block footer %}
{% if parallel and not wraprun%}
wait
{% endif %}
{% endblock %}
