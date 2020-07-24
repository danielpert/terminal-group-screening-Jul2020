{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash -l
#PBS -j oe
#PBS -l nodes={{ nn }}:ppn=16
#PBS -l walltime={{ walltime|format_timedelta }}
#PBS -q low
#PBS -m abe
#PBS -M quachcd.rahman@gmail.com

module load openmpi
module load lammps/15May15
module load gromacs/5.1.4
{% endblock %}
