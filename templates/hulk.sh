{% extends "base_script.sh" %}
{% block header %}
#!/bin/bash -l
#PBS -j oe
#PBS -l nodes={{ nn }}:ppn=16
#PBS -l walltime={{ walltime|format_timedelta }}
#PBS -q batch
#PBS -m abe
#PBS -M dpert@umich.edu
#PBS -e /home/danielpert/terminal-group-screening-Jul2020/templates/hulk_error.txt
#PBS -o /home/danielpert/terminal-group-screening-Jul2020/templates/hulk_output.txt

module load openmpi/3.1.0
module load lammps/30Jun2020
module load gromacs/2018.7
module load anaconda/5.1.0
source activate imodels37

{% endblock %}
