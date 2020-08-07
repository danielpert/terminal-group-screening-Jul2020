#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories."""
import argparse
import logging
import datetime
from copy import deepcopy
from hashlib import sha1

import signac

logging.basicConfig(filename='init.log', filemode='w', level=logging.INFO)

'''
-----------------------------
New terminal groups (11)
-----------------------------
'''
terminal_groups_new = ['cyclohexyl', 'acetylene', 'sulfhydryl']
#    ['amide', 'biphenyl', 'cyclohexyl', 'dihydroxyphenyl',
#     'formyl', 'nitroso', 'pentafluorophenyl', 'sulfhydryl',
#     'triazole', 'benzoicacid', 'isopropylbenzene']
'''
-----------------------------
Old terminal groups (19)
-----------------------------
'''
terminal_groups_old = \
    ['acetyl', 'amino', 'carboxyl', 'cyano', 'cyclopropyl',
     'difluoromethyl', 'ethylene', 'fluorophenyl', 'hydroxyl',
     'isopropyl', 'methoxy', 'methyl', 'nitro', 'nitrophenyl',
     'perfluoromethyl', 'phenol', 'phenyl', 'pyrrole', 'toluene']

    
# Initialize the project
def main(args, random_seed):
    project = signac.init_project("terminal-group-screening-Jul2020")
    logging.info("Init begin: {}".format(datetime.datetime.today()))
    logging.info("Initialized project name")
    statepoints = list()
    backbone = 'Alkylsilane'
    # generate the new top and bottom monolayers
    for replication_index in range(args.num_replicas):
        
        for tg_1 in terminal_groups_new:
            for tg_2 in terminal_groups_new:
                for tg_3 in terminal_groups_new:
                    the_statepoint = dict(
                      # Carbon backbone length
                        chainlength = args.chain_length,
                      # Backbone chemistry
                        backbone = backbone,
                      # Terminal groups
                        terminal_group_1 = tg_1,
                        terminal_group_2 = tg_2,
                        terminal_group_3 = tg_3,
                      # fractions
                        frac_1 = 0.5,
                        frac_2 = 0.5,
                      # Number of monolayer chains
                        n = 100,
                      # monolayer pattern type
                        pattern_type = args.pattern_type,
                      # Random seed
                        seed = random_seed*(replication_index + 1))
                    project.open_job(statepoint=the_statepoint).init()
                    statepoints.append(the_statepoint)
                    logging.info(msg="At the statepoint: {}".format(the_statepoint))
            
    # write statepoints to signac statepoint file
    project.write_statepoints(statepoints=statepoints)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Initialize the data space.")
    parser.add_argument(
        'random',
        type=str,
        help="A string to generate a random seed.")
    parser.add_argument(
        '-n', '--num-replicas',
        type=int,
        default=1,
        help="Initialize multiple replications.")
    parser.add_argument(
        '-c', '--chain-length',
        type=int,
        default=17,
        help="Backbone length.")
    parser.add_argument(
        '-t', '--pattern_type',
        type=str,
        default='random',
        help="Name of pattern, default \'random\'.")
    args = parser.parse_args()

    # Generate an integer from the random str.
    try:
        random_seed = int(args.random)
    except ValueError:
        random_seed = int(sha1(args.random.encode()).hexdigest(), 16) % (10 ** 8)

    logging.info("Params:\
                random: {}\
                num-replicas: {}\
                chain_length: {}\
                pattern_type: {}".format(
                    random_seed,
                    args.num_replicas,
                    args.chain_length,
                    args.pattern_type,
                    ))

    main(args, random_seed)
