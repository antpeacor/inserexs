"""
BSD 3-Clause License

Copyright (c) 2022 CNRS - Université de Strasbourg.
All rights reserved.

Author : [Antonio Pena Corredor] [antonio.penacorredor@ipcmss.unistra.fr]

This software is a reflection choice framework for Resonant Elastic X-ray Scattering.
The program consists of three different ".py" modules and a "GUI.ui" graphic interface.
- main.py: backbone, direct exchange with interface.
- intensity module.py: module for the calculation of the reflection intensities.
- sensitivity module.py: module for the calculation of the reflection sensitivities.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The program has been coded and tested on python 3.9
The version is indicated for those modules not included in the standard Python library

Module for the sensitivity calculation.
"""
import os
from pathlib import Path
import sys
import numpy as np  # v1.21.5
import pandas as pd  # v1.4.2
from itertools import combinations_with_replacement
import intensity_module as im
import promptlib  # v3.0.20


def input_generator(crystal):
    """Generate input for FDMNES."""
    chosen_atoms = [crystal.atom_list[i][0] for i in crystal.edge_checked_list]

    # List including the atomic information of the present atoms
    chosen_atomic_info = list(set([(Atomic_number(atom))
                                  for atom in chosen_atoms]))

    # chosen Z:
    chosen_atomic_numbers = [str(info[0]) for info in chosen_atomic_info]

    # Not all reflections are allowed, they should be considered:
    allowed_reflections = evaluate_ref(crystal, chosen_atomic_info)

    reflection_lines_str = '\n'.join([''.join(str(reflection[0])) + '  1  1       0.'
                                      for reflection in allowed_reflections]).replace('(', '').replace(')', '').replace(',', '')

    E_start, E_end, step = crystal.E_start, crystal.E_stop, crystal.E_step
    percentage = crystal.percent

    # To see which parameters should be refined:

    """
    The program has been designed to fetch the FDMNES program in an FDMNES
    folder located in the same home directory. This can be changed according
    to the user's preferences. However, the different working files and
    the results will be saved in this directory, so its creation is highly
    recommended before running inserexs.

    In case it does not exist, the program will ask to select where it is
    """

    global dirfdmnes  # changes variable into global type
    try:
        os.chdir(Path(home, 'FDMNES'))
        dirfdmnes = os.getcwd()
    except:
        prompter = promptlib.Files()  # calls for directory
        dirfdmnes = prompter.dir()  # changes directory
        os.chdir(Path(dirfdmnes))

    try:
        os.mkdir(f"{crystal.name}_input")
    except FileExistsError:
        pass

    Repetitions = 3

    os.chdir(Path(dirfdmnes, f"{crystal.name}_input"))

    if crystal.coupled:  # If we do not couple parameters
        coupling_matrix = [[1 + (rep - Repetitions//2)*percentage/100
                            for i in range(len(crystal.refinement_checked_list))]
                           for rep in range(Repetitions)]
    else:
        values = [1 + (rep - Repetitions//2)*percentage/100
                  for rep in range(Repetitions)]  # creates values
        coupling_matrix = list(combinations_with_replacement(values, len(crystal.refinement_checked_list)))

    locator = 0  # for the locator
    for item_set in coupling_matrix:

        list_atom_lines = []

        i = 0  # item to be modified
        for atom in range(crystal.n):

            atom_x = float(crystal.atom_list[atom][1])
            atom_y = float(crystal.atom_list[atom][2])
            atom_z = float(crystal.atom_list[atom][3])
            #atom_o = float(crystal.atom_list[atom][4])
            atom_o = 1.0

            atomic_number = Atomic_number(crystal.atom_list[atom][0])[0]

            if atom*4 in crystal.refinement_checked_list:  # x position
                atom_x = str(round(atom_x * item_set[i], 3))
                i += 1
            if atom*4 + 1 in crystal.refinement_checked_list:  # y position
                atom_y = str(round(atom_y * item_set[i], 3))
                i += 1
            if atom*4 + 2 in crystal.refinement_checked_list:  # z position
                atom_z = str(round(atom_z * item_set[i], 3))
                i += 1
            if atom*4 + 3 in crystal.refinement_checked_list:  # occupation
                atom_o = str(round(atom_o * item_set[i], 3))
                i += 1

            list_atom_lines.append(f" {atomic_number}      {atom_x}   {atom_y}   {atom_z}   {atom_o}")

        # text for atom_lines
        atom_lines = '\n'.join(list_atom_lines)

        # if the SG has a seeting, it should be considered:
        if hasattr(crystal, 'setting'):
            spacegroup = str(crystal.spacegroup) + ':' + str(crystal.setting)
        else:
            spacegroup = str(crystal.spacegroup)

        # header:
        text = ['Filout \n',
                f'FileResults\\result_{locator}',
                '\n',
                '    Range \n',
                f' {E_start}  {step}  {E_end} \n',
                '\n',
                ' Eimag',
                '  0.1 \n',
                '\n',
                'Radius \n',
                ' 5\n',
                ' \n',
                ' Green \n',
                ' Quadrupole \n',
                ' Density \n',
                ' Spherical \n',
                '\n',
                ' Z_absorber \n',
                ' '.join(chosen_atomic_numbers) + ' \n',
                '\n',
                'RXS \n',
                str(reflection_lines_str).replace('[', '').replace(']', '').replace(',', ''),
                '\n',
                'Spgroup \n',
                spacegroup,
                '\n',
                'Crystal_t \n',
                '         ' + "{:.5f}".format(crystal.a) + ' ' + "{:.5f}".format(crystal.b) + ' ' + "{:.5f}".format(crystal.c) + '  ' + str(crystal.alpha) + ' ' + str(crystal.beta) + ' ' + str(crystal.gamma) + '\n',
                str(atom_lines).replace('[', '').replace(']', '').replace(',', ''),
                '\n'
                'Convolution \n'
                '\n',
                'Estart \n',
                f'{E_start} \n',
                ' \n',
                'End \n']

        f = open(f'input_{locator}.txt', 'w')
        f.write(''.join(text))
        f.close()

        locator += 1
    return locator  # for the filenumber


def instructions(name, n):
    """Create instruction file for FDMNES."""
    text = [f'{str(n)} \n']
    for simulation in range(n):
        text += f'{name}_input\\input_{simulation}.txt \n'

    os.chdir(dirfdmnes)
    f = open('fdmfile.txt', 'w')
    f.write(''.join(text))

    # for results:
    try:
        os.mkdir('FileResults')
    except OSError:
        pass
    f.close()


def run_fdmnes():
    """Run FDMNES."""
    global dirfdmnes
    os.chdir(dirfdmnes)
    if 'win' in sys.platform:
        os.startfile('fdmnes_win64.exe')
    elif 'lin' in sys.platform:
        os.startfile('fdmnes_linux64')
    return


def Atomic_number(element):
    """Brings number from element symbol."""
    letters = [x for x in element]

    chain_no_numbers = ''.join([s for s in letters
                                if s.isdigit() == False and s.isalnum()])

    # dic for atomic numbers, Atomic weight, K edge, L1, L2, L3, M1 and M5:
    dic_atomic_numbers ={"H": (1, 1.008, 13.6, 0, 0, 0),
                         "He": (2, 4.003, 24.6, 0, 0, 0),
                         "Li": (3, 6.941, 54.8, 0, 0, 0),
                         "Be": (4, 9.012, 111, 0, 0, 0),
                         "B": (5, 10.811, 188, 0, 4.7, 4.7),
                         "C": (6, 12.011, 283.8, 0, 6.4, 6.4),
                         "N": (7, 14.007, 401.6, 0, 9.2, 9.2),
                         "O": (8, 15.999, 532, 23.7, 7.1, 7.1),
                         "F": (9, 18.998, 685.4, 31, 8.6, 8.6),
                         "Ne": (10, 20.18, 866.9, 45, 18.3, 18.3),
                         "Na": (11, 22.99, 1072.1, 63.3, 31.1, 31.1),
                         "Mg": (12, 24.305, 1305, 89.4, 51.4, 51.4),
                         "Al": (13, 26.982, 1559.6, 117.7, 73.1, 73.1),
                         "Si": (14, 28.086, 1838.9, 148.7, 99.2, 99.2),
                         "P": (15, 30.974, 2145.5, 189.3, 132.2, 132.2),
                         "S": (16, 32.066, 2472, 229.2, 164.8, 164.8),
                         "Cl": (17, 35.453, 2822.4, 270.2, 201.6, 200),
                         "Ar": (18, 39.948, 3202.9, 320, 247.3, 245.2),
                         "K": (19, 39.098, 3607.4, 377.1, 296.3, 293.6),
                         "Ca": (20, 40.078, 4038.1, 437.8, 350, 346.4),
                         "Sc": (21, 44.956, 4492.8, 500.4, 406.7, 402.2),
                         "Ti": (22, 47.88, 4966.4, 563.7, 461.5, 455.5),
                         "V": (23, 50.942, 5465.1, 628.2, 520.5, 512.9),
                         "Cr": (24, 51.996, 5989.2, 694.6, 583.7, 574.5),
                         "Mn": (25, 54.938, 6539, 769, 651.4, 640.3),
                         "Fe": (26, 55.847, 7112, 846.1, 721.1, 708.1),
                         "Co": (27, 58.933, 7708.9, 925.6, 793.8, 778.6),
                         "Ni": (28, 58.69, 8332.8, 1008.1, 871.9, 854.7),
                         "Cu": (29, 63.546, 8978.9, 1096.1, 951, 931.1),
                         "Zn": (30, 65.39, 9658.6, 1193.6, 1042.8, 1019.7),
                         "Ga": (31, 69.723, 10367.1, 1297.7, 1142.3, 1115.4),
                         "Ge": (32, 72.61, 11103.1, 1414.3, 1247.8, 1216.7),
                         "As": (33, 74.922, 11866.7, 1526.5, 1358.6, 1323.1),
                         "Se": (34, 78.96, 12657.8, 1653.9, 1476.2, 1435.8),
                         "Br": (35, 79.904, 13473.7, 1782, 1596, 1549.9),
                         "Kr": (36, 83.8, 14325.6, 1921, 1727.2, 1674.9),
                         "Rb": (37, 85.468, 15199.7, 2065.1, 1863.9, 1804.4),
                         "Sr": (38, 87.62, 16104.6, 2216.3, 2006.8, 1939.6),
                         "Y": (39, 88.906, 17038.4, 2372.5, 2155.5, 2080),
                         "Zr": (40, 91.224, 17997.6, 2531.6, 2306.7, 2222.3),
                         "Nb": (41, 92.906, 18985.6, 2697.7, 2464.7, 2370.5),
                         "Mo": (42, 95.94, 19999.5, 2865.6, 2625.1, 2520.2),
                         "Tc": (43, 98.906, 21044, 3042.5, 2793.2, 2676.9),
                         "Ru": (44, 101.07, 22117.2, 3224, 2966.9, 2837.9),
                         "Rh": (45, 102.906, 23219.9, 3411.9, 3146.1, 3003.8),
                         "Pd": (46, 106.42, 24350.3, 3604.3, 3330.3, 3173.3),
                         "Ag": (47, 107.868, 25514, 3805.8, 3523.7, 3351.1),
                         "Cd": (48, 112.411, 26711.2, 4018, 3727, 3537.5),
                         "In": (49, 114.82, 27939.9, 4237.5, 3938, 3730.1),
                         "Sn": (50, 118.71, 29200.1, 4464.7, 4156.1, 3928.8),
                         "Sb": (51, 121.75, 30491.2, 4698.3, 4380.4, 4132.2),
                         "Te": (52, 127.6, 31813.8, 4939.2, 4612, 4341.4),
                         "I": (53, 126.904, 33169.4, 5188.1, 4852.1, 4557.1),
                         "Xe": (54, 131.29, 34561.4, 5452.8, 5103.7, 4782.2),
                         "Cs": (55, 132.905, 35984.6, 5714.3, 5359.4, 5011.9),
                         "Ba": (56, 137.327, 37440.6, 5988.8, 5623.6, 5247),
                         "La": (57, 138.906, 38924.6, 6266.3, 5890.6, 5482.7),
                         "Ce": (58, 140.115, 40443, 6548.8, 6164.2, 5723.4),
                         "Pr": (59, 140.908, 41990.6, 6834.8, 6440.4, 5964.3),
                         "Nd": (60, 144.24, 43568.9, 7126, 6721.5, 6207.9),
                         "Pm": (61, 146.915, 45184, 7427.9, 7012.8, 6459.3),
                         "Sm": (62, 150.36, 46834.2, 7736.8, 7311.8, 6716.2),
                         "Eu": (63, 151.965, 48519, 8052, 7617.1, 6976.9),
                         "Gd": (64, 157.25, 50239.1, 8375.6, 7930.3, 7242.8),
                         "Tb": (65, 158.925, 51995.7, 8708, 8251.6, 7514),
                         "Dy": (66, 162.5, 53788.5, 9045.8, 8580.6, 7790.1),
                         "Ho": (67, 164.93, 55617.7, 9394.2, 8917.8, 8071.1),
                         "Er": (68, 167.26, 57485.5, 9751.3, 9264.3, 8357.9),
                         "Tm": (69, 168.934, 59389.6, 10115.7, 9616.9, 8648),
                         "Yb": (70, 173.04, 61332.3, 10486.4, 9978.2, 8943.6),
                         "Lu": (71, 174.967, 63313.8, 10870.4, 10348.6, 9244.1),
                         "Hf": (72, 178.49, 65350.8, 11270.7, 10739.4, 9560.7),
                         "Ta": (73, 180.948, 67416.4, 11681.5, 11136.1, 9881.1),
                         "W": (74, 183.85, 69525, 12099.8, 11544, 10206.8),
                         "Re": (75, 186.207, 71676.4, 12526.7, 11958.7, 10535.3),
                         "Os": (76, 190.2, 73870.8, 12968, 12385, 10870.9),
                         "Ir": (77, 192.22, 76111, 13418.5, 12824.1, 11215.2),
                         "Pt": (78, 195.08, 78394.8, 13879.9, 13272.6, 11563.7),
                         "Au": (79, 196.967, 80724.9, 14352.8, 13733.6, 11918.7),
                         "Hg": (80, 200.59, 83102.3, 14839.3, 14208.7, 12283.9),
                         "Tl": (81, 204.383, 85530.4, 15346.7, 14697.9, 12657.5),
                         "Pb": (82, 207.2, 88004.5, 15860.8, 15200, 13035.2),
                         "Bi": (83, 208.98, 90525.9, 16387.6, 15711.1, 13418.6),
                         "Po": (84, 208.982, 93105, 16939.3, 16244.3, 13813.8),
                         "At": (85, 209.987, 95729.9, 17493, 16784.7, 14213.5),
                         "Rn": (86, 222.018, 98404, 18049, 17337.1, 14619.4),
                         "Fr": (87, 223.02, 101137, 18639, 17906.5, 15031.2),
                         "Ra": (88, 226.025, 103921.9, 19236.7, 18484.3, 15444.4),
                         "Ac": (89, 227.028, 106755.3, 19840, 19083.2, 15871),
                         "Th": (90, 232.038, 109650.9, 20472.1, 19693.2, 16300.3),
                         "Pa": (91, 231.036, 112601.4, 21104.6, 20313.7, 16733.1),
                         "U": (92, 238.029, 115606.1, 21757.4, 20947.6, 17166.3),
                         "Np": (93, 237.048, 118678, 22426.8, 21600.5, 17610),
                         "Pu": (94, 244.064, 121818, 23097.2, 22266.2, 18056.8),
                         "Am": (95, 243.061, 125027, 23772.9, 22944, 18504.1),
                         "Cm": (96, 247.07, 128200, 24460, 23779, 18930),
                         "Bk": (97, 247.07, 131590, 25275, 24385, 19452),
                         "Cf": (98, 251.08, 135960, 26110, 25250, 19930),
                         "Es": (99, 252.083, 139490, 26900, 26020, 20410),
                         "Fm": (100, 257.095, 143090, 27700, 26810, 20900),
                         "Md": (101, 258.099, 146780, 28530, 27610, 21390),
                         "No": (102, 259.101, 150540, 29380, 28440, 21880),
                         "Lr": (103, 260.105, 154380, 30240, 29280, 22360)}

    return dic_atomic_numbers[chain_no_numbers]


def evaluate_ref(crystal, atomic_info_list):
    """Call to evaluate if reflections are geometrically feasible."""
    allowed_reflection_list = crystal.reflections
    for info in atomic_info_list:

        # retrieve K edge from previous info:
        wlength = 1.23984198e4/float(info[2])  # lambda in angstroms

        for reflection in crystal.reflections:
            try:
                im.angle_get(reflection[0], crystal.lattice_dimensions, wlength)
            except Exception as ex:  # if it is not allowed, removal
                allowed_reflection_list.remove(reflection)

    return allowed_reflection_list


def sensitivity_calculation(crystal):
    """Calculate sensitivity."""
    reflections_and_I = crystal.reflections

    reflection_list = [str(row[0]).replace(',', '').replace('(', '').replace(')', '')
                       for row in reflections_and_I]

    length = int((crystal.E_stop - crystal.E_start)/crystal.E_step + 1)

    intensity_matrix = np.zeros((len(reflections_and_I), length, crystal.filenumber)) # we make a matrix where we will keep the intensity values for both calculations. final range is for average

    global dirfdmnes
    os.chdir(Path(dirfdmnes, 'FileResults'))
    for repetition in range(crystal.filenumber):
        data = np.loadtxt(f'result_{repetition}_conv.txt', skiprows=1).T[2:]
        intensity_matrix[:, :, repetition] = data
    intensities = [np.mean(reflection_matrix).round(0)
                   for reflection_matrix in intensity_matrix]

    # we define two types of sensitivities
    sensitivity_matrix_N = [[[((repetition_value - np.mean(energy_value))/np.mean(energy_value))**2
                            for repetition_value in energy_value]
                            for energy_value in reflection_matrix]
                            for reflection_matrix in intensity_matrix]

    sensitivity_matrix_I = [[[(repetition_value - np.mean(energy_value))**2
                            for repetition_value in energy_value]
                            for energy_value in reflection_matrix]
                            for reflection_matrix in intensity_matrix]

    sensitivities_N = [np.mean(sensitivity_matrix)
                       for sensitivity_matrix in sensitivity_matrix_N]

    sensitivities_I = [np.mean(sensitivity_matrix).round(0)
                       for sensitivity_matrix in sensitivity_matrix_I]

    intensities = intensities/max(intensities)*100
    sensitivities_I = sensitivities_I/max(sensitivities_I)*100
    sensitivities_N = sensitivities_N/max(sensitivities_N)*100

    weights_I = [sensitivities_I[i] * intensities[i]
                 for i in range(len(sensitivities_I))]
    weights_N = [sensitivities_N[i] * intensities[i]
                 for i in range(len(sensitivities_N))]

    weights_I = weights_I/max(weights_I)*100
    weights_N = weights_N/max(weights_N)*100

    results = pd.DataFrame(list(zip(reflection_list, intensities, sensitivities_I, sensitivities_N, weights_I, weights_N)),
                           columns=['Reflections', 'Intensities', 'Sensitivities', 'Sensitivities (I normalised)', 'Weights', 'Weights (I normalised)'])

    results.to_csv('results.csv', sep = '\t')
    return results


# dictionary, we replace ' to "

home = os.getcwd()
