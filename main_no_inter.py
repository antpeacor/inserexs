"""
BSD 3-Clause License.

Copyright (c) 2022 CNRS - Université de Strasbourg.
All rights reserved.

Author : [Antonio Pena Corredor] [antonio.penacorredor@ipcms.unistra.fr]

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

Main module, backbone.
"""

import os
import sys
import numpy as np  # v1.21.5
import math
import threading  # v4.1.0
import matplotlib.pyplot as plt  # v3.5.1
import seaborn as sns #v0.12.1
from glob import glob # v11.7
import CifFile #v4.4.5

# Own modules, to be placed in same folder:
import intensity_module as im
import sensitivity_module as sm


home_cwd = os.getcwd()  # current cwd

class backbone():
    """Class that performs all the GUI operations."""

    def __init__(self, parent=None, *args, **kwargs):

        self.racine_directory = home_cwd
        self.load_cif()
        self.fetch_reflections()
        self.launchfunction()
        self.sencalcul()


    def load_cif(self):
        """Call when button "Load .cif" is called."""
        project_files = glob('*.cif')
        if (len(project_files) == 0):
            raise TypeError("No .cif files")
        elif len(project_files) == 1:
            file_cif = project_files[0]
        else:
            print('Choose .cif file: \n')
            for i in range(len(project_files)):
                print(str(i) + ' - ' + project_files[i])

            selection = input('Your selection (just the number): ')
            sel_filtered = ''.join([c for c in selection if c.isalnum()])
            file_cif = project_files[int(sel_filtered)]

        # load crystal info
        cf_object = CifFile.ReadCif(file_cif)
        get_crystal_info(crystal, cf_object)

    def fetch_reflections(self):
        """Call to fetch the reflections."""
        # to consider forbidden reflections or not

        forbidden = input('Consider forbidden reflections? [y/n]: ')
        crystal.forbidden = forbidden.lower == 'y'

        # for the max value of h, k, l:
        crystal.maxhkl = int(input('Maximum h, k or l: '))

        thread_FR = myThread(fun_fetch_reflections, crystal)

        thread_FR.start()
        thread_join(thread_FR)  # we join the thread to the rest
        #
        crystal.reflections = thread_FR.output  # we take the output away
        thread_FR.stop()
        self.refinement_checks()


    def refinement_checks(self):
        """Call to check the buttons which are checked."""
        print(' n | element |     x     |     y     |     z     | occupation')

        for i in range(len(crystal.atom_list)):
            line = crystal.atom_list[i]
            print(' ' + str(i) + ' '*(2 - len(str(i))), end='|')
            print(' '*4 + line[0] + ' '*(5 - len(line[0])), end='|')
            print(' ' + line[1] + ' '*(10 - len(line[1])), end='|')
            print(' ' + line[2] + ' '*(10 - len(line[2])), end='|')
            print(' ' + line[3] + ' '*(10 - len(line[3])), end='|')
            print(' ' + line[4], end='\n')

        print('Enter the parameters that you want to refine:')
        print('Example: 0 = x y, 1 = z occ (x and y for atom 0, occupation and z for atom 2)')

        refinement_choice = input('Your choice: ')
        refinement_items = refinement_choice.strip().lower().split(',')

        crystal.refinement_checked_list = []

        for item in refinement_items:
            item_splitted = item.split('=')
            if 'x' in item_splitted[1]:
                crystal.refinement_checked_list.append(int(item_splitted[0])*4 + 0)
            if 'y' in item_splitted[1]:
                crystal.refinement_checked_list.append(int(item_splitted[0])*4 + 1)
            if 'z' in item_splitted[1]:
                crystal.refinement_checked_list.append(int(item_splitted[0])*4 + 2)
            if 'occ' in item_splitted[1]:
                crystal.refinement_checked_list.append(int(item_splitted[0])*4 + 3)

        if len(crystal.refinement_checked_list) == 0:
            raise ValueError("No items have been chosen.")

        edge_choice = input('Enter the number of the element whose edge you wish to explore: ')
        crystal.edge_checked_list = [int(edge_choice)]

        if len(crystal.edge_checked_list) == 0:
            raise ValueError("No edges have been chosen.")

    def fetch_sensitivity(self):
        """Call to fetch the sensitivity."""
        thread_FS = myThread(fun_fetch_sensitivity, crystal)
        thread_FS.start()
        thread_join(thread_FS)  # we join the thread to the rest

        crystal.filenumber = int(thread_FS.output)  # we take the output away
        thread_FS.stop()

    def launch_fdmnes(self):
        """Call to launch FDMNES."""
        fun_run_fdmnes()

    def launchfunction(self):
        """Call to launch the function."""
        print('Energies for the simulation (eV) - 0 is the edge')
        crystal.E_start = float(input('Energy at which the simulation will start: '))
        crystal.E_stop = float(input('Energy at which the simulation will end: '))
        crystal.E_step = float(input('Step between each simulation: '))
        crystal.percent = float(input('% variation for the chosen parameter: '))

        # to see if parameters should be coupled:
        crystal.coupled = input('Do you want to couple the parameters? [y/n]').lower().strip() == 'y'
        self.fetch_sensitivity()  # launches for sensitivity calculation
        self.launch_fdmnes()  # launches FDMNES

    def sencalcul(self):
        """Call to calculate the intensity."""
        thread_SC = myThread(fun_sen_calcul, crystal)
        thread_SC.start()
        thread_join(thread_SC)  # we join the thread to the rest

        thread_SC.stop()

        self.representation()

    def representation(self):
        """Call to do final plot."""
        MPLplot()


class Crystal:
    """Class for the creation of a crystal, object."""

    def __init__(self):
        self.atom_list = []
        self.n = 0  # number of present atoms

        self.operation_list = []
        self.nsym = 0

    def add(self, element):
        """Call to add an extra element."""
        self.atom_list.append(element)

    def add_sym(self, element):
        """Call to add an symmetry element."""
        self.operation_list.append(element)


class myThread(threading.Thread):
    """Create a thread object."""

    output = ''

    def __init__(self, function, arg):
        threading.Thread.__init__(self)
        self.function = function
        self.arg = arg

        self._stop_event = threading.Event()

    def run(self):
        """Run the thread's function."""
        global threadLock
        threadLock.acquire()  # get lock to synchro threads
        self.output = self.function(self.arg)  # the function is executed
        threadLock.release()  # release

    def stop(self):
        """Stop the active thread."""
        self._stop_event.set()


def MPLplot():
    """Call this function for the representation."""
    fig, ax = plt.subplots(1, figsize=(5, 5), facecolor='white')
    plt.rcParams['font.size'] = 15
    ax.xaxis.get_label().set_fontsize(15)
    ax.yaxis.get_label().set_fontsize(15)

    x = crystal.results['Intensities']
    y = crystal.results['Sensitivities']

    # list containing the occupied coordinates. If a new reflection is too
    occ_coordinates = []

    colorgen = iter(plt.cm.rainbow(np.linspace(0, 1, len(crystal.results['Reflections'])))) # to put each label in a color

    for i in range(len(crystal.results['Reflections'])):
        x_text = crystal.results['Intensities'][i]
        y_text = crystal.results['Sensitivities'][i]

        x_text = coordinate_comeback(x_text)
        y_text = coordinate_comeback(y_text)

        color_j = next(colorgen)

        # close to an existing one, the label is not shown
        for coord in occ_coordinates:
            if math.sqrt((x_text - coord[0])**2 + (y_text - coord[1])**2) < 10:
                break
        else:
            ax.text(x_text, y_text, str(crystal.results['Reflections'][i]),
                         color=color_j, fontsize=15)
            ax.scatter(x[i], y[i], color=color_j)

        occ_coordinates.append((x_text, y_text))

    plt.xlim([-5, 105]), plt.ylim([-5, 105])
    sns.set_style('ticks')

    ax.set_xlabel('Intensity (%)', fontsize=18)
    ax.set_ylabel('Sensitivity (%)', fontsize=18)

    plt.subplots_adjust(left=0.2, bottom=0.2)
    fig.text(0,0,'.', color='white') # positioning issue

    plt.savefig(crystal.name + '.png',  dpi = 500, bbox_inches = 'tight')


def coordinate_comeback(coord):
    """Call to bring reflections into the representation."""
    if coord + 2 <= 100:
        return coord + 2
    else:
        return coord - 9


def get_crystal_info(obj, CIF):
    """
    We give crystal object, sym op object and CIF object.

    It is conceived for cif files with a single data block
    """
    # Only loads the first data block
    CIFdata = CIF[CIF.keys()[0]]

    # Loads each parameter
    obj.a = float(CIFdata['_cell_length_a'])
    obj.b = float(CIFdata['_cell_length_b'])
    obj.c = float(CIFdata['_cell_length_c'])

    obj.alpha = float(CIFdata['_cell_angle_alpha'])
    obj.beta = float(CIFdata['_cell_angle_beta'])
    obj.gamma = float(CIFdata['_cell_angle_gamma'])

    sg_info = CIFdata['_symmetry_space_group_name_H-M'].replace(' ', '')
    if ':' in sg_info:
        sg_split = sg_info.split(':')
        obj.spacegroup, obj.setting = sg_split[0], sg_split[1]
    else:
        obj.spacegroup = sg_info

    if '_space_group_symop_operation_xyz' in CIFdata.keys():
        operation_list = CIFdata['_space_group_symop_operation_xyz']
    elif '_symmetry_equiv_pos_as_xyz' in CIFdata.keys():
        operation_list = CIFdata['_symmetry_equiv_pos_as_xyz']
    else:
        raise Exception('No symmetry operations in .cif file')

    obj.operation_list = [item.split(',') for item in operation_list]

    # Atomic positions:
    atom_list = CIFdata['_atom_site_label']
    x_list = CIFdata['_atom_site_fract_x']
    y_list = CIFdata['_atom_site_fract_y']
    z_list = CIFdata['_atom_site_fract_z']

    if '_atom_site_occupancy' not in CIFdata.keys():
        occ_list = np.ones(len(x_list))
    else:
        occ_list = CIFdata['_atom_site_occupancy']

    for index in range(len(x_list)):
        obj.add((atom_list[index], x_list[index], y_list[index],
                 z_list[index], str(occ_list[index])))


    obj.n = len(atom_list)
    obj.nsym = len(obj.operation_list)
    return


def fun_fetch_reflections(obj):
    """Fecth the reflections."""
    hkl_and_Int = im.intensity_calculation(obj)

    global crystal
    crystal.reflections = hkl_and_Int
    return hkl_and_Int


def fun_fetch_sensitivity(obj):
    """Fetch the sensitivity."""
    filenumber = sm.input_generator(obj)
    return filenumber


def fun_run_fdmnes(a=0):
    """Run the FDMNES program."""
    sm.instructions(crystal.name, crystal.filenumber)
    sm.run_fdmnes()


def fun_sen_calcul(crystal):
    """Call to calculate the sensitivity."""
    input('Press enter when FDMNES is done: ')
    crystal.results = sm.sensitivity_calculation(crystal)
    return


def thread_join(thread):
    """Call this function to join the thread."""
    global threads
    threads.append(thread)
    for thr in threads:
        thr.join()
    return


if __name__ == '__main__':  # only executes the below code if  main

    threads = []  # list to save all threads
    threadLock = threading.Lock()  # for the threads
    crystal = Crystal()

    start = backbone()
    try:
        sys.exit()  # execution
    except BaseException:
        print('The plot is in the FileResults file.')
        print('Thanks for using inserexs.')
