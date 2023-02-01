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
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FC  # v3.5.1
from PyQt5 import QtCore, QtGui, uic, QtWidgets  # v5.14.4
from PyQt5.QtWidgets import QApplication, QTableWidgetItem  # v5.14.4
from pathlib import Path
import promptlib  # v3.0.20
import threading  # v4.1.0
import matplotlib.pyplot as plt  # v.3.5.1
import seaborn as sns # v0.12.1
import CifFile # v4.4.5
import pyxtal # v0.5.5

# Own modules, to be placed in same folder:
import intensity_module as im
import sensitivity_module as sm


home_cwd = os.getcwd()  # current cwd
Ui_path = Path(home_cwd, 'GUI.ui')
Ui_myWindow_base, Ui_myWindow_form = uic.loadUiType(Ui_path)


class Screen (Ui_myWindow_base, Ui_myWindow_form):
    """Main GUI."""

    def __init__(self, parent=None, *args, **kwargs):
        super(Screen, self).__init__(parent=parent, *args, **kwargs)
        self.setupUi(self)
        # self.setWindowIcon(QtGui.QIcon('Icon.ico'))  # for the windows icon
        self.racine_directory = home_cwd

        self.setWindowIcon(QtGui.QIcon('icon.ico'))

        # Widgets:
        self.loadcifButton.clicked.connect(self.load_cif)
        self.updateButton.clicked.connect(self.update_values)
        self.fetchreflectionsButton.clicked.connect(self.fetch_reflections)
        self.launchButton.clicked.connect(self.launchfunction)
        self.sensitivityButton.clicked.connect(self.sencalcul)

    def load_cif(self):
        """Call when button "Load .cif" is called."""
        prompter = promptlib.Files()  # calls for directory
        file_cif = prompter.file()  # calles for file

        crystal.name = os.path.basename(file_cif[:-4])  # for file's name

        if not file_cif.endswith('.cif'):
            raise TypeError("That is not a .cif file.")

        # write file name in correct line
        self.cifnameLine.setText(os.path.basename(file_cif))

        # load crystal info
        cf_object = CifFile.ReadCif(file_cif)
        get_crystal_info(crystal, cf_object)

        # Crystal variables:
        self.line_comb = ((self.aLine, crystal.a),
                          (self.bLine, crystal.b),
                          (self.cLine, crystal.c),
                          (self.alphaLine, crystal.alpha),
                          (self.betaLine, crystal.beta),
                          (self.gammaLine, crystal.gamma))

        # we set parameters according to values
        for (line, element) in self.line_comb:
            line.setText(str(element))

        # Table for atomic positions
        self.positionTable.setRowCount(crystal.n)
        self.edgeTable.setRowCount(crystal.n)

        for row in range(len(crystal.atom_list)):
            column_values = 0
            for column in range(len(crystal.atom_list[row])):
                item = QTableWidgetItem((crystal.atom_list[row][column]))
                item.setTextAlignment(QtCore.Qt.AlignHCenter)
                self.positionTable.setItem(row, column, item)
                column_values += 1

                if column == 0:
                    item = QTableWidgetItem((crystal.atom_list[row][column]))
                    item.setTextAlignment(QtCore.Qt.AlignHCenter)
                    self.edgeTable.setItem(row, column, item)

        # the refinement buttons are created:
        self.make_refinement_buttons()

    def update_values(self):
        """Call to update values."""
        crystal.a = self.aLine.text()
        crystal.b = self.bLine.text()
        crystal.c = self.cLine.text()
        crystal.alpha = self.alphaLine.text()
        crystal.beta = self.betaLine.text()
        crystal.gamma = self.gammaLine.text()

    def fetch_reflections(self):
        """Call to fetch the reflections."""
        # Checks if anything is checked in the edge list
        self.refinement_checks()

        # to consider forbidden reflections or not
        crystal.forbidden = self.forbiddenCheckbox.isChecked()

        # for the max value of h, k, l:
        crystal.maxhkl = int(self.maxhklLine.text())

        thread_FR = myThread(fun_fetch_reflections, crystal)

        thread_FR.start()
        thread_join(thread_FR)  # we join the thread to the rest
        #
        reflections = thread_FR.output  # we take the output away
        thread_FR.stop()

        row = 0
        self.reflectionTable.setRowCount(len(reflections))
        for reflection in reflections:
            self.reflectionTable.setItem(row, 0, QtWidgets.QTableWidgetItem
                                         (str(reflection[0]).replace(',', '')))
            self.reflectionTable.setItem(row, 1, QtWidgets.QTableWidgetItem
                                         (str(reflection[1])))
            row += 1


    def make_refinement_buttons(self):
        """Call to make refinement buttons."""
        self.groupButton = QtWidgets.QButtonGroup(self.positionTable)
        self.edgeButtons = QtWidgets.QButtonGroup(self.edgeTable)

        self.groupButton.setExclusive(False)
        self.edgeButtons.setExclusive(True)

        # to store and check ids:
        self.id_list_refinement = [[0 for column in range(5)]
                                   for column in range(len(crystal.atom_list))]
        self.id_list_edges = [0 for row in range(len(crystal.atom_list))]

        for row in range(len(crystal.atom_list)):
            for column in range(1, 5):
                ch_box = QtWidgets.QCheckBox()
                self.groupButton.addButton(ch_box)
                self.id_list_refinement[row][column] = id(ch_box)

                self.positionTable.setCellWidget(row, column, ch_box)

            # in the edge's table:
            ch_box_1 = QtWidgets.QCheckBox()
            self.edgeButtons.addButton(ch_box_1)
            self.id_list_edges[row] = id(ch_box_1)

            self.edgeTable.setCellWidget(row, 0, ch_box_1)

        # table for the storage of refinement values:
        self.refinement_checklist = np.zeros((len(crystal.atom_list), 4))

    def refinement_checks(self):
        """Call to check the buttons which are checked."""
        crystal.edge_checked_list = [i for i, button in
                                     enumerate(self.edgeButtons.buttons())
                                     if button.isChecked()]

        if hasattr(crystal, 'reflections'):
            crystal.refinement_checked_list = [i for i, button in
                                               enumerate(self.groupButton.buttons())
                                               if button.isChecked()]
            if len(crystal.refinement_checked_list) == 0:
                raise ValueError("No items have been chosen.")

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
        crystal.E_start = float(self.e_startLine.text())
        crystal.E_stop = float(self.e_endLine.text())
        crystal.E_step = float(self.e_stepLine.text())
        crystal.percent = float(self.percentLine.text())

        # to see if parameters should be coupled:
        crystal.coupled = self.coupleCheckbox.isChecked()

        # number of simulations to perform:
        crystal.nsim = float(self.nsimLine.text())

        self.refinement_checks()  # checks the checked boxes
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
        # Check if absolute sensitivity or normalised
        crystal.normalised = self.normalisedCheckbox.isChecked()

        self.graphic = MPLplot()

        self.graphicLayout.addWidget(self.graphic)


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

    def check_hex(self):
        """Check if the Miller indices should be written as hkil."""
        if self.sgnumber < 143 or self.sgnumber > 194:
            return False
        if self.sgnumber < 168 and hasattr(self, 'setting'):
            if 'h' in self.setting.lower():
                return True
            else:
                return False
        if self.sgnumber >= 168:
            return True

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


class MPLplot(FC):
    """Class for the representation."""

    def __init__(self, parent=None):
        self.fig, self.ax = plt.subplots(1, figsize=(5, 5), facecolor='white')
        plt.rcParams['font.size'] = 15
        self.ax.xaxis.get_label().set_fontsize(15)
        self.ax.xaxis.get_label().set_fontsize(15)

        super().__init__(self.fig)

        x = crystal.results['Intensities']
        if crystal.normalised:
            y = crystal.results['Sensitivities (I norm)']
        else:
            y = crystal.results['Sensitivities']

        # list containing the occupied coordinates. If a new reflection is too
        occ_coordinates = []

        colorgen = iter(plt.cm.rainbow(np.linspace(0, 1, len(crystal.results['Reflections'])))) # to put each label in a color

        for i in range(len(crystal.results['Reflections'])):
            x_text = crystal.results['Intensities'][i]

            if crystal.normalised:
                y_text = crystal.results['Sensitivities (I norm)'][i]
            else:
                y_text = crystal.results['Sensitivities'][i]

            x_text = self.coordinate_comeback(x_text)
            y_text = self.coordinate_comeback(y_text)

            color_j = next(colorgen)

            # close to an existing one, the label is not shown
            for coord in occ_coordinates:
                if math.sqrt((x_text - coord[0])**2 + (y_text - coord[1])**2) < 20:
                    break
            else:
                self.ax.text(x_text, y_text, str(crystal.results['Reflections'][i]),
                             color=color_j, fontsize=15)
                self.ax.scatter(x[i], y[i], color=color_j)

            occ_coordinates.append((x_text, y_text))

        plt.xlim([-5, 105]), plt.ylim([-5, 105])
        sns.set_style('ticks')

        self.ax.set_xlabel('Intensity (%)', fontsize=18)

        if crystal.normalised:
            self.ax.set_ylabel('Sensitivity (norm.) (%)', fontsize=18)
        else:
            self.ax.set_ylabel('Sensitivity (%)', fontsize=18)

        plt.subplots_adjust(left=0.2, bottom=0.2)
        self.fig.text(0,0,'.', color='white') # positioning issue

        plt.savefig(crystal.name + '.png',  dpi = 500, bbox_inches = 'tight')

    def coordinate_comeback(self, coord):
        """Call to bring reflections into the representation."""
        if coord + 2 <= 100:
            return coord + 2
        else:
            return coord - 9

def error_strip(a):
    """Check if a parameter has an error in it and strip it."""
    if '(' in str(a):
        return a[:-3]
    else:
        return a


def get_crystal_info(obj, CIF):
    """
    We give crystal object, sym op object and CIF object.

    It is conceived for cif files with a single data block
    """
    # Only loads the first data block
    CIFdata = CIF[CIF.keys()[0]]

    # Loads each parameter
    obj.a = float(error_strip(CIFdata['_cell_length_a']))
    obj.b = float(error_strip(CIFdata['_cell_length_b']))
    obj.c = float(error_strip(CIFdata['_cell_length_c']))

    obj.alpha = float(error_strip(CIFdata['_cell_angle_alpha']))
    obj.beta = float(error_strip(CIFdata['_cell_angle_beta']))
    obj.gamma = float(error_strip(CIFdata['_cell_angle_gamma']))

    sg_info = CIFdata['_symmetry_space_group_name_H-M'].replace(' ', '')
    if ':' in sg_info:
        sg_split = sg_info.split(':')
        obj.spacegroup, obj.setting = sg_split[0], sg_split[1]
    else:
        obj.spacegroup = sg_info

    # For the number of the space group
    if hasattr(obj, '_space_group_IT_number'):
        obj.sgnumber = int(CIF['_space_group_IT_number'])
    else: # retrieve number from space group symbol
        obj.sgnumber = int((pyxtal.symmetry.Group(obj.spacegroup)).number)

    # Check if hexagonal space group - for Miller-Bravais notation
    obj.hex = obj.check_hex()

    if '_space_group_symop_operation_xyz' in CIFdata.keys():
        operation_list = CIFdata['_space_group_symop_operation_xyz']
    elif '_symmetry_equiv_pos_as_xyz' in CIFdata.keys():
        operation_list = CIFdata['_symmetry_equiv_pos_as_xyz']
    else:
        raise Exception('No symmetry operations in .cif file')

    obj.operation_list = [item.split(',') for item in operation_list]

    # Atomic positions:
    atom_list = CIFdata['_atom_site_label']
    x_list = [error_strip(x) for x in CIFdata['_atom_site_fract_x']]
    y_list = [error_strip(y) for y in CIFdata['_atom_site_fract_y']]
    z_list = [error_strip(z) for z in CIFdata['_atom_site_fract_z']]

    if '_atom_site_occupancy' not in CIFdata.keys():
        occ_list = np.ones(len(x_list))
    else:
        occ_list = [error_strip(occ)
                    for occ in CIFdata['_atom_site_occupancy']]

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

    # We define display reflections, adding "i" if hexagonal:
    if crystal.hex == False:
        crystal.reflections_dis = crystal.reflections
    else:
        crystal.reflections_dis = [[((hkl[0], hkl[1], - hkl[0] - hkl[1], hkl[2])), I]
                                   for hkl, I in hkl_and_Int]
    return crystal.reflections_dis


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
    # main app

    app = QApplication(sys.argv)
    screen_GUI = Screen()  # we create the object

    screen_GUI.show()

    try:
        sys.exit(app.exec())  # execution
    except BaseException:
        print('Thanks for using inserexs.')
