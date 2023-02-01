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

Add-on for the generation of a cif file.
"""

import os
import sys
from PyQt5 import QtGui, uic # v5.14.4
from PyQt5.QtWidgets import QApplication, QTableWidgetItem  # v5.14.4
from pathlib import Path
import CifFile # v4.4.5
import pyxtal # v0.5.5

home_cwd = os.getcwd()  # current cwd
Ui_path = Path(home_cwd, 'builder.ui')
Ui_myWindow_base, Ui_myWindow_form = uic.loadUiType(Ui_path)

class builderScreen (Ui_myWindow_base, Ui_myWindow_form):
    """Main GUI."""

    def __init__(self, parent=None, *args, **kwargs):
        super(builderScreen, self).__init__(parent=parent, *args, **kwargs)
        self.setupUi(self)
        # self.setWindowIcon(QtGui.QIcon('Icon.ico'))  # for the windows icon
        self.racine_directory = home_cwd
        self.setWindowIcon(QtGui.QIcon('icon.ico'))
        self.addButton.clicked.connect(self.add_row)
        self.deleteButton.clicked.connect(self.delete_row)
        self.generateButton.clicked.connect(self.generate)

    def add_row(self):
        """Add new row to the table containing the atomic information."""
        rowPosition = self.tableWidget.rowCount()
        self.tableWidget.insertRow(rowPosition) #insert new row

        # new row format should follow the following grammar:
        new_row = (0, 0, 0, 1)

        for column in range(1, 5):
            self.tableWidget.setItem(rowPosition, column,
                                     QTableWidgetItem(str(new_row[column - 1])))

    def delete_row(self):
        """Delete new row to the table containing the atomic information."""
        rowPosition = self.tableWidget.rowCount()
        self.tableWidget.removeRow(rowPosition-1) #insert new row

    def get_data(self):
        """Get data from the input values."""
        crystal = Crystal()

        crystal.name = self.nameLine.text()

        crystal.a = self.aLine.text()
        crystal.b = self.bLine.text()
        crystal.c = self.cLine.text()

        crystal.alpha = self.alphaLine.text()
        crystal.beta = self.betaLine.text()
        crystal.gamma = self.gammaLine.text()

        crystal.sg = self.sgLine.text()

        crystal.atoms = [self.tableWidget.item(row, 0).text()
                     for row in range(self.tableWidget.rowCount())]
        crystal.x = [self.tableWidget.item(row, 1).text()
                     for row in range(self.tableWidget.rowCount())]
        crystal.y = [self.tableWidget.item(row, 2).text()
                     for row in range(self.tableWidget.rowCount())]
        crystal.z = [self.tableWidget.item(row, 3).text()
                     for row in range(self.tableWidget.rowCount())]
        crystal.occ = [self.tableWidget.item(row, 4).text()
                     for row in range(self.tableWidget.rowCount())]

        return crystal

    def generate(self):
        """Generate .cif file from input values."""
        obj = self.get_data()
        cf = CifFile.CifFile()

        # creation of new block
        crystal_info = CifFile.CifBlock()
        cf[obj.name] = crystal_info

        crystal_info['_cell_length_a'] = obj.a
        crystal_info['_cell_length_b'] = obj.b
        crystal_info['_cell_length_c'] = obj.c

        crystal_info['_cell_angle_alpha'] = obj.alpha
        crystal_info['_cell_angle_beta'] = obj.beta
        crystal_info['_cell_angle_gamma'] = obj.gamma

        crystal_info['_symmetry_space_group_name_H-M'] = obj.sg

        spacegroup = pyxtal.symmetry.Group(obj.sg)
        crystal_info['_space_group_symop_operation_xyz'] = str(spacegroup.Wyckoff_positions[0])
        crystal_info.CreateLoop(['_space_group_symop_operation_xyz'])

        crystal_info['_atom_site_label'] = list(obj.atoms)
        crystal_info['_atom_site_fract_x'] = list(obj.x)
        crystal_info['_atom_site_fract_y'] = list(obj.y)
        crystal_info['_atom_site_fract_z'] = list(obj.z)
        crystal_info['_atom_site_fract_occupancy'] = list(obj.occ)

        crystal_info.CreateLoop(['_atom_site_label', '_atom_site_fract_x',
                                 '_atom_site_fract_y', '_atom_site_fract_z',
                                 '_atom_site_fract_occupancy'])

        outfile = open(str(obj.name) + '.cif', 'w')
        outfile.write(cf.WriteOut())
        return(cf)

class Crystal:
    """Object that will keep all the crystal information."""

    def __init__(self):
        pass

if __name__ == 'cif builder':  # only executes the below code if  main
    app = QApplication(sys.argv)
    bscreen_GUI = builderScreen()  # we create the object
    bscreen_GUI.show()

    try:
        sys.exit(app.exec())  # execution
    except BaseException:
        print('Exit')
