#!/usr/bin/env python3
# This Python file uses the following encoding: utf-8
import sys
import os
import requests
import csv
import pprint

import PyQt5

from PyQt5.QtCore import *

from PySide2 import QtCore
from PySide2.QtWidgets import QApplication, QWidget, QTableWidget, QPushButton, QFileDialog, QListWidget
from PySide2.QtCore import QFile
from PySide2.QtUiTools import QUiLoader

class ProteoformWarehouse(QWidget):
    def __init__(self):
        super(ProteoformWarehouse, self).__init__()
        self.uisetup()
        self.gui.filenamebox.setReadOnly(True)
        self.h = 1

        self.gui.inputFileButton.clicked.connect(self.on_inputfile_clicked)
        self.gui.addAccession.clicked.connect(self.on_addAccession_clicked)
        self.gui.deleteAccession.clicked.connect(self.on_deleteAccession_clicked)
        self.gui.allAccessionGenerate.clicked.connect(self.runscript)
        self.gui.selectAccessionGenerate.clicked.connect(self.runmodifiedscript)
        self.gui.combinedxmlButton.clicked.connect(self.combinedxml)
        self.gui.combinedxmltable.clicked.connect(self.tablecombinedxml)



    def uisetup(self):
        loader = QUiLoader()
        path = os.path.join(os.path.dirname(__file__), "form.ui")
        ui_file = QFile(path)
        ui_file.reset()
        ui_file.open(QFile.ReadOnly)
        self.gui = loader.load(ui_file, self)
        ui_file.close()


    def on_inputfile_clicked(self):

        filename, _ = QFileDialog.getOpenFileName(self, "Select Input File",
                                                          os.getenv('HOME'), "*csv")

        self.gui.filenamebox.setReadOnly(False)
        self.gui.filenamebox.setText(str(filename))
        self.gui.filenamebox.setReadOnly(True)

    def on_addAccession_clicked(self):
        g = self.gui.listWidget.count()
        self.gui.listWidget.insertItem(g,"<Enter Accession Number Here>")
        entry = self.gui.listWidget.item(g)
        entry.setFlags(entry.flags() | QtCore.Qt.ItemIsEditable)

    def on_deleteAccession_clicked(self):
        if (self.gui.listWidget.count() > 1):
            self.gui.listWidget.takeItem(self.gui.listWidget.count()-1)


    def runscript(self):
        os.system('python3 updatedpro.py')

    def runmodifiedscript(self):
        command = 'python3 updatedpro.py '
        numfiles = self.gui.listWidget.count()
        for x in range(1,numfiles):
            placeholder = self.gui.listWidget.item(x)
            command = command + str(placeholder.text()) + ' '
        
        os.system(command)

    def combinedxml(self):
        goal = 'python3 combinedpro.py'
        os.system(goal)

    def tablecombinedxml(self):
        command = 'python3 combinedpro.py '
        numfiles = self.gui.listWidget.count()
        for x in range(1,numfiles):
            placeholder = self.gui.listWidget.item(x)
            command = command + str(placeholder.text()) + ' '
        
        os.system(command)



if __name__ == "__main__":
    app = QApplication([])
    widget = ProteoformWarehouse()
    widget.show()
    sys.exit(app.exec_())
