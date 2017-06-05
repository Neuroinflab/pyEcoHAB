import os
import thread
import sys
import sip
import numpy as np
#sip.setapi('QString', 2)
#sip.setapi('QVariant', 2)
from PyQt5 import QtCore, QtGui, QtWidgets

class ExperimentWindow(QtWidgets.QDialog):
    
    def __init__(self,idx,main,parent=None):
        super(ExperimentWindow, self).__init__(parent)
        self.main = main
        self.idx = idx
        self.cb = []
        self.setWindowTitle(self.main.experiments[idx]["Name"]) 
        #self.experiment = experiment
        l_vbox = QtWidgets.QVBoxLayout()
        e_vbox = QtWidgets.QVBoxLayout()
        for ix, key in enumerate(self.main.experiment_fields):
                    label = QtWidgets.QLabel(key)
                    l_vbox.addWidget(label)
                    self.cb.append(QtWidgets.QLineEdit())
                    e_vbox.addWidget(self.cb[ix])
                    self.cb[ix].setText(self.main.experiments[self.idx][key])
                    if (key in self.main.protected):
                        self.cb[ix].setEnabled(False)
    
        dialogbox = QtWidgets.QHBoxLayout()
        dialogbox.addLayout(l_vbox)
        dialogbox.addLayout(e_vbox)
        self.cbutton = QtWidgets.QPushButton('SAVE')
        self.cbutton.clicked.connect(self.save_file)
        main_box = QtWidgets.QVBoxLayout()
        main_box.addLayout(dialogbox)
        main_box.addWidget(self.cbutton)
        self.setLayout(main_box)
    def save_file(self):
        for ix, key in enumerate(self.main.experiment_fields):
            self.main.experiments[self.idx][key] = str(self.cb[ix].text())
        self.close()

