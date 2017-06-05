import os
import sys
import datetime
from PyQt5 import QtCore, QtGui,QtWidgets
from scripts.secondary_windows import ExperimentWindow
import scripts.EcoHab
from scripts.ExperimentConfigFile import ExperimentConfigFile
import pickle
from scripts.experiments_info import smells, antenna_positions

__author__ = 'JanMaka'
__version__ = "0.1"
class NewMainWindow(QtWidgets.QMainWindow):
    def closeEvent(self, event):
        quit_msg = "Are you sure you want to exit the program?"
        reply = QtWidgets.QMessageBox.question(self, 'Message',
                         quit_msg, QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

        if reply == QtWidgets.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

class _SortProxyModel(QtCore.QSortFilterProxyModel):
    """Sorting proxy model that always places folders on top."""
    def __init__(self, model):
        super(_SortProxyModel,self).__init__(None)
        self.setSourceModel(model)

    def lessThan(self, left_index, right_index):

        left_var = left_index.data(QtCore.Qt.EditRole)
        right_var = right_index.data(QtCore.Qt.EditRole)

        try:
            return float(left_var) < float(right_var)
        except (ValueError, TypeError):
            pass

        try:
            result = sorted([str(left_var), str(right_var)],key=lambda x: datetime.datetime.strptime(x,'%d.%m.%y'))
            return result[0] <result[1]
        except (ValueError, TypeError):
            pass
        return left_var < right_var


class MainWindow(NewMainWindow):

    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setWindowTitle("EcoHAB Analyser " + __version__)


        #Main variables
        self.experiments_id = {}
        self.experiments = []
        self.key_items = []
        self.experiment_fields = ["ID","Name","Mouse Type", "Date"]
        self.protected = ["ID"]
        self.selected = []
        tab_widget = QtWidgets.QTabWidget()
        self.tab1 = QtWidgets.QWidget()
        tab_widget.addTab(self.tab1, 'Data Menager')
        self.tab1.setEnabled(True)
        self.tab2 = QtWidgets.QWidget()
        self.tab2.setEnabled(False)
        tab_widget.addTab(self.tab2, 'Analisys Tab')
        self.setCentralWidget(tab_widget)
        self.initUI()

    def initUI(self):
        databuttonsbox = QtWidgets.QHBoxLayout()
        self.all_exp =  QtWidgets.QPushButton("Select all")
        self.all_exp.clicked.connect(self.select_all)
        self.non_exp = QtWidgets.QPushButton("Unselect all")
        self.non_exp.clicked.connect(self.select_none)
        databuttonsbox.addWidget(self.all_exp)
        databuttonsbox.addWidget(self.non_exp)
        
        
        runbuttonsbox = QtWidgets.QHBoxLayout()
        self.run_pre =  QtWidgets.QPushButton("RUN PREPROCESSING")
        self.run_pre.clicked.connect(self.run_preprocessing)
        runbuttonsbox.addWidget(self.run_pre)
        
        self.databox = QtWidgets.QVBoxLayout()

        self.treeView = QtWidgets.QTreeView()
        self.treeView.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.treeView.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        #self.treeView.setSelectionBehaviour(QtWidgets.QTreeView.SelectRows)
        self.treeView.clicked.connect(self.exeriment_Clicked)
        self.treeView.doubleClicked.connect(self.exeriment_DoubleClicked)
        self.model = QtGui.QStandardItemModel()
        self.treeView.setModel(self.model)
        self.model.setHorizontalHeaderLabels(self.experiment_fields)
        for el in next(os.walk('RawData/'))[1]:
            if el not in self.experiments_id.keys():
                self.experiments_id[el] = len(self.experiments_id)
                new_experiment = {}
                new_experiment["Name"] = el
                new_experiment["Mouse Type"] = el[:5]
                new_experiment["Date"] = el[-8:]
                new_experiment["Path"] = el
                try:
                    new_experiment["Smells"] = smells[el]
                except:
                    new_experiment["Smells"] = None
                try:
                    new_experiment["Antennas"] = antenna_positions
                except:
                    new_experiment["Antennas"] = None
                new_experiment["Antennas"] = None
                self.experiments.append(new_experiment)
        for idx in range(len(self.experiments)):
            idx_item = QtGui.QStandardItem(str(idx))
            idx_item.setAccessibleText(str(idx))
            self.experiments[idx]["ID"] = str(idx)
            #idx_item.setCheckState(QtCore.Qt.Unchecked)
            #idx_item.setFlags(QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
            self.key_items.append([idx_item])
            for key_val in self.experiment_fields[1:]:
                self.key_items[idx].append(QtGui.QStandardItem(self.experiments[idx][key_val]))
                self.key_items[idx][-1].setAccessibleText(self.experiments[idx][key_val])
            self.model.appendRow(self.key_items[idx])
        self.treeView.setColumnWidth(0,60)
        self.proxy_model = _SortProxyModel(self.model)
        self.treeView.setModel(self.proxy_model)
        self.treeView.setSortingEnabled(True)

        self.databox.addWidget(self.treeView)
        self.databox.addLayout(databuttonsbox)
        self.databox.addLayout(runbuttonsbox)
        mainlayout = QtWidgets.QHBoxLayout()
        mainlayout.addLayout(self.databox)
        #mainlayout.addLayout(vbox)
        #self.canvas.setFixedSize(512, 424)
        #mainlayout.addWidget(self.canvas)
        #mainlayout.addLayout(self.control_lay)
        self.tab1.setLayout(mainlayout)


    def exeriment_Clicked(self, index):
        index = self.proxy_model.mapToSource(index)
        idx = index.row()
        if idx not in self.selected:
            for i in range(len(self.experiment_fields)):
                item = self.key_items[index.row()][i]
                item.setBackground(QtGui.QColor(211,211,211))
            self.selected.append(idx)
        else:
            for i in range(len(self.experiment_fields)):
                item = self.key_items[index.row()][i]
                item.setBackground(QtGui.QColor(255,255,255))
            self.selected.remove(idx)

    def exeriment_DoubleClicked(self, index):
        print  index.row()
        dialog = ExperimentWindow(index.row(),self)#item.accessibleText()))
        dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        dialog.exec_()
        for i in range(len(self.experiment_fields)):
            item = self.key_items[index.row()][i]
            item.setText(self.experiments[index.row()][self.experiment_fields[i]])

    def select_all(self):
        for idx in range(len(self.experiments_id)):
            for i in range(len(self.experiment_fields)):
                item = self.key_items[idx][i]
                item.setBackground(QtGui.QColor(211,211,211))
                self.selected = self.experiments_id.values()

    def select_none(self):
        for idx in range(len(self.experiments_id)):
            for i in range(len(self.experiment_fields)):
                item = self.key_items[idx][i]
                item.setBackground(QtGui.QColor(255,255,255))
                self.selected = []

    def run_preprocessing(self):
        print self.selected



if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    screen = QtWidgets.QDesktopWidget().screenGeometry()
    main = MainWindow()
    x,y = 960,540
    main.setGeometry((screen.width()-x)/2,(screen.height()-y)/2, x,y)
    main.show()

    sys.exit(app.exec_())
