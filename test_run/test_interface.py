# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 18:57:43 2018

@author: Kazantsev Andrey
"""

"""


cmd = "ptimer.exe 1133.tpo -i -m -d ptimer.cfg"

os.system(cmd)
"""
import os
import sys

from PyQt5.QtWidgets import (QMainWindow, QApplication, QLabel, QRadioButton,
                             QVBoxLayout, QFrame, QWidget, QGridLayout)
from PyQt5.QtCore import Qt

class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'Ptimer'
        self.left = 20
        self.top = 30
        self.width = 505
        self.height = 670
        self.out_dir = os.getcwd() + os.sep
        self.initUI()


    def initUI(self):

        '''
        Применение основных параметров окна
        '''
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)


        self.label_workstate = QLabel('Выбор режима работы')
        self.label_workstate.setObjectName('label_workstate')
        self.label_workstate.setStyleSheet(
            'QLabel#label_workstate {font: bold;}')
        self.workstate_chron = QRadioButton('Хронометрирование')
        self.workstate_toa = QRadioButton('Генерация МПИ')


        layout_workstate = QVBoxLayout()
        #layout.setContentsMargins(5, 5, 5, 5)
        layout_workstate.addWidget(self.label_workstate,
                                          0, Qt.AlignTop | Qt.AlignCenter)
        layout_workstate.addWidget(self.workstate_chron,
                                          0, Qt.AlignTop | Qt.AlignLeft)
        layout_workstate.addWidget(self.workstate_toa,
                                          0, Qt.AlignTop | Qt.AlignLeft)


        frame_workstate = QFrame()
        frame_workstate.setObjectName('frame_workstate')
        frame_workstate.setLayout(layout_workstate)
        frame_workstate.setStyleSheet(
                'QFrame#frame_workstate {border:1px solid black;}')

        hbox = QGridLayout()
        # addWidget(QWidget, row, column, rows, columns)
        hbox.addWidget(frame_workstate, 0, 0, 2, 4)

        main_panel = QWidget()
        main_panel.setLayout(hbox)

        self.setCentralWidget(main_panel)
        self.show()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())