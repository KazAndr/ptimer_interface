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
                             QVBoxLayout, QFrame, QWidget, QGridLayout,
                             QCheckBox)
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


        '''
        Инициализация первой рабочей области с выбором режима работы
        '''
        self.label_workstate = QLabel('Выбор режима работы')
        self.label_workstate.setObjectName('label_workstate')
        self.label_workstate.setStyleSheet(
            'QLabel#label_workstate {font: bold;}')
        self.workstate_chron = QRadioButton('Хронометрирование')
        self.workstate_toa = QRadioButton('Генерация МПИ')

        '''
        Инициализация второй рабочей области с выбором опций
        '''
        self.label_option = QLabel('Выбор опций')
        self.label_option.setObjectName('label_option')
        self.label_option.setStyleSheet(
            'QLabel#label_option {font: bold;}')
        self.opt_param = QCheckBox('Не показывать начальные параметры')
        self.opt_aver = QCheckBox('Не вычитать среднее значение')
        self.opt_inp = QCheckBox(
                'Вывод уточненных параметров в файл _tim.inp')
        self.opt_apr = QCheckBox(
                'Вывод точной корреляционной матрицы в файл _tim.apr')
        self.opt_cfg = QCheckBox('Выбрать концигурационный файл')
        '''
        Компановка объектов первой рабочей области
        '''
        layout_workstate = QVBoxLayout()
        #layout.setContentsMargins(5, 5, 5, 5)
        layout_workstate.addWidget(self.label_workstate,
                                          0, Qt.AlignTop | Qt.AlignCenter)
        layout_workstate.addWidget(self.workstate_chron,
                                          0, Qt.AlignTop | Qt.AlignLeft)
        layout_workstate.addWidget(self.workstate_toa,
                                          0, Qt.AlignTop | Qt.AlignLeft)

        '''
        Компановка объектов второй рабочей области
        '''
        layout_option = QVBoxLayout()
        #layout.setContentsMargins(5, 5, 5, 5)
        layout_option.addWidget(self.label_option,
                                          0, Qt.AlignTop | Qt.AlignCenter)
        layout_option.addWidget(self.opt_param,
                                          0, Qt.AlignTop | Qt.AlignLeft)
        layout_option.addWidget(self.opt_aver,
                                          0, Qt.AlignTop | Qt.AlignLeft)
        layout_option.addWidget(self.opt_inp,
                                          0, Qt.AlignTop | Qt.AlignLeft)
        layout_option.addWidget(self.opt_apr,
                                          0, Qt.AlignTop | Qt.AlignLeft)
        layout_option.addWidget(self.opt_cfg,
                                          0, Qt.AlignTop | Qt.AlignLeft)


        '''
        Установка значений по умолчанию
        '''
        self.workstate_chron.setChecked(True)

        '''
        Подключение функций к кнопкам
        '''
        self.workstate_chron.toggled.connect(lambda:self.check_active())
        #self.workstate_toa.toggled.connect(lambda:self.check_active())
        '''
        Инициализация фрейма для первой рабочей области
        Упаковка рабочей области
        '''
        frame_workstate = QFrame()
        frame_workstate.setObjectName('frame_workstate')
        frame_workstate.setLayout(layout_workstate)
        frame_workstate.setStyleSheet(
                'QFrame#frame_workstate {border:1px solid black;}')

        '''
        Инициализация фрейма для первой рабочей области
        Упаковка рабочей области
        '''
        frame_option = QFrame()
        frame_option.setObjectName('frame_option')
        frame_option.setLayout(layout_option)
        frame_option.setStyleSheet(
                'QFrame#frame_option {border:1px solid black;}')


        '''
        Упаковка фреймов
        '''
        hbox = QGridLayout()
        # addWidget(QWidget, row, column, rows, columns)
        hbox.addWidget(frame_workstate, 0, 0)
        hbox.addWidget(frame_option, 1, 0)

        main_panel = QWidget()
        main_panel.setLayout(hbox)

        self.setCentralWidget(main_panel)
        self.check_active()
        self.show()

        '''
        Инициализация рабочих функций
        '''
    def check_active(self):
        if self.workstate_chron.isChecked():
            print('chron')
            self.opt_param.setEnabled(True)
            self.opt_aver.setEnabled(True)
            self.opt_inp.setEnabled(True)
            self.opt_apr.setEnabled(True)
        else:
            print('nechron')
            self.opt_param.setEnabled(False)
            self.opt_aver.setEnabled(False)
            self.opt_inp.setEnabled(False)
            self.opt_apr.setEnabled(False)





if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())