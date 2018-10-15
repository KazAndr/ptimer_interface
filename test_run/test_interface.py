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
                             QCheckBox, QPushButton, QFileDialog, QTextEdit,
                             QSizePolicy)
from PyQt5.QtCore import Qt, QRect


class ShowParam(QWidget):

    def __init__(self, title, file, parent=None):
        super(ShowParam, self).__init__()

        self.param_label = QTextEdit(self)
        self.title = title
        self.file = file
        self.left = 100
        self.top = 100
        self.width = 435
        self.height = 480
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        with open(self.file, 'r') as f:
            text = f.readlines()

        self.param_label.setText(str(text))
        self.param_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.param_label.setAlignment(Qt.AlignLeft)
        self.param_label.setGeometry(QRect(20, 20, 400, 480))

        workstate = QVBoxLayout()
        workstate.addWidget(self.param_label,
                                          0, Qt.AlignTop | Qt.AlignCenter)

class App(QMainWindow):

    def __init__(self):
        super().__init__()
        self.title = 'Ptimer'
        self.left = 20
        self.top = 30
        self.width = 300
        self.height = 370
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
        Кнопка, запускающая обработку
        '''
        self.procces_button = QPushButton('Начать вычисление')
        self.procces_button.setObjectName('self.procces_button')
        self.procces_button.setStyleSheet(
            'QPushButton#self.procces_button {font: bold; background-color: gray;}')

        '''
        Установка значений по умолчанию
        '''
        self.workstate_chron.setChecked(True)
        self.statusBar().showMessage('Готов')
        '''
        Подключение функций к кнопкам
        '''
        self.workstate_chron.toggled.connect(lambda:self.check_active())
        self.procces_button.clicked.connect(self.processing)
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
        hbox.addWidget(frame_workstate, 0, 0, 2, 2)
        hbox.addWidget(frame_option, 3, 0, 3, 2)
        hbox.addWidget(self.procces_button, 7, 1, 1, 1)

        main_panel = QWidget()
        main_panel.setLayout(hbox)
        main_panel.setMaximumSize(self.width, self.height)

        self.setCentralWidget(main_panel)
        self.check_active()
        self.show()

        '''
        Инициализация рабочих функций
        '''
    def check_active(self):
        if self.workstate_chron.isChecked():
            self.opt_param.setEnabled(True)
            self.opt_aver.setEnabled(True)
            self.opt_inp.setEnabled(True)
            self.opt_apr.setEnabled(True)
        else:
            self.opt_param.setEnabled(False)
            self.opt_aver.setEnabled(False)
            self.opt_inp.setEnabled(False)
            self.opt_apr.setEnabled(False)

    def processing(self):
        cmd = 'ptimer.exe '

        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog

        self.filePar, _ = QFileDialog.getOpenFileName(
                    self,'Выберите файл параметров',
                    '', 'Все файлы (*)',
                    options=options)


        if self.opt_param.isChecked():
            cmd += '-c '
        if self.opt_aver.isChecked():
            cmd += '-r '
        if self.opt_inp.isChecked():
            cmd += '-i '
        if self.opt_apr.isChecked():
            cmd += '-m '
        if self.opt_cfg.isChecked():
            self.fileConfig, _ = QFileDialog.getOpenFileName(
                    self,'Выберите конфигурационный файл',
                    '', 'Файл конфигурации (*.cfg);; Все файлы (*)',
                    options=options)
            cmd += '-d' + self.fileConfig + ' '

        cmd += self.filePar

        self.statusBar().showMessage(
                "Выполняется команда " + cmd)

        if not self.opt_param.isChecked():
            ShowParam('Параметры по умолчанию', self.filePar).show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())