# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:05:20 2024

@author: joylove5
"""
### Depedentcies
import pandas as pd
import numpy as np
import sys
#import os
from PyQt5.QtCore import Qt, QCoreApplication, QtCore
from PyQt5.QtGui import QColor, QPalette, QFont
from PyQt5.QtWidgets import (QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QTabWidget,
                             QTextEdit, QDialog, QLabel, QHBoxLayout, QLineEdit, QFileDialog)

############################# Functions for functionality #####################
### For Tab1: HetSNP List Generation
def data_prep(filename):
    name = filename.split(".vcf")[0]
    df = pd.read_csv(filename, sep="\t", comment="#", header=None)
    df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", name]
    df["ID"] = df["CHROM"] + "_" + df["POS"].astype(str)
    return df

def semi_join(data_Set1, data_Set2, joiner):
    semi = data_Set1.merge(data_Set2, on=joiner)
    data_Set1[joiner].isin(data_Set2[joiner])
    semi = data_Set1.merge(data_Set2, on=joiner)
    new_semi = data_Set1[data_Set1[joiner].isin(semi[joiner])]
    pd.DataFrame(new_semi)
    return new_semi

def filter_positions(df1, unreliable_ranges):
    exclude_mask = False
    for index, row in unreliable_ranges.iterrows():
        chromosome = row['Chromosome']
        start_position = row['Region_Start']
        end_position = row['Region_End']
        mask = (df1['CHROM'] == f'chr{chromosome}') & (df1['POS'].between(start_position, end_position))
        exclude_mask |= mask
    filtered_df = df1[~exclude_mask]
    return filtered_df

def add_header(_filtered_df, date, wd, SNP_List_Name):
    header = """##fileformat=VCFv4.1
##fileDate=""" + date + """
##source=CLC Genomics Workbench 23.0.4 build 20230515014922
##reference=/CLC_Ddrive/S288c reference genome/Reference (Genome)
##contig=<ID=chr1,length=230218>
##contig=<ID=chr2,length=813184>
##contig=<ID=chr3,length=316620>
##contig=<ID=chr4,length=1531933>
##contig=<ID=chr5,length=576874>
##contig=<ID=chr6,length=270161>
##contig=<ID=chr7,length=1090940>
##contig=<ID=chr8,length=562643>
##contig=<ID=chr9,length=439888>
##contig=<ID=chr10,length=745751>
##contig=<ID=chr11,length=666816>
##contig=<ID=chr12,length=1078177>
##contig=<ID=chr13,length=924431>
##contig=<ID=chr14,length=784333>
##contig=<ID=chr15,length=1091291>
##contig=<ID=chr16,length=948066>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total number of filtered reads used to call variants at this position in this sample.">
##FORMAT=<ID=CLCAD2,Number=R,Type=Integer,Description="Allele depth, number of filtered reads supporting the alleles, ordered as listed in REF and ALT.">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HetSNP_Variants
"""
    output_VCF = wd + SNP_List_Name
    with open(output_VCF, 'a') as vcf:
        vcf.write(header)
    final_output = _filtered_df.to_csv(output_VCF, sep="\t", mode='a', header=False, index=False)
    return final_output

def snp_main(diploid_vcf, SK1_vcf, unreliables_csv, date, wd, SNP_List_Name):
    dip_df = data_prep(diploid_vcf)
    sk1_df = data_prep(SK1_vcf)
    match_df = semi_join(dip_df, sk1_df, "ID")
    unreliable_df = pd.read_csv(unreliables_csv)
    filtered_df = filter_positions(match_df, unreliable_df)
    output_SNP_VCF = add_header(filtered_df, date, wd, SNP_List_Name)
    return output_SNP_VCF

### For Tab2: Haplotype-Aware Tables
def Make_haplotype_aware(clone_table, clone_name, path, minimum_coverage):
    # Taking in Variants  
    variant_df = pd.read_csv(clone_table, sep="\t", comment="#", header=None)
    variant_df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", clone_name]
    variant_df["ID"] = variant_df["CHROM"] + "_" + variant_df["POS"].astype(str)
    variant_df = variant_df.drop(variant_df[variant_df["ID"] == "chr9_431988"].index)
    variant_df = variant_df.drop(variant_df[variant_df["ID"] == "chr10_416768"].index)
    
    # Building Clone Variant Table
    new_df = pd.DataFrame()
    new_df["ID"] = variant_df["ID"]
    new_df["S288c_Parent_Allele"] = variant_df["REF"]
    new_df["SK1_Parent_Allele"] = variant_df["ALT"]
    new_df[clone_name + "_Variant_Allele"] = variant_df["ALT"]
    
    # Build coverage and frequency columns
    co_cov = [word.strip().split(',')[1] for word in variant_df[clone_name]]
    count, coverage = zip(*((int(num.strip().split(':')[0]), int(num.strip().split(':')[1])) for num in co_cov))
    new_df[clone_name + "_Variant_Coverage"] = coverage
    frequency = np.round([(count[i] / coverage[i]) * 100 if coverage[i] != 0 else 0 for i in range(len(count))], decimals=2)
    new_df[clone_name + "_Variant_Frequency"] = frequency

    # Genotype Determination
    zygosity = ["Homozygous for SK1" if percent > 89 else "Heterozygous" if 11 <= percent <= 89 else "Homozygous for S288c" for percent in new_df[clone_name + "_Variant_Frequency"]]
    zygosity = ["No Call Due to Low Coverage" if cov <= minimum_coverage else zygosity[i] for i, cov in enumerate(new_df[clone_name + "_Variant_Coverage"])]
    new_df[clone_name + "_Genotype"] = zygosity
    
    # Making Output CSV file
    new_df.to_csv(path + "/" + clone_name + "_Haplotype_Aware_Table(VCF).csv", index=False)
    return new_df

def Create_master_variant_table(clone_list_file, path, output_path, coverage_minimum):
    master_df = pd.DataFrame()
    with open(clone_list_file, "r") as file:
        clone_tables = [line.strip() for line in file]
    
    for clone_table in clone_tables:
        new_df = Make_haplotype_aware(clone_table, clone_table.split(' (')[0], path, coverage_minimum)
        if master_df.empty:
            master_df = new_df[["ID", "S288c_Parent_Allele", "SK1_Parent_Allele"]].copy()
        columns_to_merge = [
            "ID",
            f"{clone_table.split(' (')[0]}_Variant_Coverage",
            f"{clone_table.split(' (')[0]}_Variant_Frequency",
            f"{clone_table.split(' (')[0]}_Genotype"
        ]
        master_df = pd.merge(master_df, new_df[columns_to_merge], on="ID", how="outer")
    
    master_df.to_csv(output_path, index=False)
    return master_df


### For Tab3: Further Visualization of LOH Chromosomes



################ Classes and Functions for Graphical Interface ################
### Set-Up Window
class SetupDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Setup")
        self.setGeometry(100, 100, 300, 150)

        self.label = QLabel("Do you have a High Confidence HetSNP List for yeast mated between S288c and SK1?")
        self.label.move(100, 222) 
        self.label.setFont(QFont("Arial", 14))
        self.button1 = QPushButton("Generate HetSNP list.")
        self.button2 = QPushButton("Yes, generate haplotype-aware tables.")
        
        self.button1.setStyleSheet("background-color: darkgrey")
        self.button2.setStyleSheet("background-color: darkgrey")
        
        self.button1.clicked.connect(self.start_tab1)
        self.button2.clicked.connect(self.start_tab2)
        
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.label)
        
        self.button_layout = QHBoxLayout()
        self.button_layout.addWidget(self.button1)
        self.button_layout.addWidget(self.button2)
        
        self.layout.addLayout(self.button_layout)
        self.setLayout(self.layout)
        
        self.selected_tab = 0

    def start_tab1(self):
        self.selected_tab = 0
        self.accept()

    def start_tab2(self):
        self.selected_tab = 1
        self.accept()

### Diplay Pandas Dataframes
class PandasModel(QtCore.QAbstractTableModel):
    def __init__(self, data, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return len(self._data.values)

    def columnCount(self, parent=None):
        return self._data.columns.size

    def data(self, index, role=Qt.DisplayRole):
        if index.isValid():
            if role == Qt.DisplayRole:
                return QtCore.QVariant(str(
                    self._data.iloc[index.row()][index.column()]))
        return QtCore.QVariant()

### For Tab1: HetSNP List Generation
class Tab1(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()

        self.wd_label = QLabel("Working Directory:")
        self.wd_edit = QLineEdit()
        self.wd_button = QPushButton("Browse")
        self.wd_button.clicked.connect(self.browse_wd)

        self.diploid_vcf_label = QLabel("Diploid VCF File:")
        self.diploid_vcf_edit = QLineEdit()
        self.diploid_vcf_button = QPushButton("Browse")
        self.diploid_vcf_button.clicked.connect(self.browse_diploid_vcf)

        self.SK1_vcf_label = QLabel("SK1 Parent VCF File:")
        self.SK1_vcf_edit = QLineEdit()
        self.SK1_vcf_button = QPushButton("Browse")
        self.SK1_vcf_button.clicked.connect(self.browse_SK1_vcf)

        self.unreliables_csv_label = QLabel("Unreliable Regions CSV File:")
        self.unreliables_csv_edit = QLineEdit()
        self.unreliables_csv_button = QPushButton("Browse")
        self.unreliables_csv_button.clicked.connect(self.browse_unreliables_csv)

        self.date_label = QLabel("Date (YYYYMMDD):")
        self.date_edit = QLineEdit()

        self.snp_list_name_label = QLabel("SNP List Output Name:")
        self.snp_list_name_edit = QLineEdit()

        self.button = QPushButton("Generate HetSNP List")
        
        self.output_label = QLabel("Output:")
        self.text_edit = QTextEdit()

        self.layout.addWidget(self.wd_label)
        self.layout.addWidget(self.wd_edit)
        self.layout.addWidget(self.wd_button)
        self.layout.addWidget(self.diploid_vcf_label)
        self.layout.addWidget(self.diploid_vcf_edit)
        self.layout.addWidget(self.diploid_vcf_button)
        self.layout.addWidget(self.SK1_vcf_label)
        self.layout.addWidget(self.SK1_vcf_edit)
        self.layout.addWidget(self.SK1_vcf_button)
        self.layout.addWidget(self.unreliables_csv_label)
        self.layout.addWidget(self.unreliables_csv_edit)
        self.layout.addWidget(self.unreliables_csv_button)
        self.layout.addWidget(self.date_label)
        self.layout.addWidget(self.date_edit)
        self.layout.addWidget(self.snp_list_name_label)
        self.layout.addWidget(self.snp_list_name_edit)
        self.layout.addWidget(self.button)
        self.layout.addWidget(self.output_label)
        self.layout.addWidget(self.text_edit)
        self.setLayout(self.layout)
        
        self.button.clicked.connect(self.tab1_function)

    def browse_wd(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Working Directory")
        if directory:
            self.wd_edit.setText(directory)

    def browse_diploid_vcf(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select Diploid VCF File", "", "VCF Files (*.vcf);;All Files (*)")
        if file:
            self.diploid_vcf_edit.setText(file)

    def browse_SK1_vcf(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select SK1 VCF File", "", "VCF Files (*.vcf);;All Files (*)")
        if file:
            self.SK1_vcf_edit.setText(file)

    def browse_unreliables_csv(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select Unreliable Regions CSV File", "", "CSV Files (*.csv);;All Files (*)")
        if file:
            self.unreliables_csv_edit.setText(file)

    def tab1_function(self):
        wd = self.wd_edit.text()
        diploid_vcf = self.diploid_vcf_edit.text()
        SK1_vcf = self.SK1_vcf_edit.text()
        unreliables_csv = self.unreliables_csv_edit.text()
        date = self.date_edit.text()
        SNP_List_Name = self.snp_list_name_edit.text()

        self.text_edit.append("Executing Tab 1 Function with inputs:")
        self.text_edit.append(f"Working Directory: {wd}")
        self.text_edit.append(f"Diploid VCF File: {diploid_vcf}")
        self.text_edit.append(f"SK1 VCF File: {SK1_vcf}")
        self.text_edit.append(f"Unreliable Regions CSV File: {unreliables_csv}")
        self.text_edit.append(f"Date: {date}")
        self.text_edit.append(f"SNP List Output Name: {SNP_List_Name}")
        
        # Call the main function
        output_SNP_VCF = snp_main(diploid_vcf, SK1_vcf, unreliables_csv, date, wd, SNP_List_Name)
        self.text_edit.append(f"HetSNP List Generated: {output_SNP_VCF}")

### For Tab2: Haplotype-Aware Tables
class Tab2(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()

        self.wd_label = QLabel("Working Directory:")
        self.wd_edit = QLineEdit()
        self.wd_button = QPushButton("Browse")
        self.wd_button.clicked.connect(self.browse_wd)

        self.output_dir_label = QLabel("Output Directory:")
        self.output_dir_edit = QLineEdit()
        self.output_dir_button = QPushButton("Browse")
        self.output_dir_button.clicked.connect(self.browse_output_dir)

        self.clone_tables_list_label = QLabel("Clone Tables List:")
        self.clone_tables_list_edit = QLineEdit()
        self.clone_tables_list_button = QPushButton("Browse")
        self.clone_tables_list_button.clicked.connect(self.browse_clone_table_list)
        
        self.coverage_minimum_label = QLabel("Minimum Coverage:")
        self.coverage_minimum_edit = QLineEdit()
        
        self.master_haplotype_table_label = QLabel("Master Haplotype Table Name:")
        self.master_haplotype_table_edit = QLineEdit()

        self.button = QPushButton("Generate Haplotype-Aware Tables")
        
        self.output_label = QLabel("Output:")
        self.text_edit = QTextEdit()

        self.layout.addWidget(self.wd_label)
        self.layout.addWidget(self.wd_edit)
        self.layout.addWidget(self.wd_button)
        self.layout.addWidget(self.output_dir_label)
        self.layout.addWidget(self.output_dir_edit)
        self.layout.addWidget(self.output_dir_button)
        self.layout.addWidget(self.clone_tables_list_label)
        self.layout.addWidget(self.clone_tables_list_edit)
        self.layout.addWidget(self.clone_tables_list_button)
        self.layout.addWidget(self.coverage_minimum_label)
        self.layout.addWidget(self.coverage_minimum_edit)
        self.layout.addWidget(self.master_haplotype_table_label)
        self.layout.addWidget(self.master_haplotype_table_edit)
        self.layout.addWidget(self.button)
        self.layout.addWidget(self.output_label)
        self.layout.addWidget(self.text_edit)
        self.setLayout(self.layout)

        self.button.clicked.connect(self.tab2_function)

    def browse_wd(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Working Directory")
        if directory:
            self.wd_edit.setText(directory)

    def browse_output_dir(self):
        directory = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if directory:
            self.output_dir_edit.setText(directory)
    
    def browse_clone_table_list(self):
        file, _ = QFileDialog.getOpenFileName(self, "Select File with List of Clone Table Names", "", "TXT Files (*.txt);;All Files (*)")
        if file:
            self.clone_tables_list_edit.setText(file)
        
    def tab2_function(self):
        wd = self.wd_edit.text()
        output_dir = self.output_dir_edit.text()
        clone_tables_list = self.clone_tables_list_edit.text()
        coverage_minimum = int(self.coverage_minimum_edit.text())
        master_haplotype_table = self.master_haplotype_table_edit.text()

        """ 
        self.text_edit.append("Executing Tab 2 Function with inputs:")
        self.text_edit.append(f"Working Directory: {wd}")
        self.text_edit.append(f"Output Directory: {output_dir}")
        self.text_edit.append(f"Clone Tables List: {clone_tables_list}")
        self.text_edit.append(f"Minimum Coverage: {coverage_minimum}")
        self.text_edit.append(f"Master Haplotype Table: {master_haplotype_table}")
        """
        # Call Create_master_variant_table function
        master_df = Create_master_variant_table(clone_tables_list, wd, output_dir + "/" + master_haplotype_table, coverage_minimum)
        self.text_edit.append(f"Master Haplotype Table Created: {output_dir}/{master_haplotype_table}")
        
        master_df_str = master_df.to_string()
        self.text_edit.append(master_df_str)
    
### For Tab3: Haplotype-Aware Tables

### Main Window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LOH Caller")
        self.setGeometry(100, 100, 800, 1000)
        
        self.apply_styles()
        
        self.setup_dialog = SetupDialog()
        if self.setup_dialog.exec_() == QDialog.Accepted:
            self.selected_tab = self.setup_dialog.selected_tab
            self.setup_ui()
            
    ### Making Window Dark-Mode
    def apply_styles(self):
        dark_palette = QPalette()

        # Base colors
        dark_palette.setColor(QPalette.Window, QColor(45, 45, 45))
        dark_palette.setColor(QPalette.WindowText, Qt.white)
        dark_palette.setColor(QPalette.Base, QColor(30, 30, 30))
        dark_palette.setColor(QPalette.AlternateBase, QColor(45, 45, 45))
        dark_palette.setColor(QPalette.ToolTipBase, Qt.white)
        dark_palette.setColor(QPalette.ToolTipText, Qt.white)
        dark_palette.setColor(QPalette.Text, Qt.black)
        dark_palette.setColor(QPalette.Button, QColor(70, 70, 70))
        dark_palette.setColor(QPalette.ButtonText, Qt.white)
        dark_palette.setColor(QPalette.BrightText, Qt.red)

        # Disabled colors
        dark_palette.setColor(QPalette.Disabled, QPalette.Text, QColor(0, 0, 0))
        dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(127, 127, 127))

        # Highlight colors
        dark_palette.setColor(QPalette.Highlight, QColor(142, 45, 197).lighter())
        dark_palette.setColor(QPalette.HighlightedText, Qt.black)

        app.setPalette(dark_palette)

        self.setStyleSheet("""
            QMainWindow {
                background-color: #2D2D2D;
            }
            QTabWidget::pane {
                border: 1px solid darkgray;
                top:-1px; 
                background: rgb(30, 30, 30);; 
            } 

            QTabBar::tab {
                background: rgb(30, 30, 30); 
                border: 1px solid darkgray; 
                padding: 15px;
                } 

            QTabBar::tab:selected { 
                background: rgb(20, 20, 20); 
                margin-bottom: -1px;     
            }
            QLabel {
                color: white;
            }
            QLineEdit {
                background-color: #3E3E3E;
                color: white;
                border: 1px solid #1D1D1D;
                padding: 5px;
            }
            QTextEdit {
                background-color: #3E3E3E;
                color: white;
                border: 1px solid #1D1D1D;
                padding: 5px;
            }
            QPushButton {
                background-color: #4B0082;
                color: white;
                border: none;
                padding: 10px;
            }
            QPushButton:hover {
                background-color: #551A8B;
            }
            QPushButton:pressed {
                background-color: #3E0066;
            }
        """)

    def setup_ui(self):
        self.tabs = QTabWidget()
        self.tab1 = Tab1()
        self.tab2 = Tab2()
        
        self.tabs.addTab(self.tab1, "Generate HetSNP List")
        self.tabs.addTab(self.tab2, "Haplotype-Aware Tables Genotype Calls")
        self.tabs.setCurrentIndex(self.selected_tab)
        self.setCentralWidget(self.tabs)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())