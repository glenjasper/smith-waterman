#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import time
import argparse
import textwrap
import traceback
from collections import OrderedDict

def menu(args):
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
                                     description = "Implementation of the Smith–Waterman algorithm",
                                     epilog = textwrap.dedent("""\
         Examples of alignment:
           For amino acid sequences
             python %(prog)s -t aa -f sequences.fa -sm PAM250 -gap -1

           For nucleotide sequences
             python %(prog)s -t nt -f sequences.fa -m 2 -mi -1 -gap -2

         Thank you!"""))
    parser.add_argument("-t", "--type", choices = oalig.ARRAY_TYPE_SEQ, required = True, type = str.lower, help = oalig.mode_information(oalig.ARRAY_TYPE_SEQ, oalig.ARRAY_DESCRIPTION_SEQ))
    parser.add_argument("-sm", "--substitution_matrix", choices = oalig.ARRAY_TYPE_MATRIX, default = oalig.MATRIX_BLOSUM62, type = str.upper, help = 'Substitution Matrix type (Only for amino acid sequence) [default: BLOSUM62].')
    # parser.add_argument("-s1", "--sequence1", metavar = "SEQUENCE", required = True, help = "First sequence")
    # parser.add_argument("-s2", "--sequence2", metavar = "SEQUENCE", required = True, help = "Second sequence")
    parser.add_argument("-f", "--fasta", metavar = "FILE", required = True, help = "Fasta file")
    parser.add_argument("-m", "--match", default = 1, type = int, help = "Match value (Only for nucleotide sequence) [default: 1].")
    parser.add_argument("-mi", "--mismatch_penalty", default = 0, type = int, help = "Mismatch penalty value (Only for nucleotide sequence) [default: 0].")
    parser.add_argument("-gap", "--gap_penalty", default = 0, type = int, help = "Gap penalty value [default: 0].")
    parser.add_argument("-o", "--output", metavar = "FOLDER", help = "Output folder")
    parser.add_argument("--version", action = "version", version = "%s %s" % ('%(prog)s', oalig.VERSION))
    args = parser.parse_args()

    oalig.SEQUENCE_TYPE = args.type.lower()
    oalig.MATRIX_TYPE = args.substitution_matrix.upper()
    # oalig.SEQUENCE1 = args.sequence1
    # oalig.SEQUENCE2 = args.sequence2
    oalig.GAP_PENALTY = args.gap_penalty
    oalig.MATCH = args.match
    oalig.MISMATCH_PENALTY = args.mismatch_penalty

    _file = os.path.basename(args.fasta)
    _path = os.path.dirname(_file)
    if _path is None or _path == "":
        _path = os.getcwd().strip()
    oalig.FASTA_FILE = os.path.join(_path, _file)

    if not oalig.check_path(oalig.FASTA_FILE):
        oalig.show_print("The '%s' file doesn't exist!" % oalig.FASTA_FILE, showdate = False)
        exit()

    if args.output is not None:
        output_name = os.path.basename(args.output)
        output_path = os.path.dirname(args.output)
        if output_path is None or output_path == "":
            output_path = os.getcwd().strip()

        oalig.OUTPUT_PATH = os.path.join(output_path, output_name)
        created = oalig.create_directory(oalig.OUTPUT_PATH)
        if not created:
            oalig.show_print("%s: error: Couldn't create folder '%s'" % (os.path.basename(__file__), oalig.OUTPUT_PATH), showdate = False)
            exit()
    else:
        oalig.OUTPUT_PATH = os.getcwd().strip()

class SmithWaterman:

    def __init__(self):
        self.VERSION = 1.0
        self.ROOT = os.path.dirname(os.path.realpath(__file__))

        # Log
        self.LOG_NAME = "log_%s_%s.log" % (os.path.splitext(os.path.basename(__file__))[0], time.strftime('%Y%m%d'))
        self.LOG_FILE = None

        self.TYPE_SEQ_NT = "nt"
        self.TYPE_SEQ_AA = "aa"
        self.DESCRIPTION_SEQ_NT = "Nucleotide sequence"
        self.DESCRIPTION_SEQ_AA = "Amino acid sequence"
        self.ARRAY_TYPE_SEQ = [self.TYPE_SEQ_NT, self.TYPE_SEQ_AA]
        self.ARRAY_DESCRIPTION_SEQ = [self.DESCRIPTION_SEQ_NT, self.DESCRIPTION_SEQ_AA]

        self.MATRIX_BLOSUM45 = "BLOSUM45"
        self.MATRIX_BLOSUM50 = "BLOSUM50"
        self.MATRIX_BLOSUM62 = "BLOSUM62"
        self.MATRIX_BLOSUM80 = "BLOSUM80"
        self.MATRIX_BLOSUM90 = "BLOSUM90"
        self.MATRIX_PAM30 = "PAM30"
        self.MATRIX_PAM70 = "PAM70"
        self.MATRIX_PAM250 = "PAM250"
        self.ARRAY_TYPE_MATRIX = [self.MATRIX_BLOSUM45, self.MATRIX_BLOSUM50, self.MATRIX_BLOSUM62, self.MATRIX_BLOSUM80, self.MATRIX_BLOSUM90, self.MATRIX_PAM30, self.MATRIX_PAM70, self.MATRIX_PAM250]

        self.OUTPUT_PATH = None
        self.HEAD1 = None
        self.HEAD2 = None
        self.SEQUENCE1 = None
        self.SEQUENCE2 = None
        self.FASTA_FILE = None
        self.MATRIX_FILE = 'alignment_matrix.txt'

        self.MATCH = None
        self.MISMATCH_PENALTY = None
        self.GAP_PENALTY = None

        self.ZERO = 0
        self.MAX_SCORE = {'score': 0, 'row': 0, 'col': 0}

        self.MATRIX_TYPE = None
        self.SUBSTITUTION_MATRIX = None

        self.param_diagonal = 'd'
        self.param_up = 'u'
        self.param_left = 'l'
        self.param_zero = 'z'
        self.param_maximum = 'm'
        self.param_direction = 'a'

        self.direct_diagonal = 'diagonal'
        self.direct_up = 'up'
        self.direct_left = 'left'
        self.direct_zero = 'zero'

    def mode_information(self, array1, array2):
        _information = ["%s: %s" % (i, j) for i, j in zip(array1, array2)]
        return " | ".join(_information)

    def show_print(self, message, logs = None, showdate = True, font = None):
        msg_print = message
        msg_write = message

        if font is not None:
            msg_print = "%s%s%s" % (font, msg_print, self.END)

        if showdate is True:
            _time = time.strftime('%Y-%m-%d %H:%M:%S')
            msg_print = "%s %s" % (_time, msg_print)
            msg_write = "%s %s" % (_time, message)

        print(msg_print)
        if logs is not None:
            for log in logs:
                if log is not None:
                    with open(log, 'a', encoding = 'utf-8') as f:
                        f.write("%s\n" % msg_write)
                        f.close()

    def start_time(self):
        return time.time()

    def finish_time(self, start, message = None):
        finish = time.time()
        runtime = time.strftime("%H:%M:%S", time.gmtime(finish - start))
        if message is None:
            return runtime
        else:
            return "%s: %s" % (message, runtime)

    def check_path(self, path):
        if len(path) > 0 and os.path.exists(path):
            return True
        else:
            return False

    def create_directory(self, path):
        output = True
        try:
            if len(path) > 0 and not os.path.exists(path):
                os.makedirs(path)
        except Exception as e:
            output = False
        return output

    def get_num_lines(self, input_file):
        num_lines = sum(1 for line in open(input_file))
        return num_lines

    def init_matrix(self, rows, cols):
        opt = {self.param_diagonal: 0,
               self.param_up: 0,
               self.param_left: 0,
               self.param_zero: 0,
               self.param_maximum: 0,
               self.param_direction: []}
        arr = [opt for i in range(cols + 1)]
        matrix = [arr.copy() for i in range(rows + 1)]

        # First column
        for i in range(rows + 1):
            matrix[i][0] = self.ZERO

        # First row
        for j in range(cols + 1):
            matrix[0][j] = self.ZERO

        return matrix

    def delta(self, si, sj):
        si = si.upper()
        sj = sj.upper()

        value = None
        if self.SEQUENCE_TYPE == self.TYPE_SEQ_NT:
            if si == sj:
                value = self.MATCH
            else:
                value = self.MISMATCH_PENALTY
        elif self.SEQUENCE_TYPE == self.TYPE_SEQ_AA:
            try:
                value = self.SUBSTITUTION_MATRIX[si + sj]
            except Exception as e:
                value = self.SUBSTITUTION_MATRIX[sj + si]

        return value

    def init_substitution_matrix(self):
        self.SUBSTITUTION_MATRIX = {}
        if self.MATRIX_TYPE == self.MATRIX_BLOSUM45:
            self.SUBSTITUTION_MATRIX = {'AA': 5, 'AR': -2, 'AN': -1, 'AD': -2, 'AC': -1, 'AQ': -1, 'AE': -1, 'AG': 0, 'AH': -2, 'AI': -1, 'AL': -1, 'AK': -1,
                                        'AM': -1, 'AF': -2, 'AP': -1, 'AS': 1, 'AT': 0, 'AW': -2, 'AY': -2, 'AV': 0, 'RR': 7, 'RN': 0, 'RD': -1, 'RC': -3,
                                        'RQ': 1, 'RE': 0, 'RG': -2, 'RH': 0, 'RI': -3, 'RL': -2, 'RK': 3, 'RM': -1, 'RF': -2, 'RP': -2, 'RS': -1, 'RT': -1,
                                        'RW': -2, 'RY': -1, 'RV': -2, 'NN': 6, 'ND': 2, 'NC': -2, 'NQ': 0, 'NE': 0, 'NG': 0, 'NH': 1, 'NI': -2, 'NL': -3,
                                        'NK': 0, 'NM': -2, 'NF': -2, 'NP': -2, 'NS': 1, 'NT': 0, 'NW': -4, 'NY': -2, 'NV': -3, 'DD': 7, 'DC': -3, 'DQ': 0,
                                        'DE': 2, 'DG': -1, 'DH': 0, 'DI': -4, 'DL': -3, 'DK': 0, 'DM': -3, 'DF': -4, 'DP': -1, 'DS': 0, 'DT': -1, 'DW': -4,
                                        'DY': -2, 'DV': -3, 'CC': 12, 'CQ': -3, 'CE': -3, 'CG': -3, 'CH': -3, 'CI': -3, 'CL': -2, 'CK': -3, 'CM': -2, 'CF': -2,
                                        'CP': -4, 'CS': -1, 'CT': -1, 'CW': -5, 'CY': -3, 'CV': -1, 'QQ': 6, 'QE': 2, 'QG': -2, 'QH': 1, 'QI': -2, 'QL': -2,
                                        'QK': 1, 'QM': 0, 'QF': -4, 'QP': -1, 'QS': 0, 'QT': -1, 'QW': -2, 'QY': -1, 'QV': -3, 'EE': 6, 'EG': -2, 'EH': 0,
                                        'EI': -3, 'EL': -2, 'EK': 1, 'EM': -2, 'EF': -3, 'EP': 0, 'ES': 0, 'ET': -1, 'EW': -3, 'EY': -2, 'EV': -3, 'GG': 7,
                                        'GH': -2, 'GI': -4, 'GL': -3, 'GK': -2, 'GM': -2, 'GF': -3, 'GP': -2, 'GS': 0, 'GT': -2, 'GW': -2, 'GY': -3, 'GV': -3,
                                        'HH': 10, 'HI': -3, 'HL': -2, 'HK': -1, 'HM': 0, 'HF': -2, 'HP': -2, 'HS': -1, 'HT': -2, 'HW': -3, 'HY': 2, 'HV': -3,
                                        'II': 5, 'IL': 2, 'IK': -3, 'IM': 2, 'IF': 0, 'IP': -2, 'IS': -2, 'IT': -1, 'IW': -2, 'IY': 0, 'IV': 3, 'LL': 5,
                                        'LK': -3, 'LM': 2, 'LF': 1, 'LP': -3, 'LS': -3, 'LT': -1, 'LW': -2, 'LY': 0, 'LV': 1, 'KK': 5, 'KM': -1, 'KF': -3,
                                        'KP': -1, 'KS': -1, 'KT': -1, 'KW': -2, 'KY': -1, 'KV': -2, 'MM': 6, 'MF': 0, 'MP': -2, 'MS': -2, 'MT': -1, 'MW': -2,
                                        'MY': 0, 'MV': 1, 'FF': 8, 'FP': -3, 'FS': -2, 'FT': -1, 'FW': 1, 'FY': 3, 'FV': 0, 'PP': 9, 'PS': -1, 'PT': -1,
                                        'PW': -3, 'PY': -3, 'PV': -3, 'SS': 4, 'ST': 2, 'SW': -4, 'SY': -2, 'SV': -1, 'TT': 5, 'TW': -3, 'TY': -1, 'TV': 0,
                                        'WW': 15, 'WY': 3, 'WV': -3, 'YY': 8, 'YV': -1, 'VV': 5}
        elif self.MATRIX_TYPE == self.MATRIX_BLOSUM50:
            self.SUBSTITUTION_MATRIX = {'AA': 5, 'AR': -2, 'AN': -1, 'AD': -2, 'AC': -1, 'AQ': -1, 'AE': -1, 'AG': 0, 'AH': -2, 'AI': -1, 'AL': -2, 'AK': -1,
                                        'AM': -1, 'AF': -3, 'AP': -1, 'AS': 1, 'AT': 0, 'AW': -3, 'AY': -2, 'AV': 0, 'RR': 7, 'RN': -1, 'RD': -2, 'RC': -4,
                                        'RQ': 1, 'RE': 0, 'RG': -3, 'RH': 0, 'RI': -4, 'RL': -3, 'RK': 3, 'RM': -2, 'RF': -3, 'RP': -3, 'RS': -1, 'RT': -1,
                                        'RW': -3, 'RY': -1, 'RV': -3, 'NN': 7, 'ND': 2, 'NC': -2, 'NQ': 0, 'NE': 0, 'NG': 0, 'NH': 1, 'NI': -3, 'NL': -4,
                                        'NK': 0, 'NM': -2, 'NF': -4, 'NP': -2, 'NS': 1, 'NT': 0, 'NW': -4, 'NY': -2, 'NV': -3, 'DD': 8, 'DC': -4, 'DQ': 0,
                                        'DE': 2, 'DG': -1, 'DH': -1, 'DI': -4, 'DL': -4, 'DK': -1, 'DM': -4, 'DF': -5, 'DP': -1, 'DS': 0, 'DT': -1, 'DW': -5,
                                        'DY': -3, 'DV': -4, 'CC': 13, 'CQ': -3, 'CE': -3, 'CG': -3, 'CH': -3, 'CI': -2, 'CL': -2, 'CK': -3, 'CM': -2, 'CF': -2,
                                        'CP': -4, 'CS': -1, 'CT': -1, 'CW': -5, 'CY': -3, 'CV': -1, 'QQ': 7, 'QE': 2, 'QG': -2, 'QH': 1, 'QI': -3, 'QL': -2,
                                        'QK': 2, 'QM': 0, 'QF': -4, 'QP': -1, 'QS': 0, 'QT': -1, 'QW': -1, 'QY': -1, 'QV': -3, 'EE': 6, 'EG': -3, 'EH': 0,
                                        'EI': -4, 'EL': -3, 'EK': 1, 'EM': -2, 'EF': -3, 'EP': -1, 'ES': -1, 'ET': -1, 'EW': -3, 'EY': -2, 'EV': -3, 'GG': 8,
                                        'GH': -2, 'GI': -4, 'GL': -4, 'GK': -2, 'GM': -3, 'GF': -4, 'GP': -2, 'GS': 0, 'GT': -2, 'GW': -3, 'GY': -3, 'GV': -4,
                                        'HH': 10, 'HI': -4, 'HL': -3, 'HK': 0, 'HM': -1, 'HF': -1, 'HP': -2, 'HS': -1, 'HT': -2, 'HW': -3, 'HY': 2, 'HV': -4,
                                        'II': 5, 'IL': 2, 'IK': -3, 'IM': 2, 'IF': 0, 'IP': -3, 'IS': -3, 'IT': -1, 'IW': -3, 'IY': -1, 'IV': 4, 'LL': 5,
                                        'LK': -3, 'LM': 3, 'LF': 1, 'LP': -4, 'LS': -3, 'LT': -1, 'LW': -2, 'LY': -1, 'LV': 1, 'KK': 6, 'KM': -2, 'KF': -4,
                                        'KP': -1, 'KS': 0, 'KT': -1, 'KW': -3, 'KY': -2, 'KV': -3, 'MM': 7, 'MF': 0, 'MP': -3, 'MS': -2, 'MT': -1, 'MW': -1,
                                        'MY': 0, 'MV': 1, 'FF': 8, 'FP': -4, 'FS': -3, 'FT': -2, 'FW': 1, 'FY': 4, 'FV': -1, 'PP': 10, 'PS': -1, 'PT': -1,
                                        'PW': -4, 'PY': -3, 'PV': -3, 'SS': 5, 'ST': 2, 'SW': -4, 'SY': -2, 'SV': -2, 'TT': 5, 'TW': -3, 'TY': -2, 'TV': 0,
                                        'WW': 15, 'WY': 2, 'WV': -3, 'YY': 8, 'YV': -1, 'VV': 5}
        elif self.MATRIX_TYPE == self.MATRIX_BLOSUM62:
            self.SUBSTITUTION_MATRIX = {'AA': 4, 'AR': -1, 'AN': -2, 'AD': -2, 'AC': 0, 'AQ': -1, 'AE': -1, 'AG': 0, 'AH': -2, 'AI': -1, 'AL': -1, 'AK': -1,
                                        'AM': -1, 'AF': -2, 'AP': -1, 'AS': 1, 'AT': 0, 'AW': -3, 'AY': -2, 'AV': 0, 'RR': 5, 'RN': 0, 'RD': -2, 'RC': -3,
                                        'RQ': 1, 'RE': 0, 'RG': -2, 'RH': 0, 'RI': -3, 'RL': -2, 'RK': 2, 'RM': -1, 'RF': -3, 'RP': -2, 'RS': -1, 'RT': -1,
                                        'RW': -3, 'RY': -2, 'RV': -3, 'NN': 6, 'ND': 1, 'NC': -3, 'NQ': 0, 'NE': 0, 'NG': 0, 'NH': 1, 'NI': -3, 'NL': -3,
                                        'NK': 0, 'NM': -2, 'NF': -3, 'NP': -2, 'NS': 1, 'NT': 0, 'NW': -4, 'NY': -2, 'NV': -3, 'DD': 6, 'DC': -3, 'DQ': 0,
                                        'DE': 2, 'DG': -1, 'DH': -1, 'DI': -3, 'DL': -4, 'DK': -1, 'DM': -3, 'DF': -3, 'DP': -1, 'DS': 0, 'DT': -1, 'DW': -4,
                                        'DY': -3, 'DV': -3, 'CC': 9, 'CQ': -3, 'CE': -4, 'CG': -3, 'CH': -3, 'CI': -1, 'CL': -1, 'CK': -3, 'CM': -1, 'CF': -2,
                                        'CP': -3, 'CS': -1, 'CT': -1, 'CW': -2, 'CY': -2, 'CV': -1, 'QQ': 5, 'QE': 2, 'QG': -2, 'QH': 0, 'QI': -3, 'QL': -2,
                                        'QK': 1, 'QM': 0, 'QF': -3, 'QP': -1, 'QS': 0, 'QT': -1, 'QW': -2, 'QY': -1, 'QV': -2, 'EE': 5, 'EG': -2, 'EH': 0,
                                        'EI': -3, 'EL': -3, 'EK': 1, 'EM': -2, 'EF': -3, 'EP': -1, 'ES': 0, 'ET': -1, 'EW': -3, 'EY': -2, 'EV': -2, 'GG': 6,
                                        'GH': -2, 'GI': -4, 'GL': -4, 'GK': -2, 'GM': -3, 'GF': -3, 'GP': -2, 'GS': 0, 'GT': -2, 'GW': -2, 'GY': -3, 'GV': -3,
                                        'HH': 8, 'HI': -3, 'HL': -3, 'HK': -1, 'HM': -2, 'HF': -1, 'HP': -2, 'HS': -1, 'HT': -2, 'HW': -2, 'HY': 2, 'HV': -3,
                                        'II': 4, 'IL': 2, 'IK': -3, 'IM': 1, 'IF': 0, 'IP': -3, 'IS': -2, 'IT': -1, 'IW': -3, 'IY': -1, 'IV': 3, 'LL': 4,
                                        'LK': -2, 'LM': 2, 'LF': 0, 'LP': -3, 'LS': -2, 'LT': -1, 'LW': -2, 'LY': -1, 'LV': 1, 'KK': 5, 'KM': -1, 'KF': -3,
                                        'KP': -1, 'KS': 0, 'KT': -1, 'KW': -3, 'KY': -2, 'KV': -2, 'MM': 5, 'MF': 0, 'MP': -2, 'MS': -1, 'MT': -1, 'MW': -1,
                                        'MY': -1, 'MV': 1, 'FF': 6, 'FP': -4, 'FS': -2, 'FT': -2, 'FW': 1, 'FY': 3, 'FV': -1, 'PP': 7, 'PS': -1, 'PT': -1,
                                        'PW': -4, 'PY': -3, 'PV': -2, 'SS': 4, 'ST': 1, 'SW': -3, 'SY': -2, 'SV': -2, 'TT': 5, 'TW': -2, 'TY': -2, 'TV': 0,
                                        'WW': 11, 'WY': 2, 'WV': -3, 'YY': 7, 'YV': -1, 'VV': 4}
        elif self.MATRIX_TYPE == self.MATRIX_BLOSUM80:
            self.SUBSTITUTION_MATRIX = {'AA': 5, 'AR': -2, 'AN': -2, 'AD': -2, 'AC': -1, 'AQ': -1, 'AE': -1, 'AG': 0, 'AH': -2, 'AI': -2, 'AL': -2, 'AK': -1,
                                        'AM': -1, 'AF': -3, 'AP': -1, 'AS': 1, 'AT': 0, 'AW': -3, 'AY': -2, 'AV': 0, 'RR': 6, 'RN': -1, 'RD': -2, 'RC': -4,
                                        'RQ': 1, 'RE': -1, 'RG': -3, 'RH': 0, 'RI': -3, 'RL': -3, 'RK': 2, 'RM': -2, 'RF': -4, 'RP': -2, 'RS': -1, 'RT': -1,
                                        'RW': -4, 'RY': -3, 'RV': -3, 'NN': 6, 'ND': 1, 'NC': -3, 'NQ': 0, 'NE': -1, 'NG': -1, 'NH': 0, 'NI': -4, 'NL': -4,
                                        'NK': 0, 'NM': -3, 'NF': -4, 'NP': -3, 'NS': 0, 'NT': 0, 'NW': -4, 'NY': -3, 'NV': -4, 'DD': 6, 'DC': -4, 'DQ': -1,
                                        'DE': 1, 'DG': -2, 'DH': -2, 'DI': -4, 'DL': -5, 'DK': -1, 'DM': -4, 'DF': -4, 'DP': -2, 'DS': -1, 'DT': -1, 'DW': -6,
                                        'DY': -4, 'DV': -4, 'CC': 9, 'CQ': -4, 'CE': -5, 'CG': -4, 'CH': -4, 'CI': -2, 'CL': -2, 'CK': -4, 'CM': -2, 'CF': -3,
                                        'CP': -4, 'CS': -2, 'CT': -1, 'CW': -3, 'CY': -3, 'CV': -1, 'QQ': 6, 'QE': 2, 'QG': -2, 'QH': 1, 'QI': -3, 'QL': -3,
                                        'QK': 1, 'QM': 0, 'QF': -4, 'QP': -2, 'QS': 0, 'QT': -1, 'QW': -3, 'QY': -2, 'QV': -3, 'EE': 6, 'EG': -3, 'EH': 0,
                                        'EI': -4, 'EL': -4, 'EK': 1, 'EM': -2, 'EF': -4, 'EP': -2, 'ES': 0, 'ET': -1, 'EW': -4, 'EY': -3, 'EV': -3, 'GG': 6,
                                        'GH': -3, 'GI': -5, 'GL': -4, 'GK': -2, 'GM': -4, 'GF': -4, 'GP': -3, 'GS': -1, 'GT': -2, 'GW': -4, 'GY': -4, 'GV': -4,
                                        'HH': 8, 'HI': -4, 'HL': -3, 'HK': -1, 'HM': -2, 'HF': -2, 'HP': -3, 'HS': -1, 'HT': -2, 'HW': -3, 'HY': 2, 'HV': -4,
                                        'II': 5, 'IL': 1, 'IK': -3, 'IM': 1, 'IF': -1, 'IP': -4, 'IS': -3, 'IT': -1, 'IW': -3, 'IY': -2, 'IV': 3, 'LL': 4,
                                        'LK': -3, 'LM': 2, 'LF': 0, 'LP': -3, 'LS': -3, 'LT': -2, 'LW': -2, 'LY': -2, 'LV': 1, 'KK': 5, 'KM': -2, 'KF': -4,
                                        'KP': -1, 'KS': -1, 'KT': -1, 'KW': -4, 'KY': -3, 'KV': -3, 'MM': 6, 'MF': 0, 'MP': -3, 'MS': -2, 'MT': -1, 'MW': -2,
                                        'MY': -2, 'MV': 1, 'FF': 6, 'FP': -4, 'FS': -3, 'FT': -2, 'FW': 0, 'FY': 3, 'FV': -1, 'PP': 8, 'PS': -1, 'PT': -2,
                                        'PW': -5, 'PY': -4, 'PV': -3, 'SS': 5, 'ST': 1, 'SW': -4, 'SY': -2, 'SV': -2, 'TT': 5, 'TW': -4, 'TY': -2, 'TV': 0,
                                        'WW': 11, 'WY': 2, 'WV': -3, 'YY': 7, 'YV': -2, 'VV': 4}
        elif self.MATRIX_TYPE == self.MATRIX_BLOSUM90:
            self.SUBSTITUTION_MATRIX = {'AA': 5, 'AR': -2, 'AN': -2, 'AD': -3, 'AC': -1, 'AQ': -1, 'AE': -1, 'AG': 0, 'AH': -2, 'AI': -2, 'AL': -2, 'AK': -1,
                                        'AM': -2, 'AF': -3, 'AP': -1, 'AS': 1, 'AT': 0, 'AW': -4, 'AY': -3, 'AV': -1, 'RR': 6, 'RN': -1, 'RD': -3, 'RC': -5,
                                        'RQ': 1, 'RE': -1, 'RG': -3, 'RH': 0, 'RI': -4, 'RL': -3, 'RK': 2, 'RM': -2, 'RF': -4, 'RP': -3, 'RS': -1, 'RT': -2,
                                        'RW': -4, 'RY': -3, 'RV': -3, 'NN': 7, 'ND': 1, 'NC': -4, 'NQ': 0, 'NE': -1, 'NG': -1, 'NH': 0, 'NI': -4, 'NL': -4,
                                        'NK': 0, 'NM': -3, 'NF': -4, 'NP': -3, 'NS': 0, 'NT': 0, 'NW': -5, 'NY': -3, 'NV': -4, 'DD': 7, 'DC': -5, 'DQ': -1,
                                        'DE': 1, 'DG': -2, 'DH': -2, 'DI': -5, 'DL': -5, 'DK': -1, 'DM': -4, 'DF': -5, 'DP': -3, 'DS': -1, 'DT': -2, 'DW': -6,
                                        'DY': -4, 'DV': -5, 'CC': 9, 'CQ': -4, 'CE': -6, 'CG': -4, 'CH': -5, 'CI': -2, 'CL': -2, 'CK': -4, 'CM': -2, 'CF': -3,
                                        'CP': -4, 'CS': -2, 'CT': -2, 'CW': -4, 'CY': -4, 'CV': -2, 'QQ': 7, 'QE': 2, 'QG': -3, 'QH': 1, 'QI': -4, 'QL': -3,
                                        'QK': 1, 'QM': 0, 'QF': -4, 'QP': -2, 'QS': -1, 'QT': -1, 'QW': -3, 'QY': -3, 'QV': -3, 'EE': 6, 'EG': -3, 'EH': -1,
                                        'EI': -4, 'EL': -4, 'EK': 0, 'EM': -3, 'EF': -5, 'EP': -2, 'ES': -1, 'ET': -1, 'EW': -5, 'EY': -4, 'EV': -3, 'GG': 6,
                                        'GH': -3, 'GI': -5, 'GL': -5, 'GK': -2, 'GM': -4, 'GF': -5, 'GP': -3, 'GS': -1, 'GT': -3, 'GW': -4, 'GY': -5, 'GV': -5,
                                        'HH': 8, 'HI': -4, 'HL': -4, 'HK': -1, 'HM': -3, 'HF': -2, 'HP': -3, 'HS': -2, 'HT': -2, 'HW': -3, 'HY': 1, 'HV': -4,
                                        'II': 5, 'IL': 1, 'IK': -4, 'IM': 1, 'IF': -1, 'IP': -4, 'IS': -3, 'IT': -1, 'IW': -4, 'IY': -2, 'IV': 3, 'LL': 5,
                                        'LK': -3, 'LM': 2, 'LF': 0, 'LP': -4, 'LS': -3, 'LT': -2, 'LW': -3, 'LY': -2, 'LV': 0, 'KK': 6, 'KM': -2, 'KF': -4,
                                        'KP': -2, 'KS': -1, 'KT': -1, 'KW': -5, 'KY': -3, 'KV': -3, 'MM': 7, 'MF': -1, 'MP': -3, 'MS': -2, 'MT': -1, 'MW': -2,
                                        'MY': -2, 'MV': 0, 'FF': 7, 'FP': -4, 'FS': -3, 'FT': -3, 'FW': 0, 'FY': 3, 'FV': -2, 'PP': 8, 'PS': -2, 'PT': -2,
                                        'PW': -5, 'PY': -4, 'PV': -3, 'SS': 5, 'ST': 1, 'SW': -4, 'SY': -3, 'SV': -2, 'TT': 6, 'TW': -4, 'TY': -2, 'TV': -1,
                                        'WW': 11, 'WY': 2, 'WV': -3, 'YY': 8, 'YV': -3, 'VV': 5}
        elif self.MATRIX_TYPE == self.MATRIX_PAM30:
            self.SUBSTITUTION_MATRIX = {'AA': 6, 'AR': -7, 'AN': -4, 'AD': -3, 'AC': -6, 'AQ': -4, 'AE': -2, 'AG': -2, 'AH': -7, 'AI': -5, 'AL': -6, 'AK': -7,
                                        'AM': -5, 'AF': -8, 'AP': -2, 'AS': 0, 'AT': -1, 'AW': -13, 'AY': -8, 'AV': -2, 'RR': 8, 'RN': -6, 'RD': -10, 'RC': -8,
                                        'RQ': -2, 'RE': -9, 'RG': -9, 'RH': -2, 'RI': -5, 'RL': -8, 'RK': 0, 'RM': -4, 'RF': -9, 'RP': -4, 'RS': -3, 'RT': -6,
                                        'RW': -2, 'RY': -10, 'RV': -8, 'NN': 8, 'ND': 2, 'NC': -11, 'NQ': -3, 'NE': -2, 'NG': -3, 'NH': 0, 'NI': -5, 'NL': -7,
                                        'NK': -1, 'NM': -9, 'NF': -9, 'NP': -6, 'NS': 0, 'NT': -2, 'NW': -8, 'NY': -4, 'NV': -8, 'DD': 8, 'DC': -14, 'DQ': -2,
                                        'DE': 2, 'DG': -3, 'DH': -4, 'DI': -7, 'DL': -12, 'DK': -4, 'DM': -11, 'DF': -15, 'DP': -8, 'DS': -4, 'DT': -5, 'DW': -15,
                                        'DY': -11, 'DV': -8, 'CC': 10, 'CQ': -14, 'CE': -14, 'CG': -9, 'CH': -7, 'CI': -6, 'CL': -15, 'CK': -14, 'CM': -13,
                                        'CF': -13, 'CP': -8, 'CS': -3, 'CT': -8, 'CW': -15, 'CY': -4, 'CV': -6, 'QQ': 8, 'QE': 1, 'QG': -7, 'QH': 1, 'QI': -8,
                                        'QL': -5, 'QK': -3, 'QM': -4, 'QF': -13, 'QP': -3, 'QS': -5, 'QT': -5, 'QW': -13, 'QY': -12, 'QV': -7, 'EE': 8, 'EG': -4,
                                        'EH': -5, 'EI': -5, 'EL': -9, 'EK': -4, 'EM': -7, 'EF': -14, 'EP': -5, 'ES': -4, 'ET': -6, 'EW': -17, 'EY': -8, 'EV': -6,
                                        'GG': 6, 'GH': -9, 'GI': -11, 'GL': -10, 'GK': -7, 'GM': -8, 'GF': -9, 'GP': -6, 'GS': -2, 'GT': -6, 'GW': -15, 'GY': -14,
                                        'GV': -5, 'HH': 9, 'HI': -9, 'HL': -6, 'HK': -6, 'HM': -10, 'HF': -6, 'HP': -4, 'HS': -6, 'HT': -7, 'HW': -7, 'HY': -3,
                                        'HV': -6, 'II': 8, 'IL': -1, 'IK': -6, 'IM': -1, 'IF': -2, 'IP': -8, 'IS': -7, 'IT': -2, 'IW': -14, 'IY': -6, 'IV': 2,
                                        'LL': 7, 'LK': -8, 'LM': 1, 'LF': -3, 'LP': -7, 'LS': -8, 'LT': -7, 'LW': -6, 'LY': -7, 'LV': -2, 'KK': 7, 'KM': -2,
                                        'KF': -14, 'KP': -6, 'KS': -4, 'KT': -3, 'KW': -12, 'KY': -9, 'KV': -9, 'MM': 11, 'MF': -4, 'MP': -8, 'MS': -5, 'MT': -4,
                                        'MW': -13, 'MY': -11, 'MV': -1, 'FF': 9, 'FP': -10, 'FS': -6, 'FT': -9, 'FW': -4, 'FY': 2, 'FV': -8, 'PP': 8, 'PS': -2,
                                        'PT': -4, 'PW': -14, 'PY': -13, 'PV': -6, 'SS': 6, 'ST': 0, 'SW': -5, 'SY': -7, 'SV': -6, 'TT': 7, 'TW': -13, 'TY': -6,
                                        'TV': -3, 'WW': 13, 'WY': -5, 'WV': -15, 'YY': 10, 'YV': -7, 'VV': 7}
        elif self.MATRIX_TYPE == self.MATRIX_PAM70:
            self.SUBSTITUTION_MATRIX = {'AA': 5, 'AR': -4, 'AN': -2, 'AD': -1, 'AC': -4, 'AQ': -2, 'AE': -1, 'AG': 0, 'AH': -4, 'AI': -2, 'AL': -4, 'AK': -4,
                                        'AM': -3, 'AF': -6, 'AP': 0, 'AS': 1, 'AT': 1, 'AW': -9, 'AY': -5, 'AV': -1, 'RR': 8, 'RN': -3, 'RD': -6, 'RC': -5,
                                        'RQ': 0, 'RE': -5, 'RG': -6, 'RH': 0, 'RI': -3, 'RL': -6, 'RK': 2, 'RM': -2, 'RF': -7, 'RP': -2, 'RS': -1, 'RT': -4,
                                        'RW': 0, 'RY': -7, 'RV': -5, 'NN': 6, 'ND': 3, 'NC': -7, 'NQ': -1, 'NE': 0, 'NG': -1, 'NH': 1, 'NI': -3, 'NL': -5,
                                        'NK': 0, 'NM': -5, 'NF': -6, 'NP': -3, 'NS': 1, 'NT': 0, 'NW': -6, 'NY': -3, 'NV': -5, 'DD': 6, 'DC': -9, 'DQ': 0,
                                        'DE': 3, 'DG': -1, 'DH': -1, 'DI': -5, 'DL': -8, 'DK': -2, 'DM': -7, 'DF': -10, 'DP': -4, 'DS': -1, 'DT': -2, 'DW': -10,
                                        'DY': -7, 'DV': -5, 'CC': 9, 'CQ': -9, 'CE': -9, 'CG': -6, 'CH': -5, 'CI': -4, 'CL': -10, 'CK': -9, 'CM': -9, 'CF': -8,
                                        'CP': -5, 'CS': -1, 'CT': -5, 'CW': -11, 'CY': -2, 'CV': -4, 'QQ': 7, 'QE': 2, 'QG': -4, 'QH': 2, 'QI': -5, 'QL': -3,
                                        'QK': -1, 'QM': -2, 'QF': -9, 'QP': -1, 'QS': -3, 'QT': -3, 'QW': -8, 'QY': -8, 'QV': -4, 'EE': 6, 'EG': -2, 'EH': -2,
                                        'EI': -4, 'EL': -6, 'EK': -2, 'EM': -4, 'EF': -9, 'EP': -3, 'ES': -2, 'ET': -3, 'EW': -11, 'EY': -6, 'EV': -4, 'GG': 6,
                                        'GH': -6, 'GI': -6, 'GL': -7, 'GK': -5, 'GM': -6, 'GF': -7, 'GP': -3, 'GS': 0, 'GT': -3, 'GW': -10, 'GY': -9, 'GV': -3,
                                        'HH': 8, 'HI': -6, 'HL': -4, 'HK': -3, 'HM': -6, 'HF': -4, 'HP': -2, 'HS': -3, 'HT': -4, 'HW': -5, 'HY': -1, 'HV': -4,
                                        'II': 7, 'IL': 1, 'IK': -4, 'IM': 1, 'IF': 0, 'IP': -5, 'IS': -4, 'IT': -1, 'IW': -9, 'IY': -4, 'IV': 3, 'LL': 6, 'LK': -5,
                                        'LM': 2, 'LF': -1, 'LP': -5, 'LS': -6, 'LT': -4, 'LW': -4, 'LY': -4, 'LV': 0, 'KK': 6, 'KM': 0, 'KF': -9, 'KP': -4,
                                        'KS': -2, 'KT': -1, 'KW': -7, 'KY': -7, 'KV': -6, 'MM': 10, 'MF': -2, 'MP': -5, 'MS': -3, 'MT': -2, 'MW': -8, 'MY': -7,
                                        'MV': 0, 'FF': 8, 'FP': -7, 'FS': -4, 'FT': -6, 'FW': -2, 'FY': 4, 'FV': -5, 'PP': 7, 'PS': 0, 'PT': -2, 'PW': -9,
                                        'PY': -9, 'PV': -3, 'SS': 5, 'ST': 2, 'SW': -3, 'SY': -5, 'SV': -3, 'TT': 6, 'TW': -8, 'TY': -4, 'TV': -1, 'WW': 13,
                                        'WY': -3, 'WV': -10, 'YY': 9, 'YV': -5, 'VV': 6}
        elif self.MATRIX_TYPE == self.MATRIX_PAM250:
            self.SUBSTITUTION_MATRIX = {'AA': 2, 'AR': -2, 'AN': 0, 'AD': 0, 'AC': -2, 'AQ': 0, 'AE': 0, 'AG': 1, 'AH': -1, 'AI': -1, 'AL': -2, 'AK': -1,
                                        'AM': -1, 'AF': -3, 'AP': 1, 'AS': 1, 'AT': 1, 'AW': -6, 'AY': -3, 'AV': 0, 'RR': 6, 'RN': 0, 'RD': -1, 'RC': -4,
                                        'RQ': 1, 'RE': -1, 'RG': -3, 'RH': 2, 'RI': -2, 'RL': -3, 'RK': 3, 'RM': 0, 'RF': -4, 'RP': 0, 'RS': 0, 'RT': -1,
                                        'RW': 2, 'RY': -4, 'RV': -2, 'NN': 2, 'ND': 2, 'NC': -4, 'NQ': 1, 'NE': 1, 'NG': 0, 'NH': 2, 'NI': -2, 'NL': -3,
                                        'NK': 1, 'NM': -2, 'NF': -3, 'NP': 0, 'NS': 1, 'NT': 0, 'NW': -4, 'NY': -2, 'NV': -2, 'DD': 4, 'DC': -5, 'DQ': 2,
                                        'DE': 3, 'DG': 1, 'DH': 1, 'DI': -2, 'DL': -4, 'DK': 0, 'DM': -3, 'DF': -6, 'DP': -1, 'DS': 0, 'DT': 0, 'DW': -7,
                                        'DY': -4, 'DV': -2, 'CC': 12, 'CQ': -5, 'CE': -5, 'CG': -3, 'CH': -3, 'CI': -2, 'CL': -6, 'CK': -5, 'CM': -5, 'CF': -4,
                                        'CP': -3, 'CS': 0, 'CT': -2, 'CW': -8, 'CY': 0, 'CV': -2, 'QQ': 4, 'QE': 2, 'QG': -1, 'QH': 3, 'QI': -2, 'QL': -2,
                                        'QK': 1, 'QM': -1, 'QF': -5, 'QP': 0, 'QS': -1, 'QT': -1, 'QW': -5, 'QY': -4, 'QV': -2, 'EE': 4, 'EG': 0, 'EH': 1,
                                        'EI': -2, 'EL': -3, 'EK': 0, 'EM': -2, 'EF': -5, 'EP': -1, 'ES': 0, 'ET': 0, 'EW': -7, 'EY': -4, 'EV': -2, 'GG': 5,
                                        'GH': -2, 'GI': -3, 'GL': -4, 'GK': -2, 'GM': -3, 'GF': -5, 'GP': 0, 'GS': 1, 'GT': 0, 'GW': -7, 'GY': -5, 'GV': -1,
                                        'HH': 6, 'HI': -2, 'HL': -2, 'HK': 0, 'HM': -2, 'HF': -2, 'HP': 0, 'HS': -1, 'HT': -1, 'HW': -3, 'HY': 0, 'HV': -2,
                                        'II': 5, 'IL': 2, 'IK': -2, 'IM': 2, 'IF': 1, 'IP': -2, 'IS': -1, 'IT': 0, 'IW': -5, 'IY': -1, 'IV': 4, 'LL': 6,
                                        'LK': -3, 'LM': 4, 'LF': 2, 'LP': -3, 'LS': -3, 'LT': -2, 'LW': -2, 'LY': -1, 'LV': 2, 'KK': 5, 'KM': 0, 'KF': -5,
                                        'KP': -1, 'KS': 0, 'KT': 0, 'KW': -3, 'KY': -4, 'KV': -2, 'MM': 6, 'MF': 0, 'MP': -2, 'MS': -2, 'MT': -1, 'MW': -4,
                                        'MY': -2, 'MV': 2, 'FF': 9, 'FP': -5, 'FS': -3, 'FT': -3, 'FW': 0, 'FY': 7, 'FV': -1, 'PP': 6, 'PS': 1, 'PT': 0,
                                        'PW': -6, 'PY': -5, 'PV': -1, 'SS': 2, 'ST': 1, 'SW': -2, 'SY': -3, 'SV': -1, 'TT': 3, 'TW': -5, 'TY': -3, 'TV': 0,
                                        'WW': 17, 'WY': 0, 'WV': -6, 'YY': 10, 'YV': -2, 'VV': 4}

    def smith_waterman(self, seq1, seq2):
        arr_seq1 = [i for i in seq1]
        arr_seq2 = [i for i in seq2]
        n = len(arr_seq1)
        m = len(arr_seq2)
        matrix = self.init_matrix(m, n)

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                options_d = matrix[i - 1][j - 1]
                options_u = matrix[i - 1][j]
                options_l = matrix[i][j - 1]
                options = matrix[i][j].copy()

                value_d = options_d[self.param_maximum] if isinstance(options_d, dict) else options_d
                value_u = options_u[self.param_maximum] if isinstance(options_u, dict) else options_u
                value_l = options_l[self.param_maximum] if isinstance(options_l, dict) else options_l

                match = value_d + self.delta(arr_seq1[j - 1], arr_seq2[i - 1])
                delete = value_l + self.GAP_PENALTY
                insert = value_u + self.GAP_PENALTY
                maximum = max(match, delete, insert, self.ZERO)

                direction = options[self.param_direction].copy()
                if match < 0 and delete < 0 and insert < 0:
                    direction.append(self.direct_zero)
                if maximum == match:
                    direction.append(self.direct_diagonal)
                if maximum == delete:
                    direction.append(self.direct_left)
                if maximum == insert:
                    direction.append(self.direct_up)
                options[self.param_direction] = direction.copy()

                options[self.param_diagonal] = match
                options[self.param_left] = delete
                options[self.param_up] = insert
                options[self.param_zero] = self.ZERO
                options[self.param_maximum] = maximum
                matrix[i][j] = options

                if maximum > self.MAX_SCORE['score']:
                    self.MAX_SCORE['score'] = maximum
                    self.MAX_SCORE['row'] = i
                    self.MAX_SCORE['col'] = j

        # Traceback
        align1 = ''
        align2 = ''
        i = self.MAX_SCORE['row']
        j = self.MAX_SCORE['col']
        while i > 0 and j > 0:
            options_current = matrix[i][j]
            direction = options_current[self.param_direction]

            if options_current[self.param_maximum] == 0:
                break

            if self.direct_diagonal in direction:
                # Optimal alignment
                # Diagonal
                align1 += arr_seq1[j - 1]
                align2 += arr_seq2[i - 1]
                i -= 1
                j -= 1
            elif self.direct_left in direction:
                # Horizontal
                align1 += arr_seq1[j - 1]
                align2 += '-'
                j -= 1
            elif self.direct_up in direction:
                # Vertical
                align1 += '-'
                align2 += arr_seq2[i - 1]
                i -= 1

        align1 = align1[::-1]
        align2 = align2[::-1]

        return align1, align2, self.MAX_SCORE['score'], matrix

    def matrix_format(self, matrix, arr_seq1, arr_seq2):
        matrix_sum = []
        matrix_row = []
        for row in matrix:
            for cell in row:
                if isinstance(cell, dict):
                    direction = cell[self.param_direction]
                    maximum = cell[self.param_maximum]

                    if self.direct_zero in direction:
                        matrix_row.append('%s|%s' % (maximum, 'z'))
                    elif self.direct_diagonal in direction:
                        matrix_row.append('%s|%s' % (maximum, 'd'))
                    elif self.direct_left in direction:
                        matrix_row.append('%s|%s' % (maximum, 'l'))
                    elif self.direct_up in direction:
                        matrix_row.append('%s|%s' % (maximum, 'u'))

            if matrix_row:
                matrix_sum.append(matrix_row)
                matrix_row = []

        nchars = 0
        matrix_full = []
        for row in matrix_sum:
            for cell in row:
                nchars = len(cell) if len(cell) > nchars else nchars
        for index, row in enumerate(matrix_sum, start = 1):
            if index == 1:
                head_row = ['-', '-']
                first_row = ['-', '0']
                for idx, cell in enumerate(row, start = 1):
                    head_row.append(arr_seq1[idx - 1].rjust(nchars, ' '))
                    first_row.append('0'.rjust(nchars, ' '))
                matrix_full.append(head_row)
                matrix_full.append(first_row)
            rows = [arr_seq2[index - 1], '0']
            for cell in row:
                rows.append(cell.rjust(nchars, ' '))
            matrix_full.append(rows)

        return matrix_full

    def get_sequences(self):
        collect_seq = OrderedDict()

        if self.get_num_lines(self.FASTA_FILE) == 0:
            self.show_print("The '%s' file is empty" % (os.path.basename(self.FASTA_FILE)), [self.LOG_FILE])
            exit()
        else:
            head = ''
            sequence = ''
            num_seq = 1
            with open(self.FASTA_FILE, 'r') as fr:
                for line in fr:
                    line = line.strip()

                    if line.startswith('>'):
                        if sequence:
                            collect_seq.update({head: sequence})
                            head = ''
                            sequence = ''
                            num_seq += 1

                        if num_seq > 2: # Num sequences
                            break

                        head = line.split(' ')[0]
                    else:
                        sequence += line

                if sequence:
                    collect_seq.update({head: sequence})
            fr.close()

            self.HEAD1, self.HEAD2 = list(collect_seq.keys())
            self.HEAD1 = self.HEAD1.replace('>', '')
            self.HEAD2 = self.HEAD2.replace('>', '')
            self.SEQUENCE1, self.SEQUENCE2 = list(collect_seq.values())

        return collect_seq

    def get_alignment_characters(self, alignment1, alignment2):
        merged_list = [(alignment1[i], alignment2[i]) for i in range(0, len(alignment1))]

        align1 = ''
        align2 = ''
        characters = ''
        for tupla in merged_list:
            si = tupla[0]
            sj = tupla[1]

            char = '•'
            if si.upper() == sj.upper():
                char = '|'
            elif si == '-' or sj == '-':
                char = ' '

            align1 += si
            align2 += sj
            characters += char

        size_left = 0
        if len(self.HEAD1) > size_left:
            size_left = len(self.HEAD1)
        if len(self.HEAD2) > size_left:
            size_left = len(self.HEAD2)
        size_left += 4
        size_subseq = 60

        alignment = []
        if len(align1) <= size_subseq:
            _alignment = []
            _alignment.append('%s%s' % (self.HEAD1.ljust(size_left), align1))
            _alignment.append('%s%s' % (''.ljust(size_left), characters))
            _alignment.append('%s%s' % (self.HEAD2.ljust(size_left), align2))
            alignment.append(_alignment)
        else:
            total = len(align1)
            ini = 0
            end = size_subseq
            while total > 0:
                _alignment = []
                _alignment.append('%s%s' % (self.HEAD1.ljust(size_left), align1[ini:end]))
                _alignment.append('%s%s' % (''.ljust(size_left), characters[ini:end]))
                _alignment.append('%s%s' % (self.HEAD2.ljust(size_left), align2[ini:end]))
                alignment.append(_alignment)
                ini = end
                end += size_subseq
                total = total - size_subseq

        return alignment

    def save_matrix(self, matrix):
        arr_seq1 = [i for i in self.SEQUENCE1]
        arr_seq2 = [i for i in self.SEQUENCE2]

        with open(self.MATRIX_FILE, 'w') as fw:
            for row in self.matrix_format(matrix, arr_seq1, arr_seq2):
                fw.write('%s\n' % ' '.join(row))
        fw.close()

def main(args):
    try:
        start = oalig.start_time()
        menu(args)
        oalig.create_directory(oalig.OUTPUT_PATH)
        oalig.LOG_FILE = os.path.join(oalig.OUTPUT_PATH, oalig.LOG_NAME)
        oalig.MATRIX_FILE = os.path.join(oalig.OUTPUT_PATH, oalig.MATRIX_FILE)

        oalig.show_print("################################################################", [oalig.LOG_FILE])
        oalig.show_print("################### Smith–Waterman Algorithm ###################", [oalig.LOG_FILE])
        oalig.show_print("################################################################", [oalig.LOG_FILE])

        oalig.show_print("Input:", [oalig.LOG_FILE])
        oalig.show_print("  Fasta file: %s" % oalig.FASTA_FILE, [oalig.LOG_FILE])
        # oalig.show_print("  Sequence 1: %s" % oalig.SEQUENCE1, [oalig.LOG_FILE])
        # oalig.show_print("  Sequence 2: %s" % oalig.SEQUENCE2, [oalig.LOG_FILE])
        oalig.show_print("", [oalig.LOG_FILE])

        oalig.show_print("Parameters:", [oalig.LOG_FILE])
        if oalig.SEQUENCE_TYPE == oalig.TYPE_SEQ_NT:
            oalig.show_print("  Match: %s" % oalig.MATCH, [oalig.LOG_FILE])
            oalig.show_print("  Mismatch penalty: %s" % oalig.MISMATCH_PENALTY, [oalig.LOG_FILE])
        elif oalig.SEQUENCE_TYPE == oalig.TYPE_SEQ_AA:
            oalig.init_substitution_matrix()
            oalig.show_print("  Matrix: %s" % oalig.MATRIX_TYPE, [oalig.LOG_FILE])

        oalig.show_print("  Gap penalty: %s" % oalig.GAP_PENALTY, [oalig.LOG_FILE])
        oalig.show_print("", [oalig.LOG_FILE])

        dict_sequences = oalig.get_sequences()
        oalig.show_print("Alignment:", [oalig.LOG_FILE])

        alignment1, alignment2, score, matrix = oalig.smith_waterman(oalig.SEQUENCE1, oalig.SEQUENCE2)
        alignment = oalig.get_alignment_characters(alignment1, alignment2)
        oalig.save_matrix(matrix)

        oalig.show_print("  Score: %s" % score, [oalig.LOG_FILE])
        oalig.show_print("", [oalig.LOG_FILE])
        for part in alignment:
            for align in part:
                oalig.show_print("  %s" % align, [oalig.LOG_FILE])
            oalig.show_print("", [oalig.LOG_FILE])

        oalig.show_print("Matrix file: %s" % oalig.MATRIX_FILE, [oalig.LOG_FILE])
        oalig.show_print("", [oalig.LOG_FILE])

        oalig.show_print(oalig.finish_time(start, "Elapsed time"), [oalig.LOG_FILE])
        oalig.show_print("Done.", [oalig.LOG_FILE])
    except Exception as e:
        oalig.show_print("\n%s" % traceback.format_exc(), [oalig.LOG_FILE])
        oalig.show_print(oalig.finish_time(start, "Elapsed time"), [oalig.LOG_FILE])
        oalig.show_print("Done.", [oalig.LOG_FILE])

if __name__ == '__main__':
    oalig = SmithWaterman()
    main(sys.argv)
