#!/usr/bin/env python
'''
@author: Dazhi Tan
@created: September 19th, 2017
'''
import os, re
import pandas as pd

class AmberLib():

    def __init__(self, input_lib_file):
        self.inputfile = input_lib_file
        self.LoadLib()

    def LoadLib(self):
        self.atoms_dict = dict()
        self.atomspertinfo_dict = dict()
        entry_pattern = re.compile("entry\.(\w+)\.unit\.(\w+)")

        with open(self.inputfile, 'r') as fh:
            self._raw_string = fh.read()
            self.work_list = self._raw_string.split('!')

        self.work_list = filter(None, self.work_list) # Filter out empty strings

        self.atoms_edit_idx = []
        self.atomspertinfo_edit_idx = []
        for n, block in enumerate(self.work_list):
            entry_match = entry_pattern.match(block)
            if entry_match:
                self.work_list[n] = '!' + block
                resi_name, entry_name = entry_match.groups()
                if entry_name == 'atoms':
                    self.atoms_edit_idx.append(n)
                    self.atoms_dict[resi_name] = pd.DataFrame(
                        columns=['name', 'type', 'typex', 'resx', 'flag', 
                                 'seq', 'elmnt', 'charge'])
                    atom_list = block.split('\n')[1:-1]
                    for m, al in enumerate(atom_list):
                        self.atoms_dict[resi_name].loc[m] = al.split(' ')[1:]

                elif entry_name == 'atomspertinfo':
                    self.atomspertinfo_edit_idx.append(n)
                    self.atomspertinfo_dict[resi_name] = pd.DataFrame(
                        columns=['pname', 'ptype', 'ptypex', 'pelmnt', 'pchg'])
                    atom_list = block.split('\n')[1:-1]
                    for m, al in enumerate(atom_list):
                        self.atomspertinfo_dict[resi_name].loc[m] = al.split(' ')[1:]
                else:
                    continue
            else:
                self.work_list[n] = '!!' + block

    def PrintLib(self, output_fh):
        for n, block in enumerate(self.work_list):
            if n in self.atoms_edit_idx:
                output_lines = block.split('\n')[:-1]
                print >> output_fh, output_lines[0]
                resname = re.search("entry\.(\w+)\.unit", block).groups()[0]
                resi_params = self.atoms_dict[resname]
                for _, row in resi_params.iterrows():
                    print >> output_fh, ' %s %s %s %s %s %s %s %s' % tuple(row)
            elif n in self.atomspertinfo_edit_idx:
                output_lines = block.split('\n')[:-1]
                print >> output_fh, output_lines[0]
                resname = re.search("entry\.(\w+)\.unit", block).groups()[0]
                resi_params = self.atomspertinfo_dict[resname]
                for _, row in resi_params.iterrows():
                    print >> output_fh, ' %s %s %s %s %s' % tuple(row)
                pass
            else:
                print >> output_fh, block, 

    def SetAtomType(self, resi_name, atom_name, new_type):
        if resi_name in self.atoms_dict.keys():
            atoms_df = self.atoms_dict[resi_name]
            atomspertinfo_df = self.atomspertinfo_dict[resi_name]

            i1 = atoms_df[atoms_df.name==atom_name].index
            i2 = atomspertinfo_df[atomspertinfo_df.pname==atom_name].index

            atoms_df.loc[i1, 'type'] = new_type
            atomspertinfo_df.loc[i2, 'ptype'] = new_type
        else:
            print '%s does not exist in %s' % (resi_name, self.inputfile)
            pass

    def SetAtomCharge(self, resi_name, atom_name, new_charge):
        if resi_name in self.atoms_dict.keys():

            atoms_df = self.atoms_dict[resi_name]

            i = atoms_df[atoms_df.name==atom_name].index
            atoms_df.loc[i, 'charge'] = new_charge
        else:
            print '%s does not exist in %s' % (resi_name, self.inputfile)
            pass
