#!/usr/bin/env python
'''
@author: Dazhi Tan
@created: September 19th, 2017
'''
import os, re
import pandas as pd

from UserDict import IterableUserDict
from StringIO import StringIO

class AmberLib(IterableUserDict):
    """
    A class describing the .lib file in AmberFFCombo.
    """
    def __init__(self, lib_file_path):
        IterableUserDict.__init__(self)
        self.libpath = lib_file_path
        self.title = os.path.basename(self.libpath)
        self.loadLib()

    def __repr__(self):
        return '<AmberLib %s>' % self.title

    def loadLib(self):
        with open(self.libpath, 'r') as fh:
            raw_string = fh.read()
        blocklist = raw_string.split('!')
        blocklist = filter(None, blocklist)
        self.blocklist = blocklist

        self.residue_list = blocklist[0].split('\n')[1:-1]
        self.residue_list = [x.replace(' ', '') for x in self.residue_list]
        self.residue_list = [x.replace('"', '') for x in self.residue_list]

        for resname in self.residue_list:
            self.data[resname] = dict()

        header_pattern = re.compile("entry\.(\w+)\.unit\.(\w+)")
        for block in blocklist[1:]:
            resname, entryname = \
                header_pattern.match(block).groups()
            temp_list = block.split('\n')[0].split(' ')
            if temp_list[1] == 'table':
                header_string_list = block.split('\n')[0].split('  ')
                #data_names = [re.search('\s(\w+)', x).groups()[0] for x in
                #              header_string_list[1:]]
                block_csv = StringIO(block)
                info_df = pd.read_csv(block_csv, sep='\s+', skiprows=[0], 
                                      keep_default_na=False,
                                      names=header_string_list[1:])
                self.data[resname][entryname] = info_df
            if temp_list[1] == 'array' or temp_list[1] == 'single':
                self.data[resname][entryname] = block.split('\n')[:-1]

    def printLib(self, out_file_path):
        fh = open(out_file_path, 'w')
        print >> fh, '!!index array str'
        for r in self.residue_list:
            print >> fh, ' "%s"' % r

        for r in self.residue_list:

            print >> fh, '!entry.%s.unit.atoms table' % r,
            for c in self.data[r]['atoms'].columns:
                print >> fh, ' %s' % c, 
            print >> fh, ''
            for _, row in self.data[r]['atoms'].iterrows():
                print >> fh, ' "%s" "%s" %d %d %d %d %d %f' % tuple(row)

            print >> fh, '!entry.%s.unit.atomspertinfo table' % r,
            for c in self.data[r]['atomspertinfo'].columns:
                print >> fh, ' %s' % c, 
            print >> fh, ''
            for _, row in self.data[r]['atomspertinfo'].iterrows():
                print >> fh, ' "%s" "%s" %d %d %.1f' % tuple(row)

            for n, elmnt in enumerate(self.data[r]['boundbox']):
                if n == 0:
                    print >> fh, '!%s' % elmnt
                elif n > 1:
                    print >> fh, ' %.1f' % float(elmnt) 
                else:
                    print >> fh, ' %f' % float(elmnt) 

            for n, elmnt in enumerate(self.data[r]['childsequence']):
                if n == 0:
                    print >> fh, '!%s' % elmnt
                else:
                    print >> fh, ' %d' % int(elmnt) 

            for n, elmnt in enumerate(self.data[r]['connect']):
                if n == 0:
                    print >> fh, '!%s' % elmnt
                else:
                    print >> fh, ' %d' % int(elmnt) 

            print >> fh, '!entry.%s.unit.connectivity table' % r,
            for c in self.data[r]['connectivity'].columns:
                print >> fh, ' %s' % c, 
            print >> fh, ''
            for _, row in self.data[r]['connectivity'].iterrows():
                print >> fh, ' %d %d %d' % tuple(row)

            print >> fh, '!entry.%s.unit.hierarchy table' % r,
            for c in self.data[r]['hierarchy'].columns:
                print >> fh, ' %s' % c, 
            print >> fh, ''
            for _, row in self.data[r]['hierarchy'].iterrows():
                print >> fh, ' "%s" %d "%s" %d' % tuple(row)

            for n, elmnt in enumerate(self.data[r]['name']):
                if n == 0:
                    print >> fh, '!%s' % elmnt
                else:
                    print >> fh, ' %s' % elmnt.replace(' ', '') 

            print >> fh, '!entry.%s.unit.positions table' % r,
            for c in self.data[r]['positions'].columns:
                print >> fh, ' %s' % c, 
            print >> fh, ''
            for _, row in self.data[r]['positions'].iterrows():
                print >> fh, ' %f %f %.6e' % tuple(row)

            print >> fh, '!entry.%s.unit.residueconnect table' % r,
            for c in self.data[r]['residueconnect'].columns:
                print >> fh, ' %s' % c, 
            print >> fh, ''
            for _, row in self.data[r]['residueconnect'].iterrows():
                print >> fh, ' %d %d %d %d %d %d' % tuple(row)

            print >> fh, '!entry.%s.unit.residues table' % r,
            for c in self.data[r]['residues'].columns:
                print >> fh, ' %s' % c, 
            print >> fh, ''
            for _, row in self.data[r]['residues'].iterrows():
                print >> fh, ' "%s" %d %d %d "%s" %d' % tuple(row)

            for n, elmnt in enumerate(self.data[r]['residuesPdbSequenceNumber']):
                if n == 0:
                    print >> fh, '!%s' % elmnt
                else:
                    print >> fh, ' %d' % int(elmnt) 

            for n, elmnt in enumerate(self.data[r]['solventcap']):
                if n == 0:
                    print >> fh, '!%s' % elmnt
                elif n > 2:
                    print >> fh, ' %.1f' % float(elmnt) 
                else:
                    print >> fh, ' %f' % float(elmnt) 

            print >> fh, '!entry.%s.unit.velocities table' % r,
            for c in self.data[r]['velocities'].columns:
                print >> fh, ' %s' % c, 
            print >> fh, ''
            for _, row in self.data[r]['velocities'].iterrows():
                print >> fh, ' %.1f %.1f %.1f' % tuple(row)

        print "Printed %s to %s" % (self, out_file_path)

    def setAtomType(self, resname, atomname, newtype):
        resname = resname.replace(' ', '')
        atomname = atomname.replace(' ', '')
        newtype = newtype.replace(' ', '')
        assert 0 < len(newtype) < 3, \
            "Atomtype must be a non-empty string with no more than two characters."
        if resname in self.data.keys():
            atm_df = self.data[resname]['atoms']
            atm_idx = atm_df.index[atm_df['str name'] == atomname]
            if len(atm_idx) == 0:
                print '%s does not exist in residue %s in %s, no atomtype set' \
                    % (atomname, resname, self)
                return
            atm_df.loc[atm_idx, 'str type'] = newtype

            atmpert_df = self.data[resname]['atomspertinfo']
            atmpert_idx = atmpert_df.index[atmpert_df['str pname'] == atomname]
            atmpert_df.loc[atmpert_idx, 'str ptype'] = newtype
            print 'The type of %s in residue %s is set to be %s' \
                % (atomname, resname, newtype)
        else:
            print '%s does not exist in %s, no atomtype set' % (resname, self)

    def setAtomCharge(self, resname, atomname, newcharge):
        resname = resname.replace(' ', '')
        atomname = atomname.replace(' ', '')
        assert isinstance(newcharge, (float, int)), \
            "Charge must be a number."
        if resname in self.data.keys():
            atm_df = self.data[resname]['atoms']
            atm_idx = atm_df.index[atm_df['str name'] == atomname]
            if len(atm_idx) == 0:
                print '%s does not exist in residue %s in %s, no charge set' \
                    % (atomname, resname, self)
                return
            atm_df.loc[atm_idx, 'dbl chg'] = newcharge

            print 'The charge of %s in residue %s is set to be %f' \
                % (atomname, resname, newcharge)
        else:
            print '%s does not exist in %s, no atomtype set' % (resname, self)

    def calcNetCharge(self, resname):
        resname = resname.replace(' ', '')
        if resname in self.data.keys():
            atm_df = self.data[resname]['atoms']
            net_charge = atm_df['dbl chg'].sum()
            return net_charge
        else:
            print '%s does not exist in %s, no net charge calculated' % (resname, self)
            return None
