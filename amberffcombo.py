#!/usr/bin/env python
'''
@author: Dazhi Tan
@created: September 19th, 2017
'''
import os, re
import pandas as pd

from StringIO import StringIO

def _check_atomtype_match(df, input_types, labels):
    #assert len(input_types) == len(labels), \
    #    "The number of input_types has to be the same as that of labels you" \
    #     " want to match"
    intypesf = input_types
    intypesr = list(reversed(input_types))
    idx_f = df.index
    idx_r = df.index

    for itf, l in zip(intypesf, labels):
        idx_match = df.index[df[l] == itf]
        idx_f = idx_f.intersection(idx_match)

    for itr, l in zip(intypesr, labels):
        idx_match = df.index[df[l] == itr]
        idx_r = idx_r.intersection(idx_match)

    return idx_f, idx_r

class AmberDat():
    """
    A class describing the .dat file in an AmberFFCombo.
    """
    def __init__(self, dat_file_path=None):
        if dat_file_path:
            self.datpath = dat_file_path
            self.loadDat()
        else:
            self.createDat()
    
    def __repr__(self):
        return "<AmberDat '%s'>" % self.title


    def createDat(self):
        self.title = 'New AmberDat'
        self.mass_df = pd.DataFrame(columns=['type', 'mass', 'pol'])
        self.stretch_df = pd.DataFrame(columns=['type1', 'type2', 'fc', 'r0'])
        self.angle_df = pd.DataFrame(columns=['type1', 'type2', 'type3', 'fc',
                                              'theta0'])
        self.proper_df = pd.DataFrame(columns=['type1', 'type2', 'type3', 
                                               'type4', 'divider', 'fc', 
                                               'theta0', 'periodicity'])
        self.improper_df = pd.DataFrame(columns=['type1', 'type2', 'type3', 
                                                 'type4', 'fc', 
                                                 'theta0', 'periodicity'])
        self.vdw_df = pd.DataFrame(columns=['type', 'half_rmin', 'epsilon'])

        self.hydrophilic_line = \
            'C   H   HO  N   NA  NB  NC  N2  NT  N2  N3  N*  O   OH  OS  P   O2 '        
        self.fast_water_flag = \
            '  HW  OW  0000.     0000.                                4.  flag for fast water'
        self.equivalent_vdwtype = \
            'N   NA  N2  N*  NC  NB  NT  NY\nC*  CA  CB  CC  CD  CK  CM  CN  CQ  CR  CV  CW  CY  CZ  CP  CS'
        self.mod = 'MOD4      RE'
        self.end = 'END'

    def loadDat(self):
        with open(self.datpath, 'r') as fh:
            raw_string = fh.read()
        blocklist = raw_string.split('\n\n')

        assert len(blocklist) == 9, \
                "Please check the format of the .dat file."

        # Load title and mass
        self.title = blocklist[0].split('\n')[0] # First line
        mass_block_csv = StringIO(blocklist[0])
        self.mass_df = pd.read_csv(mass_block_csv, sep='\s+', 
                                   skip_blank_lines=True,
                                   keep_default_na=False,
                                   header=None, names=['type', 'mass', 'pol'],
                                   usecols=[0, 1, 2], skiprows=[0])
        self.mass_df.replace(r'\s*([^\s]+)\s*', r'\1', inplace=True, regex=True)

        # Load stretch terms
        self.hydrophilic_line = blocklist[1].split('\n')[0]
        stretch_block_csv = StringIO(blocklist[1])
        self.stretch_df = pd.read_csv(stretch_block_csv, sep='\s\s+|-', 
                                      engine='python', header=None,
                                      skip_blank_lines=True,
                                      keep_default_na=False,
                                      names=['type1', 'type2', 'fc', 'r0'],
                                      usecols=[0, 1, 2, 3], skiprows=[0])
        self.stretch_df.replace(r'\s*([^\s]+)\s*', r'\1', inplace=True, regex=True)

        # Load angle terms
        angle_block_csv = StringIO(blocklist[2])
        self.angle_df = pd.read_csv(angle_block_csv, sep='\s\s+|-', 
                                    engine='python', header=None,
                                    skip_blank_lines=True,
                                    keep_default_na=False,
                                    names=['type1', 'type2', 'type3', 'fc', 'theta0'],
                                    usecols=[0, 1, 2, 3, 4])
        self.angle_df.replace(r'\s*([^\s]+)\s*', r'\1', inplace=True, regex=True)

        # Load proper torsion terms
        proper_block_csv = StringIO(blocklist[3])
        data_df = pd.read_csv(proper_block_csv, sep='\s\s+', 
                              engine='python', header=None,
                              skip_blank_lines=True,
                              keep_default_na=False,
                              names=['divider', 'fc', 
                                     'theta0','periodicity'],
                              usecols=[1, 2, 3, 4])

        proper_block_csv = StringIO(blocklist[3])
        type_df = pd.read_csv(proper_block_csv, sep='\s\s+|-', header=None,
                              engine='python',
                              skip_blank_lines=True, keep_default_na=False,
                              names=['type1', 'type2', 'type3', 'type4'],
                              usecols=[0, 1, 2, 3])
        self.proper_df = pd.concat([type_df, data_df], axis=1)
        self.proper_df.replace(r'\s*([^\s]+)\s*', r'\1', inplace=True, regex=True)

        # Load improper torsion terms
        improper_block_csv = StringIO(blocklist[4])
        data_df = pd.read_csv(improper_block_csv, sep='\s\s+', 
                              engine='python', header=None,
                              skip_blank_lines=True,
                              keep_default_na=False,
                              names=['fc', 'theta0', 'periodicity'],
                              usecols=[1, 2, 3])

        improper_block_csv = StringIO(blocklist[4])
        type_df = pd.read_csv(improper_block_csv, sep='\s\s+|-', 
                              engine='python', header=None,
                              skip_blank_lines=True,
                              keep_default_na=False,
                              names=['type1', 'type2', 'type3', 'type4'],
                              usecols=[0, 1, 2, 3])

        self.improper_df = pd.concat([type_df, data_df], axis=1)
        self.improper_df.replace(r'\s*([^\s]+)\s*', r'\1', inplace=True, regex=True)

        # Load other stuff
        self.fast_water_flag = blocklist[5]
        self.equivalent_vdwtype = blocklist[6]

        # Load vdw terms
        self.mod = blocklist[7].split('\n')[0]
        vdw_block_csv = StringIO(blocklist[7])
        self.vdw_df = pd.read_csv(vdw_block_csv, sep='\s+', 
                                  engine='python', header=None,
                                  skip_blank_lines=True,
                                  keep_default_na=False,
                                  names=['type', 'half_rmin', 'epsilon'],
                                  usecols=[0, 1, 2], skiprows=[0])
        #self.vdw_df['sigma'] = self.vdw_df.half_rmin/(2**(1./6))*2.
        self.vdw_df.replace(r'\s*([^\s]+)\s*', r'\1', inplace=True, regex=True)

        # Load last line
        self.end = blocklist[8]

    def printDat(self, out_file_path):
        out_fh = open(out_file_path, 'w')
        print >> out_fh, self.title

        # print mass block
        for n in range(len(self.mass_df)):
            row = self.mass_df.iloc[n]
            if n == 0:
                print >> out_fh, "%-2s%6s%14.3f  !" % tuple(row)
            else:
                print >> out_fh, "%-2s%6s%14.3f" % tuple(row)

        print >> out_fh, ''
        print >> out_fh, self.hydrophilic_line

        # print stretch block
        for n in range(len(self.stretch_df)):
            row = self.stretch_df.iloc[n]
            if n == 0:
                print >> out_fh, "%-2s-%-2s%7.1f%9.3f     !" % tuple(row)
            else:
                print >> out_fh, "%-2s-%-2s%7.1f%9.3f" % tuple(row)

        print >> out_fh, ''

        # print angle block
        for n in range(len(self.angle_df)):
            row = self.angle_df.iloc[n]
            print >> out_fh, "%-2s-%-2s-%-2s%8.1f%12.2f" % tuple(row)

        print >> out_fh, ''

        # print proper block
        for n in range(len(self.proper_df)):
            row = self.proper_df.iloc[n]
            print >> out_fh, "%-2s-%-2s-%-2s-%-2s%4d%8.3f%13.1f%14.0f." % tuple(row)

        print >> out_fh, ''

        # print improper block
        for n in range(len(self.improper_df)):
            row = self.improper_df.iloc[n]
            print >> out_fh, "%-2s-%-2s-%-2s-%-2s         %-13.1f%-14s%.0f." % tuple(row)

        print >> out_fh, ''
        print >> out_fh, self.fast_water_flag
        print >> out_fh, ''
        print >> out_fh, self.equivalent_vdwtype
        print >> out_fh, ''
        print >> out_fh, self.mod

        # print vdw block
        for n in range(len(self.vdw_df)):
            row = self.vdw_df.iloc[n]
            if n == 0:
                print >> out_fh, "%4s%16.4f%8.4f            !" % tuple(row)
            else:
                print >> out_fh, "%4s%16.4f%8.4f" % tuple(row)

        print >> out_fh, ''
        print >> out_fh, self.end

    def setTitle(self, new_title):
        self.title = new_title

    def addMass(self, atomtype, new_mass):
        assert new_mass >= 0, \
            'Mass of an atom cannot be negative.'
        atomtype = atomtype.replace(' ', '')
        assert 0 < len(atomtype) < 3, \
            'Atomtype must be an non-empty string with no more than two characters.'

        idx, _ = _check_atomtype_match(self.mass_df, [atomtype], ['type'])
        if len(idx) == 1:
            print "Atomtype %s already exists in %s, no mass added." \
                    % (atomtype, self)
        elif len(idx) == 0:
            n = len(self.mass_df)
            self.mass_df.loc[n] = [atomtype, new_mass, 0.000]
            print "Atomtype %s with mass %.3f added to %s." \
                    % (atomtype, new_mass, self)
        else:
            print "Something is wrong with %s, no mass set." \
                    % self
    def setMass(self, atomtype, new_mass):
        assert new_mass >= 0, \
            'Mass of an atom cannot be negative.'

        idx, _ = _check_atomtype_match(self.mass_df, [atomtype], ['type'])
        if len(idx) == 1:
            self.mass_df.loc[idx, 'mass'] = new_mass
            print "Mass of atomtype %s has been set to %.3f." \
                    % (atomtype, new_mass)
        elif len(idx) == 0:
            print "Atomtype %s does not exist in %s, no mass set." \
                    % (atomtype, self)
        else:
            print "Something is wrong with %s, no mass set." \
                    % self

    def addStretch(self, atomtype1, atomtype2, new_fc, new_r0):
        assert new_fc >= 0, \
            "Force constant (fc) of a stretch term cannot be negative."
        atomtype1 = atomtype1.replace(' ', '')
        atomtype2 = atomtype2.replace(' ', '')
        assert 0 < len(atomtype1) < 3 and 0 < len(atomtype2) < 3, \
            'Atomtype must be an non-empty string with no more than two characters.'

        idf, idr = _check_atomtype_match(self.stretch_df, 
                       [atomtype1, atomtype2], ['type1', 'type2'])

        if len(idf) == 1 and len(idr) == 0:
            print "Stretch term %-2s-%-2s already exists in %s, no term added." \
                    % (atomtype1, atomtype2, self)

        elif len(idf) == 0 and len(idr) == 1:
            print "Stretch term %-2s-%-2s already exists in %s, no term added." \
                    % (atomtype2, atomtype1, self)

        elif len(idf) == 0 and len(idr) == 0:
            n = len(self.stretch_df)
            self.stretch_df.loc[n] = [atomtype1, atomtype2, new_fc, new_r0]
            print "Stretch term %-2s-%-2s with fc %.1f and r0 %.3f added to %s," \
                  % (atomtype1, atomtype2, new_fc, new_r0, self)

        else:
            print "Something is wrong with %s, no stretch term set." \
                    % self

    def setStretch(self, atomtype1, atomtype2, new_fc, new_r0):
        assert new_fc >= 0, \
            "Force constant (fc) of a stretch term cannot be negative."

        idf, idr = _check_atomtype_match(self.stretch_df, 
                       [atomtype1, atomtype2], ['type1', 'type2'])
        if len(idf) == 1 and len(idr) == 0:
            self.stretch_df.loc[idf, 'fc'] = new_fc
            self.stretch_df.loc[idf, 'r0'] = new_r0
            print "fc of stretch term %-2s-%-2s has been set to %.1f." \
                    % (atomtype1, atomtype2, new_fc)
            print "r0 of stretch term %-2s-%-2s has been set to %.3f." \
                    % (atomtype1, atomtype2, new_r0)

        elif len(idf) == 0 and len(idr) == 1:
            self.stretch_df.loc[idr, 'fc'] = new_fc
            self.stretch_df.loc[idr, 'r0'] = new_r0
            print "fc of stretch term %-2s-%-2s has been set to %.1f." \
                    % (atomtype2, atomtype1, new_fc)
            print "r0 of stretch term %-2s-%-2s has been set to %.3f." \
                    % (atomtype2, atomtype1, new_r0)

        elif len(idf) == 0 and len(idr) == 0:
            print "Stretch term %-2s-%-2s or %-2s-%-2s does not exist in %s," \
                  % (atomtype1, atomtype2, atomtype2, atomtype1, self),
            print "no stretch term set." 

        else:
            print "Something is wrong with %s, no stretch term set." \
                    % self

    def addAngle(self, atomtype1, atomtype2, atomtype3, new_fc, new_theta0):
        assert new_fc >= 0, \
            "Force constant (fc) of a angle term cannot be negative."
        atomtype1 = atomtype1.replace(' ', '')
        atomtype2 = atomtype2.replace(' ', '')
        atomtype3 = atomtype3.replace(' ', '')
        assert 0 < len(atomtype1) < 3 and 0 < len(atomtype2) < 3 and \
            0 < len(atomtype3) < 3, \
            'Atomtype must be an non-empty string with no more than two characters.'

        idf, idr = _check_atomtype_match(self.angle_df, 
                       [atomtype1, atomtype2, atomtype3], 
                       ['type1', 'type2', 'type3'])

        if len(idf) == 1 and len(idr) == 0:
            print "Angle term %-2s-%-2s-%-2s already exists in %s, no term added." \
                    % (atomtype1, atomtype2, atomtype3, self)

        elif len(idf) == 0 and len(idr) == 1:
            print "Angle term %-2s-%-2s-%-2s already exists in %s, no term added." \
                    % (atomtype3, atomtype2, atomtype1, self)

        elif len(idf) == 0 and len(idr) == 0:
            n = len(self.angle_df)
            self.angle_df.loc[n] = [atomtype1, atomtype2, atomtype3, new_fc, new_theta0]
            print "Angle term %-2s-%-2s-%-2s with fc %.1f and theta0 %.3f added to %s," \
                  % (atomtype1, atomtype2, atomtype3, new_fc, new_theta0, self)

        else:
            print "Something is wrong with %s, no angle term set." \
                    % self

    def setAngle(self, atomtype1, atomtype2, atomtype3, new_fc, new_theta0):
        assert new_fc >= 0, \
            "Force constant (fc) of an angle term cannot be negative."

        idf, idr = _check_atomtype_match(self.angle_df, 
                       [atomtype1, atomtype2, atomtype3], 
                       ['type1', 'type2', 'type3'])

        if len(idf) == 1 and len(idr) == 0:
            self.angle_df.loc[idf, 'fc'] = new_fc
            self.angle_df.loc[idf, 'theta0'] = new_theta0
            print "fc of angle term %-2s-%-2s-%-2s has been set to %.1f." \
                    % (atomtype1, atomtype2, atomtype3, new_fc)
            print "theta0 of angle term %-2s-%-2s-%-2s has been set to %.2f." \
                    % (atomtype1, atomtype2, atomtype3, new_theta0)

        elif len(idf) == 0 and len(idr) == 1:
            self.angle_df.loc[idr, 'fc'] = new_fc
            self.angle_df.loc[idr, 'theta0'] = new_theta0
            print "fc of angle term %-2s-%-2s-%-2s has been set to %.1f." \
                    % (atomtype3, atomtype2, atomtype1, new_fc)
            print "theta0 of angle term %-2s-%-2s-%-2s has been set to %.2f." \
                    % (atomtype3, atomtype2, atomtype1, new_theta0)

        elif len(idf) == 0 and len(idr) == 0:
            print "Angle term %-2s-%-2s-%-2s or %-2s-%-2s-%-2s does not exist in %s," \
                  % (atomtype1, atomtype2, atomtype3, 
                     atomtype3, atomtype2, atomtype1, self),
            print "no angle term set." 

        else:
            print "Something is wrong with %s, no angle term set." \
                    % self

    def addProper(self, atomtype1, atomtype2, atomtype3, atomtype4, 
                  new_divider, new_fc, new_theta0, new_periodicity):
        
        assert new_periodicity in [x for x in range(-6, 7) if x != 0], \
                "periodicity of a proper torsion term has to be a non-zero " \
                "integer between -6 and 6"
        atomtype1 = atomtype1.replace(' ', '')
        atomtype2 = atomtype2.replace(' ', '')
        atomtype3 = atomtype3.replace(' ', '')
        atomtype4 = atomtype4.replace(' ', '')
        assert 0 < len(atomtype1) < 3 and 0 < len(atomtype2) < 3 and \
            0 < len(atomtype3) < 3 and 0 < len(atomtype4) < 3, \
            'Atomtype must be an non-empty string with no more than two characters.'

        n = len(self.proper_df)
        self.proper_df.loc[n] = [atomtype1, atomtype2, atomtype3, atomtype4,
                                 new_divider, new_fc, new_theta0, 
                                 new_periodicity]
        
        self._refreshProper(atomtype1, atomtype2, atomtype3, atomtype4)
        print "Proper term %-2s-%-2s-%-2s-%-2s with divider %d fc %.3f" \
              % (atomtype1, atomtype2, atomtype3, atomtype4, 
                 new_divider, new_fc),
        print "theta0 %.1f and periodicity %d added. " \
              % (new_theta0, new_periodicity)

    def _refreshProper(self, atomtype1, atomtype2, atomtype3, atomtype4):
        idf, idr = _check_atomtype_match(self.proper_df, 
                       [atomtype1, atomtype2, atomtype3, atomtype4], 
                       ['type1', 'type2', 'type3', 'type4'])

        if len(idf) > 0 and len(idr) == 0:
            for i in idf[:-1]:
                if self.proper_df.loc[i, 'periodicity'] > 0:
                    self.proper_df.loc[i, 'periodicity'] *= -1
            if self.proper_df.loc[idf[-1], 'periodicity'] < 0:
                self.proper_df.loc[idf[-1], 'periodicity'] *= -1

        elif len(idf) == 0 and len(idr) > 0:
            for i in idr[:-1]:
                if self.proper_df.loc[i, 'periodicity'] > 0:
                    self.proper_df.loc[i, 'periodicity'] *= -1
            if self.proper_df.loc[idr[-1], 'periodicity'] < 0:
                self.proper_df.loc[idr[-1], 'periodicity'] *= -1

        elif len(idf) == 0 and len(idr) == 0:
            pass

        else:
            print "Something is wrong with %s, no proper torsion term set." \
                    % self

    def addImproper(self, atomtype1, atomtype2, atomtype3, atomtype4, 
                    new_fc, new_theta0, new_periodicity):
        assert new_fc >= 0 and new_theta0 >= 0, \
            "fc and theta0 of an improper term cannot be negative."
        assert new_periodicity in [x for x in range(-6, 7) if x != 0], \
                "periodicity of an improper torsion term has to be a non-zero " \
                "integer between -6 and 6"

        atomtype1 = atomtype1.replace(' ', '')
        atomtype2 = atomtype2.replace(' ', '')
        atomtype3 = atomtype3.replace(' ', '')
        atomtype4 = atomtype4.replace(' ', '')
        assert 0 < len(atomtype1) < 3 and 0 < len(atomtype2) < 3 and \
            0 < len(atomtype3) < 3 and 0 < len(atomtype4) < 3, \
            'Atomtype must be an non-empty string with no more than two characters.'

        idf, _ = _check_atomtype_match(self.improper_df, 
                       [atomtype1, atomtype2, atomtype3, atomtype4], 
                       ['type1', 'type2', 'type3', 'type4'])

        if len(idf) == 1:
            print "Improper term %-2s-%-2s-%-2s-%-2s already exists in %s, no term added." \
                    % (atomtype1, atomtype2, atomtype3, atomtype4, self)

        elif len(idf) == 0:
            n = len(self.improper_df)
            self.improper_df.loc[n] = [atomtype1, atomtype2, atomtype3, atomtype4, 
                                       new_fc, new_theta0, new_periodicity]
            print "Improper term %-2s-%-2s-%-2s-%-2s with fc %.1f theta0 %.1f periodicity %d added to %s," \
                  % (atomtype1, atomtype2, atomtype3, atomtype4, new_fc, 
                     new_theta0, new_periodicity, self)

        else:
            print "Something is wrong with %s, no improper term set." \
                    % self

    def setImproper(self, atomtype1, atomtype2, atomtype3, atomtype4, 
                    new_fc, new_theta0, new_periodicity):
        assert new_fc >= 0 and new_theta0 >= 0, \
            "fc and theta0 of an improper term cannot be negative."
        assert new_periodicity in [x for x in range(-6, 7) if x != 0], \
                "periodicity of an improper torsion term has to be a non-zero " \
                "integer between -6 and 6"

        atomtype1 = atomtype1.replace(' ', '')
        atomtype2 = atomtype2.replace(' ', '')
        atomtype3 = atomtype3.replace(' ', '')
        atomtype4 = atomtype4.replace(' ', '')
        assert 0 < len(atomtype1) < 3 and 0 < len(atomtype2) < 3 and \
            0 < len(atomtype3) < 3 and 0 < len(atomtype4) < 3, \
            'Atomtype must be an non-empty string with no more than two characters.'

        idf, _ = _check_atomtype_match(self.improper_df, 
                       [atomtype1, atomtype2, atomtype3, atomtype4], 
                       ['type1', 'type2', 'type3', 'type4'])

        if len(idf) == 1:
            self.improper_df.loc[idf, 'fc'] = new_fc
            self.improper_df.loc[idf, 'theta0'] = new_theta0
            self.improper_df.loc[idf, 'periodicity'] = new_periodicity
            print "fc of improper term %-2s-%-2s-%-2s-%-2s set to %.1f" \
                    % (atomtype1, atomtype2, atomtype3, atomtype4, new_fc)
            print "theta0 of improper term %-2s-%-2s-%-2s-%-2s set to %.1f" \
                    % (atomtype1, atomtype2, atomtype3, atomtype4, new_theta0)
            print "periodicity of improper term %-2s-%-2s-%-2s-%-2s set to %d" \
                    % (atomtype1, atomtype2, atomtype3, atomtype4, new_periodicity)

        elif len(idf) == 0:
            print "Improper term %-2s-%-2s-%-2s-%-2s or %-2s-%-2s-%-2s-%-2s does not exist in %s," \
                  % (atomtype1, atomtype2, atomtype3, atomtype4,
                     atomtype4, atomtype3, atomtype2, atomtype1, self),
            print "no improper term set." 

        else:
            print "Something is wrong with %s, no improper term set." \
                    % self

    def setVdw(self, atomtype, new_half_rmin, new_epsilon):
        if atomtype in self.vdw_df.type.values:
            index = self.vdw_df.index[self.vdw_df.type == atomtype]
            self.vdw_df.loc[index, 'half_rmin'] = new_half_rmin
            self.vdw_df.loc[index, 'epsilon'] = new_epsilon
            print "Half Rmin of atomtype %s has been set to %.4f" \
                    % (atomtype, new_half_rmin)
            print "Epsilon of atomtype %s has been set to %.4f" \
                    % (atomtype, new_epsilon)
        else:
            print "Atomtype %s does not exist in %s, no vdw params set" \
                    % (atomtype, self)

class AmberLib():

    def __init__(self, input_lib_file):
        self.inputfile = input_lib_file
        self.LoadLib()

    def LoadLib(self):
        self.atoms_dict = dict()
        self.atomspertinfo_dict = dict()
        entry_pattern = re.compile("entry\.(\w+)\.unit\.(\w+)")

        with open(self.inputfile, 'r') as fh:
            raw_string = fh.read()
            self.work_list = raw_string.split('!')

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
