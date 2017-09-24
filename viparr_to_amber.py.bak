#!/usr/bin/env python
'''
@author: Dazhi Tan
@created: July 14th, 2017
'''
import os, re, json, argparse
from amberlib import AmberLib
from helper import load_viparr_dir
import pandas as pd

def viparr_template_to_amberlib(viparr_dir, in_amberlib, out_amberlib):

    viparr_json_dict = load_viparr_dir(viparr_dir)
    template_amberlib = AmberLib(in_amberlib)
    #-------------------------#
    # Process templates files #
    #-------------------------#
    template_jsons = []
    for key, json in viparr_json_dict.iteritems():
        if 'template' in key:
            template_jsons.append(json)
    for tj in template_jsons:
        for resiname, resiparam in tj.iteritems():
            for atomparam in resiparam['atoms']:
                (atom_name, atomic_number, charge, 
                 [bonded_type, nonbonded_type]) = atomparam

                #----------#
                # Clean up #
                #----------#
                atom_name = atom_name.replace("&", "*")
                bonded_type = bonded_type.replace("&", "*")
                nonbonded_type = nonbonded_type.replace("&", "*")

                atom_name = '"%s"' % atom_name
                bonded_type = '"%s"' % bonded_type
                nonbonded_type = '"%s"' % nonbonded_type
                charge = '%.6f' % charge

                # Set properties
                template_amberlib.SetAtomType(resiname, atom_name, bonded_type) 
                template_amberlib.SetAtomCharge(resiname, atom_name, charge) 

    template_amberlib.PrintLib(open(out_amberlib, 'w'))


def viparr_params_to_amberdat(viparr_dir, out_amberdat):
    viparr_json_dict = load_viparr_dir(viparr_dir)
    fh = open(out_amberdat, 'w')
    print >> fh, 'STX-AMBER'

    # Process mass
    for n, item in enumerate(viparr_json_dict['mass']):
	if not isinstance(item['type'], list):
	    item['type'] = item['type'].split(' ')
        if n == 0:
            print >> fh, "%-2s %5s%14.3f  !            %s" \
                         % (item['type'][0], item['params']['amu'],
                            0.616, item['memo'])
        else:
            print >> fh, "%-2s %5s%14.3f               %s" \
                         % (item['type'][0], item['params']['amu'],
                            0.616, item['memo'])

    #Process bond
    print >> fh, ''
    print >> fh, 'C   H   HO  N   NA  NB  NC  N2  NT  N2  N3  N*  O   OH  OS  P   O2'
    for n, stretch_term in enumerate(viparr_json_dict['stretch_harm']):
        k, length = stretch_term['params']['fc'], stretch_term['params']['r0']
	if not isinstance(stretch_term['type'], list):
		stretch_term['type'] = stretch_term['type'].split(' ')
        type1, type2 = stretch_term['type']
        if n == 0:
            print >>fh, "%-2s-%-2s%7.1f%9.3f     ! STX" % (type1, type2, k, length)
        else:
            print >>fh, "%-2s-%-2s%7.1f%9.3f       STX" % (type1, type2, k, length)

    #Process angle
    print >> fh, ''
    for angle_term in viparr_json_dict['angle_harm']:
        k, theta = angle_term['params']['fc'], angle_term['params']['theta0']
	if not isinstance(angle_term['type'], list):
		angle_term['type'] = angle_term['type'].split(' ')
        type1, type2, type3 = angle_term['type']
        print >>fh, "%-2s-%-2s-%-2s%8.1f%12.2f    %s" % (type1, type2, type3, k, theta, 'STX')

    #Process proper torsion
    print >> fh, ''

    amber_torsion_list = []
    for proper_term in viparr_json_dict['dihedral_trig']:
        # Clean up
	if not isinstance(proper_term['type'], list):
		proper_term['type'] = proper_term['type'].split(' ')
        for n, t in enumerate(proper_term['type']):
            if t == '*':
                proper_term['type'][n] = 'X'

        type1, type2, type3, type4 = proper_term['type']
        phase_angle = proper_term['params']['phi0']

        #-----------------------------------------------#
        # keys are phi0, fc0, fc1, ... fc6. Ignore fc0. #
        #-----------------------------------------------#
        for key, value in proper_term['params'].iteritems():
            if '0' in key : 
                continue  
            #elif value == 0:
            #    continue
            elif value < 0:
                periodicity = int(key[-1])
                k = value
                #amber_torsion_list.append([type1, type2, type3, type4, 1, -k, phase_angle+180, periodicity])
                amber_torsion_list.append([type1, type2, type3, type4, 1, k, phase_angle, periodicity])
            else:
                periodicity = int(key[-1])
                k = value
                amber_torsion_list.append([type1, type2, type3, type4, 1, k, phase_angle, periodicity])

    for m, item in enumerate(amber_torsion_list[:-1]):
        type1, type2, type3, type4, _, k, phase_angle, periodicity = item
        if item[:4] == amber_torsion_list[m+1][:4]:
            print >> fh, "%-2s-%-2s-%-2s-%-2s%4d%8.2f%13.1f%14d.         %s" \
                         % (type1, type2, type3, type4, 1, k, phase_angle, periodicity*-1, 'STX')
        else:
            print >> fh, "%-2s-%-2s-%-2s-%-2s%4d%8.2f%13.1f%14d.         %s" \
                         % (type1, type2, type3, type4, 1, k, phase_angle, periodicity, 'STX')
    print >> fh, "%-2s-%-2s-%-2s-%-2s%4d%9.3f%12.1f%14d." % tuple(amber_torsion_list[-1])


    #Process improper torsion
    print >> fh, ''

    for improper_term in viparr_json_dict['improper_trig']:
        # Clean up
	if not isinstance(improper_term['type'], list):
		improper_term['type'] = improper_term['type'].split(' ')
        for n, t in enumerate(improper_term['type']):
            if t == '*':
                improper_term['type'][n] = 'X'

        type1, type2, type3, type4 = improper_term['type']
        phase_angle = improper_term['params']['phi0']

        #-----------------------------------------------#
        # keys are phi0, fc0, fc1, ... fc6. Ignore fc0. #
        #-----------------------------------------------#
        for key, value in improper_term['params'].iteritems():
            if '0' in key : 
                continue  
            elif value == 0:
                continue
            elif value < 0:
                periodicity = int(key[-1])
                k = value
                #print >> fh, "%-2s-%-2s-%-2s-%-2s         %-4.1f%12d.%11d.           %s" \
                #             % (type1, type2, type3, type4, k, phase_angle, periodicity, 'STX')
                print >> fh, "%-2s-%-2s-%-2s-%-2s         %-4.1f%12d.%11d.           %s" \
                             % (type1, type2, type3, type4, -k, phase_angle+180, periodicity, 'STX')
            else:
                periodicity = int(key[-1])
                k = value
                print >> fh, "%-2s-%-2s-%-2s-%-2s         %-4.1f%12d.%11d.           %s" \
                             % (type1, type2, type3, type4, k, phase_angle, periodicity, 'STX')
    # Process vdw
    print >> fh, ''
    print >> fh, "  HW  OW  0000.     0000.                                4.  flag for fast water"
    print >> fh, ''
    print >> fh, "N   NA  N2  N*  NC  NB  NT  NY"
    print >> fh, "C*  CA  CB  CC  CD  CK  CM  CN  CQ  CR  CV  CW  CY  CZ  CP  CS"
    print >> fh, ''
    print >> fh, "MOD4      RE"
    for n, vdw_term in enumerate(viparr_json_dict['vdw1']):
        # Clean up
        if not isinstance(vdw_term['type'], list):
            vdw_term['type'] = vdw_term['type'].split(' ')
        atom_type = vdw_term['type'][0]
        sigma, epsilon = vdw_term['params']['sigma'], vdw_term['params']['epsilon']
        Rmin_2 = sigma * 2**(1./6) / 2.
        if n == 0:
            print >> fh, "%4s%16.4f%8.4f%13s%s" % (atom_type, Rmin_2, epsilon, '!', 'STX')
        else:
            print >> fh, "%4s%16.4f%8.4f%13s%s" % (atom_type, Rmin_2, epsilon, '', 'STX')
    print >> fh, ''
    print >> fh, 'END'

def viparr_params_to_amberfrcmod(viparr_dir, out_amberfrcmod):
    viparr_json_dict = load_viparr_dir(viparr_dir)
    fh = open(out_amberfrcmod, 'w')
    print >> fh, 'DES-AMBER'

    # Process mass
    print >> fh, 'MASS'
    for item in viparr_json_dict['mass']:
	if not isinstance(item['type'], list):
	    item['type'] = item['type'].split(' ')
        print >> fh, "%-2s %5s%14.3f               %s" \
                     % (item['type'][0], item['params']['amu'],
                        0.616, item['memo'])

    #Process bond
    print >> fh, ''
    print >> fh, 'BOND'
    for stretch_term in viparr_json_dict['stretch_harm']:
        k, length = stretch_term['params']['fc'], stretch_term['params']['r0']
	if not isinstance(stretch_term['type'], list):
		stretch_term['type'] = stretch_term['type'].split(' ')
        type1, type2 = stretch_term['type']
        print >>fh, "%-2s-%-2s%7.1f%10.3f" % (type1, type2, k, length)

    #Process angle
    print >> fh, ''
    print >> fh, 'ANGL'
    for angle_term in viparr_json_dict['angle_harm']:
        k, theta = angle_term['params']['fc'], angle_term['params']['theta0']
	if not isinstance(angle_term['type'], list):
		angle_term['type'] = angle_term['type'].split(' ')
        type1, type2, type3 = angle_term['type']
        print >>fh, "%-2s-%-2s-%-2s%8.1f%12.2f" % (type1, type2, type3, k, theta)

    #Process proper torsion
    print >> fh, ''
    print >> fh, 'DIHE'

    amber_torsion_list = []
    for proper_term in viparr_json_dict['dihedral_trig']:
        # Clean up
	if not isinstance(proper_term['type'], list):
		proper_term['type'] = proper_term['type'].split(' ')
        for n, t in enumerate(proper_term['type']):
            if t == '*':
                proper_term['type'][n] = 'X'

        type1, type2, type3, type4 = proper_term['type']
        phase_angle = proper_term['params']['phi0']

        #-----------------------------------------------#
        # keys are phi0, fc0, fc1, ... fc6. Ignore fc0. #
        #-----------------------------------------------#
        for key, value in proper_term['params'].iteritems():
            if '0' in key : 
                continue  
            elif value == 0:
                continue
            elif value < 0:
                periodicity = int(key[-1])
                k = value
                #amber_torsion_list.append([type1, type2, type3, type4, 1, -k, phase_angle+180, periodicity])
                amber_torsion_list.append([type1, type2, type3, type4, 1, k, phase_angle, periodicity])
            else:
                periodicity = int(key[-1])
                k = value
                amber_torsion_list.append([type1, type2, type3, type4, 1, k, phase_angle, periodicity])

        for m, item in enumerate(amber_torsion_list[:-1]):
            type1, type2, type3, type4, _, k, phase_angle, periodicity = item
            if item[:4] == amber_torsion_list[m+1][:4]:
                print >> fh, "%-2s-%-2s-%-2s-%-2s%4d%9.3f%12.1f%14d." \
                             % (type1, type2, type3, type4, 1, k, phase_angle, periodicity*-1)
            else:
                print >> fh, "%-2s-%-2s-%-2s-%-2s%4d%9.3f%12.1f%14d." \
                             % (type1, type2, type3, type4, 1, k, phase_angle, periodicity)
        print >> fh, "%-2s-%-2s-%-2s-%-2s%4d%9.3f%12.1f%14d." % tuple(amber_torsion_list[-1])


    #Process improper torsion
    print >> fh, ''
    print >> fh, 'IMPR'

    for improper_term in viparr_json_dict['improper_trig']:
        # Clean up
	if not isinstance(improper_term['type'], list):
		improper_term['type'] = improper_term['type'].split(' ')
        for n, t in enumerate(improper_term['type']):
            if t == '*':
                improper_term['type'][n] = 'X'

        type1, type2, type3, type4 = improper_term['type']
        phase_angle = improper_term['params']['phi0']

        #-----------------------------------------------#
        # keys are phi0, fc0, fc1, ... fc6. Ignore fc0. #
        #-----------------------------------------------#
        for key, value in improper_term['params'].iteritems():
            if '0' in key : 
                continue  
            elif value == 0:
                continue
            elif value < 0:
                periodicity = int(key[-1])
                k = value
                print >> fh, "%-2s-%-2s-%-2s-%-2s         %-4.1f%12d.%10d." \
                             % (type1, type2, type3, type4, k, phase_angle, periodicity)
                             #% (type1, type2, type3, type4, -k, phase_angle+180, periodicity)
            else:
                periodicity = int(key[-1])
                k = value
                print >> fh, "%-2s-%-2s-%-2s-%-2s         %-4.1f%12d.%10d." \
                             % (type1, type2, type3, type4, k, phase_angle, periodicity)
    # Process vdw
    print >> fh, ''
    print >> fh, 'NONB'
    for vdw_term in viparr_json_dict['vdw1']:
        # Clean up
        if not isinstance(vdw_term['type'], list):
            vdw_term['type'] = vdw_term['type'].split(' ')
        atom_type = vdw_term['type'][0]
        sigma, epsilon = vdw_term['params']['sigma'], vdw_term['params']['epsilon']
        Rmin_2 = sigma * 2**(1./6) / 2.
        print >> fh, "%4s%16.4f%8.4f" % (atom_type, Rmin_2, epsilon)
if __name__ == '__main__':
    usage = '%prog [opts] \n'
    arg = argparse.ArgumentParser( usage )
    arg.add_argument('--viparr_dir', default='viparr_dirs/stx-amber_protein' )
    arg.add_argument('--outdat', default='test/stx-amber-proteins-new.dat')
    args = arg.parse_args()
    viparr_params_to_amberdat(args.viparr_dir, args.outdat)
    #viparr_template_to_amberlib(args.viparr_dir, 'STx_amber_ffs/amino_stxamber.lib', 'test/amino_stxamber.lib')
    #viparr_template_to_amberlib(args.viparr_dir, 'STx_amber_ffs/aminont_stxamber.lib', 'test/aminont_stxamber.lib')
    #viparr_template_to_amberlib(args.viparr_dir, 'STx_amber_ffs/aminoct_stxamber.lib', 'test/aminoct_stxamber.lib')
