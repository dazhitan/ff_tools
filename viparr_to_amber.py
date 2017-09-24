#!/usr/bin/env python
'''
@author: Dazhi Tan
@created: July 14th, 2017
'''
import os, re, json, argparse
import amberlib
import amberdat
import desmondff
from amberlib import AmberLib
from amberdat import AmberDat
from desmondff import DesmondFF

def desmond_template_to_amberlib(viparr_dir, in_amberlib, out_amberlib):

    desmond_ff = DesmondFF(viparr_dir)
    template_amberlib = AmberLib(in_amberlib)

    #-------------------------#
    # Process templates files #
    #-------------------------#
    template_jsons = []
    for key, json in desmond_ff.iteritems():
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

                #atom_name = '"%s"' % atom_name
                #bonded_type = '"%s"' % bonded_type
                #nonbonded_type = '"%s"' % nonbonded_type
                #charge = '%.6f' % charge

                # Set properties
                template_amberlib.setAtomType(resiname, atom_name, bonded_type) 
                template_amberlib.setAtomCharge(resiname, atom_name, charge) 

    template_amberlib.printLib(out_amberlib)


def desmond_params_to_amberdat(viparr_dir, out_amberdat):

    desmond_ff = DesmondFF(viparr_dir)
    amber_dat = AmberDat()

    amber_dat.setTitle('STX-AMBER')

    # Process mass
    for massparam in desmond_ff['mass']:
        # Clean up
        for n, t in enumerate(massparam['type']):
            if t == '*':
                massparam['type'][n] = 'X'
            else:
                t = t.replace('&', '*')
                massparam['type'][n] = t
        atomtype = massparam['type'][0]
        mass = float(massparam['params']['amu'])
        amber_dat.addMass(atomtype, mass)

    #Process stretch
    for stretchparam in desmond_ff['stretch_harm']:
        # Clean up
        for n, t in enumerate(stretchparam['type']):
            if t == '*':
                stretchparam['type'][n] = 'X'
            else:
                t = t.replace('&', '*')
                stretchparam['type'][n] = t
        atomtype1, atomtype2 = stretchparam['type']
        fc = stretchparam['params']['fc']
        r0 = stretchparam['params']['r0']
        amber_dat.addStretch(atomtype1, atomtype2, fc, r0)

    #Process angle
    for angleparam in desmond_ff['angle_harm']:
        # Clean up
        for n, t in enumerate(angleparam['type']):
            if t == '*':
                angleparam['type'][n] = 'X'
            else:
                t = t.replace('&', '*')
                angleparam['type'][n] = t
        atomtype1, atomtype2, atomtype3, = angleparam['type']
        fc = angleparam['params']['fc']
        theta0 = angleparam['params']['theta0']
        amber_dat.addAngle(atomtype1, atomtype2, atomtype3, fc, theta0)

    #Process proper torsion

    amber_torsion_list = []
    for properparam in desmond_ff['dihedral_trig']:
        # Clean up
        for n, t in enumerate(properparam['type']):
            if t == '*':
                properparam['type'][n] = 'X'
            else:
                t = t.replace('&', '*')
                properparam['type'][n] = t

        atomtype1, atomtype2, atomtype3, atomtype4, = properparam['type']
        phase_angle = properparam['params']['phi0']

        #-----------------------------------------------#
        # keys are phi0, fc0, fc1, ... fc6. Ignore fc0. #
        #-----------------------------------------------#
        for key, value in properparam['params'].iteritems():
            if '0' in key : 
                continue  
            else:
                periodicity = int(key[-1])
                fc = value
                amber_torsion_list.append(
                    [atomtype1, atomtype2, atomtype3, atomtype4, 1, fc, phase_angle, periodicity])
    for torsion_term in amber_torsion_list:
        amber_dat.addProper(*torsion_term)

    #Process improper torsion
    for improperparam in desmond_ff['improper_trig']:
        # Clean up
        for n, t in enumerate(improperparam['type']):
            if t == '*':
                improperparam['type'][n] = 'X'
            else:
                t = t.replace('&', '*')
                improperparam['type'][n] = t

        atomtype1, atomtype2, atomtype3, atomtype4 = improperparam['type']
        phase_angle = improperparam['params']['phi0']

        #-----------------------------------------------#
        # keys are phi0, fc0, fc1, ... fc6. Ignore fc0. #
        #-----------------------------------------------#
        for key, value in improperparam['params'].iteritems():
            if '0' in key : 
                continue  
            elif value == 0:
                continue
            elif value < 0:
                periodicity = int(key[-1])
                fc = value
                amber_dat.addImproper(atomtype1, atomtype2, atomtype3, atomtype4, -fc, phase_angle+180, periodicity)
            else:
                periodicity = int(key[-1])
                fc = value
                amber_dat.addImproper(atomtype1, atomtype2, atomtype3, atomtype4, fc, phase_angle, periodicity)
    # Process vdw
    for vdwparam in desmond_ff['vdw1']:
        atomtype = vdwparam['type'][0].replace('&', '*')
        sigma, epsilon = vdwparam['params']['sigma'], vdwparam['params']['epsilon']
        Rmin_2 = sigma * 2**(1./6) / 2.
        amber_dat.addVdw(atomtype, Rmin_2, epsilon)

    amber_dat.printDat(out_amberdat)

if __name__ == '__main__':
    usage = '%prog [opts] \n'
    arg = argparse.ArgumentParser( usage )
    arg.add_argument('--viparr_dir', default='stx_amber_protein_viparr' )
    arg.add_argument('--outdat', default='stx_amber_protein_new.dat')
    args = arg.parse_args()
    desmond_params_to_amberdat(args.viparr_dir, args.outdat)
    #desmond_template_to_amberlib(args.viparr_dir, 'allamino_stxamber.lib', 'allamino_stxamber_new.lib')
    #viparr_template_to_amberlib(args.viparr_dir, 'STx_amber_ffs/aminont_stxamber.lib', 'test/aminont_stxamber.lib')
    #viparr_template_to_amberlib(args.viparr_dir, 'STx_amber_ffs/aminoct_stxamber.lib', 'test/aminoct_stxamber.lib')
