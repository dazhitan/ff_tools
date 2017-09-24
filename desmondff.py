#!/usr/bin/env python
'''
@author: Dazhi Tan
@created: September 22th, 2017
'''
import os, re, json
from glob import glob
from UserDict import IterableUserDict

class DesmondFF(IterableUserDict):
    """
    A class describing the Desmond force field format.
    """
    def __init__(self, ff_dir_path):
        IterableUserDict.__init__(self)
        self.ff_dir = ff_dir_path
        self.title = os.path.basename(self.ff_dir)
        self.loadDir()
    
    def loadDir(self):
        ff_file_list = glob(os.path.join(self.ff_dir, '*'))
        for file in ff_file_list:
            filename = os.path.basename(file)
            with open(file, 'r') as fh:
                json_object = json.load(fh)
                self.data[filename] = json_object
        self._cleantype()

    def _cleantype(self):
        for filename, json_object in self.data.iteritems():
            if isinstance(json_object, list):
                for elmnt in json_object:
                    if not isinstance(elmnt['type'], (tuple, list)):
                        elmnt['type'] = elmnt['type'].split(' ')
            else:
                continue

    def _find_atomtype_inconsistency(self):
        for filename in self.data.keys():
            if 'template' in filename:
                for resname, resparam in self.data[filename].iteritems():
                    for atom in resparam['atoms']:
                        if atom[-1][0] != atom[-1][1]:
                            print resname, atom[0], atom[-1][0], atom[-1][1]

"""
def lookup_lj_param(nonbonded_type, vdw_json):
    for item in vdw_json:
	if not isinstance(item['type'], list):
		item['type'] = item['type'].split(' ')
        if item['type'][0] == nonbonded_type:
            return item['params']['sigma']*0.1, \
                   item['params']['epsilon']*4.184

def lookup_mass(bonded_type, mass_json):
    for item in mass_json:
	if not isinstance(item['type'], list):
		item['type'] = item['type'].split(' ')
        if item['type'][0] == bonded_type:
            return item['params']['amu']

def lookup_element(atomic_number):
    element_dict = {1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C', 7:'N', 8:'O',
                    9:'F', 10:'Ne', 11:'Na', 12:'Mg', 13:'Al', 14:'Si', 15:'P',
                    16:'S', 17:'Cl', 18:'Ar', 19:'K', 20:'Ca', 25:'Mn', 
                    26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 34:'Se', 
                    35:'Br', 37:'Rb', 47:'Ag', 53:'I', 55:'Cs', 79:'Au', 80:'Hg'} 
    return element_dict[atomic_number]

def viparr_to_xml(viparr_dir, xml_file_name):

    viparr_json_dict = load_viparr_dir(viparr_dir)
    ff = ForceField()
    #-------------------------#
    # Process templates files #
    #-------------------------#
    template_jsons = []
    for key, json in viparr_json_dict.iteritems():
        if 'template' in key:
            template_jsons.append(json)
    for tj in template_jsons:
        for resiname, resiparam in tj.iteritems():
            ff.addResidue(resiname)
            for atomparam in resiparam['atoms']:
                (atom_name, atomic_number, charge, 
                 [bonded_type, nonbonded_type]) = atomparam

                sigma, epsilon = \
                lookup_lj_param(nonbonded_type, viparr_json_dict['vdw1'])

                mass = \
                lookup_mass(bonded_type, viparr_json_dict['mass'])

                element = lookup_element(atomic_number)

                #------------------------#
                # Clean up wierd symbols #
                #------------------------#
                atom_name = atom_name.replace("'", "pri")
                add_atom = ff.addAtom(resiname, atom_name)
                atom_type = add_atom.attrib['type']

                bonded_type = bonded_type.replace("&", "amp")
                bonded_type = bonded_type.replace("*", "")
                nonbonded_type = nonbonded_type.replace("&", "amp")

                ff.addAtom(resiname, atom_name)
                ff.addType(atom_type, bonded_type, element, mass) 
                ff.addNonbondedTerm(atom_type, charge, sigma, epsilon)

            if resiparam.has_key('bonds'): 
                for atom_name1, atom_name2 in resiparam['bonds']:
                    atom_name1 = atom_name1.replace("'", "pri")
                    atom_name2 = atom_name2.replace("'", "pri")
                    if '$' in atom_name1:
                        ff.addExternalBond(resiname, atom_name2)
                    elif '$' in atom_name2:
                        ff.addExternalBond(resiname, atom_name1)
                    else:
                        ff.addBond(resiname, atom_name1, atom_name2)

    #---------------#
    # Process rules #
    #---------------#
    es_scale, lj_scale = viparr_json_dict['rules']['es_scale'][-1], \
                         viparr_json_dict['rules']['lj_scale'][-1]
    ff.setNonbondedRules(es_scale, lj_scale)

    #----------------------#
    # Process stretch_harm #
    #----------------------#
    for stretch_term in viparr_json_dict['stretch_harm']:
        k, length = stretch_term['params']['fc']*4.184*10**2*2, \
                    stretch_term['params']['r0']*0.1
	if not isinstance(stretch_term['type'], list):
		stretch_term['type'] = stretch_term['type'].split(' ')
        bonded_type1, bonded_type2 = stretch_term['type']
        bonded_type1 = bonded_type1.replace("&", "amp")
        bonded_type1 = bonded_type1.replace("*", "")
        bonded_type2 = bonded_type2.replace("&", "amp")
        bonded_type2 = bonded_type2.replace("*", "")
        
        ff.addHarmonicBondForce(bonded_type1, bonded_type2, length, k) 

    #--------------------#
    # Process angle_harm # 
    #--------------------#
    for angle_term in viparr_json_dict['angle_harm']:
        k, theta = angle_term['params']['fc']*4.184*2, \
                   radians(angle_term['params']['theta0'])
	if not isinstance(angle_term['type'], list):
		angle_term['type'] = angle_term['type'].split(' ')
        bonded_type1, bonded_type2, bonded_type3 = angle_term['type']
        bonded_type1 = bonded_type1.replace("&", "amp")
        bonded_type1 = bonded_type1.replace("*", "")
        bonded_type2 = bonded_type2.replace("&", "amp")
        bonded_type2 = bonded_type2.replace("*", "")
        bonded_type3 = bonded_type3.replace("&", "amp")
        bonded_type3 = bonded_type3.replace("*", "")
        
        ff.addHarmonicAngleForce(bonded_type1, bonded_type2, bonded_type3, 
                                 theta, k) 

    #-----------------------#
    # Process dihedral_trig #
    #-----------------------#
    existing_atom_list = []
    initializer = 1
    for proper_term in viparr_json_dict['dihedral_trig']:
	if not isinstance(proper_term['type'], list):
		proper_term['type'] = proper_term['type'].split(' ')
        bonded_type1, bonded_type2, bonded_type3, bonded_type4 = proper_term['type']
        bonded_type1 = bonded_type1.replace("&", "amp")
        bonded_type1 = bonded_type1.replace("*", "")
        bonded_type2 = bonded_type2.replace("&", "amp")
        bonded_type2 = bonded_type2.replace("*", "")
        bonded_type3 = bonded_type3.replace("&", "amp")
        bonded_type3 = bonded_type3.replace("*", "")
        bonded_type4 = bonded_type4.replace("&", "amp")
        bonded_type4 = bonded_type4.replace("*", "")

        current_atom_list = [bonded_type1, bonded_type2, bonded_type3, bonded_type4]
        
        if current_atom_list != existing_atom_list:
            existing_atom_list = current_atom_list
            initializer = 1
        else:
            initializer += 1
        kwargs = dict()
        phase_angle = radians(proper_term['params']['phi0'])
        counter = initializer
        #-----------------------------------------------#
        # keys are phi0, fc0, fc1, ... fc6. Ignore fc0. #
        #-----------------------------------------------#
        for key, value in proper_term['params'].iteritems():
            if '0' in key : 
                continue  
            elif value > 0:
                periodicity = int(key[-1])
                kwargs.update({'periodicity%s' % counter:periodicity})
                kwargs.update({'phase%s' % counter:phase_angle})
                kwargs.update({'k%s' % counter:value*4.184})
                counter += 1
            #----------------------------------------------------------------#
            # Make sure all k values are positive. if fc is negative, add pi #
            # to the phase angle.                                            #
            #----------------------------------------------------------------#
            elif value < 0:
                periodicity = int(key[-1])
                kwargs.update({'periodicity%s' % counter:periodicity})
                kwargs.update({'phase%s' % counter:phase_angle+pi})
                kwargs.update({'k%s' % counter:value*4.184*-1})
                counter += 1
        #------------------------------------------------------------------#
        # First encouter of the 4-atom list. If counter == initializer, it #
        # means this torsion term is empty, and we ignore it.              #
        #------------------------------------------------------------------#
        if initializer == 1 and counter != initializer: 
            ff.addProperPeriodicTorsionForce(
            bonded_type1, bonded_type2, bonded_type3, bonded_type4, **kwargs)
        #---------------------------------------------------------------------#
        # If initializer > 1, it means that the same 4-atom list has multiple #
        # entries in the dihedral_trig file. So instead of adding a new term, #
        # we modify the existing term by adding more attributes.              #
        #---------------------------------------------------------------------#
        if initializer > 1 and counter != initializer:  
            ff.setProperPeriodicTorsionForce(
            bonded_type1, bonded_type2, bonded_type3, bonded_type4, **kwargs)

    #------------------------#
    # Process improper_trig  # 
    #------------------------#
    for improper_term in viparr_json_dict['improper_trig']:
	if not isinstance(improper_term['type'], list):
		improper_term['type'] = improper_term['type'].split(' ')
        bonded_type1, bonded_type2, bonded_type3, bonded_type4 = improper_term['type']
        bonded_type1 = bonded_type1.replace("&", "amp")
        bonded_type1 = bonded_type1.replace("*", "")
        bonded_type2 = bonded_type2.replace("&", "amp")
        bonded_type2 = bonded_type2.replace("*", "")
        bonded_type3 = bonded_type3.replace("&", "amp")
        bonded_type3 = bonded_type3.replace("*", "")
        bonded_type4 = bonded_type4.replace("&", "amp")
        bonded_type4 = bonded_type4.replace("*", "")

        kwargs = dict()
        phase_angle = radians(improper_term['params']['phi0'])
        counter = 1
        for key, value in improper_term['params'].iteritems():
            if '0' in key : 
                continue  
            elif value > 0:
                periodicity = int(key[-1])
                kwargs.update({'periodicity%s' % counter:periodicity})
                kwargs.update({'phase%s' % counter:phase_angle})
                kwargs.update({'k%s' % counter:value*4.184})
                counter += 1
            elif value < 0:
                periodicity = int(key[-1])
                kwargs.update({'periodicity%s' % counter:periodicity})
                kwargs.update({'phase%s' % counter:phase_angle+pi})
                kwargs.update({'k%s' % counter:value*4.184*-1})
                counter += 1
        if counter != 1:
            ff.addImproperPeriodicTorsionForce(
            bonded_type3, bonded_type1, bonded_type2, bonded_type4, **kwargs)
        
    #-------------------#
    # Write output file #
    #-------------------#
    doc = ElementTree(ff)
    doc.write(open(xml_file_name, 'w'), pretty_print=True)

if __name__ == '__main__':
    usage = '%prog [opts] \n'
    arg = argparse.ArgumentParser( usage )
    arg.add_argument('viparr_dir',  
                     help='Viparr directory')
    arg.add_argument('-o', '--output', 
                     help='File name of the output xml file')
    args = arg.parse_args()
    viparr_to_xml(args.viparr_dir, args.output)
"""
