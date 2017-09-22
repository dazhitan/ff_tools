#!/usr/bin/env python
'''
@author: Dazhi Tan
@created: July 14th, 2017
'''
import os, re, json, warnings
from lxml import etree
from lxml.etree import ElementBase
from lxml.etree import XPath
from glob import glob

def check_duplicates(parent, tag, **kwargs):
    """
    Check whether children with tag=tag and attributes=kwargs exist under the 
    parent. 
    """
    xpath_string = "./%s" % tag
    for key, value in kwargs.iteritems():
        xpath_string += "[@%s='%s']" % (key, value)
    if len(parent.xpath(xpath_string)) > 0:
        return parent.xpath(xpath_string)
    else:
        return False

def load_viparr_dir(viparr_dir):

    viparr_file_list = glob(os.path.join(viparr_dir, '*'))
    viparr_json_dict = dict()
    for vf in viparr_file_list:
        abspath = os.path.abspath(vf)
        basename = os.path.basename(abspath)
        json_object = json.load(open(abspath, 'r'))
        viparr_json_dict[basename] = json_object
    return viparr_json_dict
