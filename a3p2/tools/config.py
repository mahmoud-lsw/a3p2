#===========================================================================
# Imports

import logging
import pyfits

#===========================================================================
# Functions & classes

#---------------------------------------------------------------------------
class ConfigBase(object) :
    def __init__(self, defaults=None) :
        self.__dict__.update(**defaults)
    def check_dict(self, dict_) :
        set_default = set(self.__dict__)
        set_dict = set(dict_)
        return (set_default.intersection(set_dict), set_dict - set_default)
    def init_from_dict(self, dict_, verbose=True) :
        update_keys, unrecognized_keys = self.check_dict(dict_)
        if verbose :
            if len(unrecognized_keys) :
                logging.warning('Unrecognized keys: {0}'.format(str(list(unrecognized_keys))))
            if len(update_keys) :
                logging.info('Updating entries for keys: {0}'.format(str(list(update_keys))))
        for key in update_keys :
            self.__dict__[key] = dict_[key]
    def print_(self) :
        maxkeylen = str(max([len(k) for k in self.__dict__.keys()]))
        for key, item in sorted(self.__dict__.iteritems(), key=lambda (k, v) : k) :
            logging.info(('{0:' + maxkeylen + '} : {1}').format(key, item))
    def write_to_fits_header(self, fh, exclude=[], rename={}) :
        for key, item in sorted(self.__dict__.iteritems(), key=lambda (k, v) : k) :
            if key in exclude :
                continue
            if key in rename.keys() :
                key = rename[key]
            if item == None :
                item = ''
            try :
                len(item)
                item = str(item)
            except :
                pass
            if len(key) > 8 :
                logging.warning('Shorten key {0} to {1} for FITS header.'.format(key, key[:8].upper()))
            fh.update(key[:8].upper(), item)

#===========================================================================
#===========================================================================
