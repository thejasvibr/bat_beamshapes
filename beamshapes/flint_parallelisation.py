'''
Flint parallelisation
=====================
Holds convenience functions which help in parallelising FLINT
code. 

The idea is to:
    * convert acb objects into strings to allow pickling
    * convert string --> acb to allow computation in function
    * convert acb back to string and export results

'''

import flint 
acb = flint.acb


def conv_numeric_to_str(entry):
    '''
    

    '''
    try:
        str_out = conv_acb_to_str(entry)
    except:
        str_out = str(entry)
    return str_out

def conv_acb_to_str(entry):
    return entry.mid().str()

def conv_str_acb(entry):
    if not isinstance(entry, str):
        raise TypeError(f'Input  {entry} is not string type.')
    else:
        # split into two parts 
        real, imag = entry.split('] + [')
        real = real.replace('[','')
        for character in [']','j']:
            imag = imag.replace(character, '')

    return acb(real, imag)
        

def interchange_params_and_str(params, to_str=True):
    '''
    Converts the midpoint of the entries into strings.
    '''
    if to_str:
        return  { key: conv_numeric_to_str(entry) for key, entry in params.items()}
    else:
        return  { key: acb(entry) for key, entry in params.items()}
