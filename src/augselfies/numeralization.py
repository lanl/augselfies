import re
import copy
import logging 
from typing import Union 
import selfies as sf 
import selfies.constants
import selfies.grammar_rules
import group_selfies.grammar_rules
import group_selfies.constants
import augselfies.constants
from augselfies.constants import INT_ALPHABET

augselfies_logger = logging.getLogger(__name__)

'''
Routines for handing numeric data, particularly ring and branch indices, in SELFIES and numSELFIES. 
'''

def is_branch(token_str:str)->bool:
    '''Returns true if token is a branch token, false otherwise'''
    return "branch" in token_str.lower()


def is_ring(token_str:str)->bool:
    '''Returns true if token is a ring token, false otherwise'''
    return "ring" in token_str.lower()

def is_group(token_str:str,group_regex = re.compile("\[\:[\s\S]?\d+.*\]"))->bool:
    if group_regex.match(token_str):
        return True
    else:
        return False 



def int_to_str(val:int,base_int:16):
    if base_int == 10:
        return str(val)
    elif base_int == 16:
        if val >= 0:
            return hex(val)[2:]
        else:
            hex_str =  hex(val)
            return hex_str[0]+hex_str[3:]
    else:
        if base_int > len(INT_ALPHABET):
            raise ValueError(f"Base {base_int} too large for alphabet {INT_ALPHABET}")
        
        str_digits = []

        rem_val = int(abs(val)) 
        power = 0
        while rem_val > 0:
            remainder = rem_val%(base_int)
            str_digits.append(INT_ALPHABET[remainder])
            rem_val = int(rem_val/base_int)
        str_digits.reverse()
        if val < 0:
            str_digits.insert(0,"-")
        return "".join(str_digits)
    
def natural_number_to_tokens(natural_number:int,single_token = True,base_int: int = 16)->list[str]:
    '''Converts integer to list of =tokens [0] to [f]'''
    if natural_number < 0 :
        raise ValueError(f"{natural_number} is negative and not a natural number")
    token_str = int_to_str(natural_number,base_int=base_int)
    tokens = []
    if (single_token):
        return [f"[{token_str}]"]
    else:
        for char in token_str:
            tokens.append(f"[{char}]")
        return tokens 


def tokens_to_natural_number(numeric_tokens:list[str],base_int:int = 16)->int:
    numeric_str = "".join(numeric_tokens).replace("[","").replace("]","")
    return int(numeric_str,base_int)




def get_token_list(selfies_like_string:str)->list[str]:
    '''Returns list of tokens for string tokenized by [a][b][c]'''
    if not selfies_like_string:
        token_list = []
    else:    
        token_list = selfies_like_string.split("][")
        token_list[0] = token_list[0][1:]
        token_list[-1] = token_list[-1][:-1]
        token_list = list(map(lambda x: f"[{x}]",token_list))
    return token_list

def get_n_sf_indices_from_token(indexing_regex: re.Pattern,token: str)->int:
    '''Returns the number of SELFIES index tokens following a graph (e.g., branch, ring) token.'''
    n_index_tokens = indexing_regex.match(token.lower()).groups()[0]
    if n_index_tokens.isdigit():
        n_index_tokens = int(n_index_tokens)
        n_index_tokens = max(n_index_tokens,1)
    else:
        n_index_tokens = 1
    return n_index_tokens

def get_n_ring_or_branch_indices(token:str,branch_regex=re.compile(r".*branch(\d*).*"),ring_regex=re.compile(r".*ring(\d*).*")):
    '''Returns the number of indexing tokens following a token. 0 if not a branch or ring'''
    n_index_tokens = 0
    if is_branch(token):
        n_index_tokens = get_n_sf_indices_from_token(branch_regex,token)
    elif is_ring(token):
        n_index_tokens = get_n_sf_indices_from_token(ring_regex,token)
    return n_index_tokens

def get_n_ring_or_group_indices(token:str,ring_regex=re.compile(r".*ring(\d*).*")):
    '''Returns the number of indexing tokens following a token. 0 if not a branch or ring'''
    n_index_tokens = 0
    if is_ring(token):
        n_index_tokens = get_n_sf_indices_from_token(ring_regex,token)
    if is_group(token):
        n_index_tokens = 1
    return n_index_tokens

SELFIES_TOKEN_DIGIT_REGEX = re.compile(r"\[(\D*)\d*\]")

def selfies_to_num_selfies(selfies_string:str,base_int:int = 16,is_gselfies = False,single_numeric_token:bool = True)->str:
    '''
    Converts a SELFIES string to a numSELFIES string. numSELFIES are essentially the same as SELFIES,
    but they cannot be directly parsed by the SELFIES grammar
    initial selfies must have each character surrounded byur_ square brackets
    Args:
        selfies_string (str): SELFIES string
    Returns:
        str: numSELFIES string
    '''
    if is_gselfies == True:
        raise NotImplementedError
    # get each elemental token, remove remaining brackets, put them back 
    try:
        token_list = get_token_list(selfies_string)
    except:
        token_list = []
    num_selfies_token_list = [] 
    idx = 0
    token_no_digit_regex = SELFIES_TOKEN_DIGIT_REGEX
    try:
        while idx < len(token_list):
            num_selfies_tokens = []
            cur_token = token_list[idx]
            n_index_tokens = get_n_ring_or_branch_indices(cur_token)
            if n_index_tokens >= 1:
                #NOTE: this uses the indexing code that the `selfies` installation uses or `group_selfies` for is_gselfies = True
                selfies_args = token_list[(idx+1):(idx+1+n_index_tokens)]
                index_value = selfies.grammar_rules.get_index_from_selfies(*selfies_args)
                base_token = token_no_digit_regex.match(cur_token).groups()[0]
                numeric_tokens = natural_number_to_tokens(index_value,single_token=single_numeric_token,base_int=base_int)
                if single_numeric_token:
                    num_selfies_tokens.append(f"[{base_token}]")
                else:
                    num_selfies_tokens.append(f"[{base_token}{len(numeric_tokens)}]")
                num_selfies_tokens.extend(numeric_tokens)
                idx+=n_index_tokens+1
            else:
                num_selfies_tokens = [f"{cur_token}"]
                idx+=1
            num_selfies_token_list.extend(num_selfies_tokens)
    except Exception as err: #will return empty string
        augselfies_logger.log(logging.ERROR,f"{err}")
        num_selfies_token_list = []
    return "".join(num_selfies_token_list)

def num_selfies_to_selfies(num_selfies_string:str,is_gselfies:bool = False,base_int:int = 16)->str:
    '''
    Converts a numSELFIES string back to a SELFIES string. Information should be preserved roundtrip.
    Uses indexing alphabet as defined in SELFIES. 
    NOTE: Do not use this for NumGroupSELFIES
    Args:
        num_selfies_string (str): numSELFIES string
    Returns:
        str: SELFIES string
    '''
    if is_gselfies == True:
        raise NotImplementedError
    try:
        token_list = get_token_list(num_selfies_string)
    except:
        token_list = []
    selfies_token_list = [] 
    idx = 0
    token_no_digit_regex = SELFIES_TOKEN_DIGIT_REGEX
    try:
        while idx < len(token_list):
            selfies_tokens = []
            cur_token = token_list[idx]
            n_index_tokens = get_n_ring_or_branch_indices(cur_token)
            if n_index_tokens >= 1:
                index_value = tokens_to_natural_number(token_list[(idx+1):(idx+1+n_index_tokens)],base_int=base_int)
                new_index_tokens = selfies.grammar_rules.get_selfies_from_index(index_value)
                base_token = token_no_digit_regex.match(cur_token).groups()[0]
                selfies_tokens.append(f"[{base_token}{len(new_index_tokens)}]")
                idx+=n_index_tokens+1
                selfies_tokens.extend(new_index_tokens)  
            else:
                selfies_tokens = [f"{cur_token}"]
                idx+=1
            selfies_token_list.extend(selfies_tokens)
    except Exception as err: #will return empty string
        augselfies_logger.log(logging.ERROR,f"{err}")
        selfies_token_list = ""
    return "".join(selfies_token_list)

def num_group_selfies_to_group_selfies(num_group_selfies:str,base_int = 16,
                              index_alphabet_to_return= augselfies.constants.DEFAULT_SELFIES_INDEX_ALPHABET)->str:
    '''
    NOTE: The reverse transformation from GroupSELFIES to NumGroupSELFIES does not exist. Rather, one should create GroupSELFIES directly from SMILES
    but with an amended index alphabet.
    Args:
        num_selfies (str):
        base_int (int, optional):  Defaults to 16.
        index_alphabet_to_return (_type_, optional): Defaults to augselfies.constants.DEFAULT_SELFIES_INDEX_ALPHABET.

    Raises:
        ValueError:

    Returns:
        str:
    '''
    if base_int > len(INT_ALPHABET):
        raise ValueError(f"Base {base_int} too large for alphabet {INT_ALPHABET} of size {len(INT_ALPHABET)}")
    int_token_regex = re.compile(f"\[([{INT_ALPHABET[:base_int]}]+\])")
    token_list = get_token_list(num_group_selfies)
    return_token_list = []
    try:
        for i, token in enumerate(token_list):
            int_match = int_token_regex.match(token)
            if int_match is None:
                return_token_list.append(token)
            else:
                num_str = int_match.groups()[0]
                num_value = tokens_to_natural_number(num_str,base_int=base_int)
                new_tokens = group_selfies.grammar_rules.get_selfies_from_index(num_value,index_alphabet_to_return)
                return_token_list.extend(new_tokens)
    except Exception as err: #will return empty string
        augselfies_logger.log(logging.ERROR,f"{err}")
        return_token_list = [] 
    return "".join(return_token_list) 

def smiles_to_num_selfies(smiles_string:str)->str:
    '''Wrapper for converting SMILES directly to numSELFIES'''
    return selfies_to_num_selfies(sf.encoder(smiles_string))

def num_selfies_to_smiles(num_selfies_string:str)->str:
    '''Wrapper for convering numSELFIES directly to SMILES'''
    return sf.decoder(num_selfies_to_selfies(num_selfies_string))

def get_tokens_handler(token_str:Union[str,list[str]])->list[str]:
    if isinstance(token_str,str):
        token_list = get_token_list(token_str)
    else:
        token_list = token_str
    return token_list

def num_selfies_to_text_num_encodings(num_selfies:Union[str,list[str]],num_token:str= "[NUM]",base_int:int = 16,omit_zero_token:bool = False,merge_numeric_tokens:bool = False)->tuple[list[str],list[float]]:
    '''
    Separates out numSELFIES into a pair of text (with "[NUM]" tokens) and continuous numerical encodings.
    Note that only nonnegative integer hex strings are supported as those are the numerical values present in numSELFIES
    as indexing tokens 
    For an overview of continuous numerical embeddings, see https://doi.org/10.48550/arXiv.2310.02989
    '''
    natural_int_regex = re.compile(r"\[([\da-f]+)\]")
    default_value = 1
    if merge_numeric_tokens:
        num_selfies = merge_numeric_token_integers(num_selfies,numeric_token_regex=natural_int_regex,base_int=base_int)
    if isinstance(num_selfies,str):
        token_list = get_token_list(num_selfies)
    elif isinstance(num_selfies,list):
        token_list = num_selfies
    else:
        raise ValueError(f"Invalid type of numselfies: {type(num_selfies)}")
    num_matches = list(map(lambda x: natural_int_regex.match(x),token_list))
    num_values = list(map(lambda x: default_value if x is None else int(x.groups()[0],base_int),num_matches))
    new_tokens = list(map(lambda x,y: x if y is None else num_token,token_list,num_matches))
    if omit_zero_token:
        zero_indices = [idx for idx, token in enumerate(token_list) if token == "[0]"]
        for index in zero_indices:
            num_values[index] = default_value
            new_tokens[index] = token_list[index]        
    return new_tokens, num_values



def merge_numeric_token_integers(token_str:Union[str,list[str]],numeric_token_regex = re.compile(r"\[([\da-f]+)\]"),base_int:int = 16):
    token_list = get_tokens_handler(token_str)
    token_iter = iter(token_list)
    return_tokens = []
    while True:
        try:
            cur_token = next(token_iter)
            if re.match(numeric_token_regex,cur_token):
                num_tokens = [cur_token]
                while True:
                    try:
                        next_token = next(token_iter)
                        if re.match(numeric_token_regex,next_token):
                            num_tokens.append(next_token)
                        else:
                            break 
                    except StopIteration:
                        next_token = ""
                        break
                int_val = 0
                for i,num_token in enumerate(num_tokens):
                    int_index = int(numeric_token_regex.match(num_token).groups()[0],base= base_int)
                    if int_index >= base_int:
                        raise ValueError(f"Numeric tokens {num_tokens} are too large to be converted to a single number for base {base_int}")
                    int_val+=int_index*base_int**(len(num_tokens)-i-1)
                combined_num_token = f"[{int_to_str(int_val,base_int=base_int)}]"
                return_tokens.extend([combined_num_token,next_token])
            else:
                return_tokens.append(cur_token)
        except StopIteration:
            break  
    return "".join(return_tokens)

HEX_INDEX_ALPHABET = augselfies.constants.DEFAULT_NUM_SELFIES_INDEX_ALPHABET
