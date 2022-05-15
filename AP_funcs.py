#!/usr/bin/env python
# coding: utf-8
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Microsoft VS header
#--------------------------------------------------#
import os 
import sys
import os.path
from sys import platform
from pathlib import Path
#--------------------------------------------------#
#if os.name == 'nt' or platform == 'win32':
#    print("Running on Windows")
#    if 'ptvsd' in sys.modules:
#        print("Running in Visual Studio")
#        try:
#            os.chdir(os.path.dirname(__file__))
#            print('CurrentDir: ', os.getcwd())
#        except:
#            pass
##--------------------------------------------------#
#    else:
#        print("Running outside Visual Studio")
#        try:
#            if not 'workbookDir' in globals():
#                workbookDir = os.getcwd()
#                print('workbookDir: ' + workbookDir)
#                os.chdir(workbookDir)
#        except:
#            pass
#--------------------------------------------------#
if os.name != 'nt' and platform != 'win32':
    print("Not Running on Windows")
#--------------------------------------------------#
from rdkit import Chem
from rdkit.Chem import AllChem
#--------------------------------------------------#
import random
import itertools
from pathlib import Path

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def replace_n(str1, n, str2):
    letters = (
        str2 if i == n else char
        for i, char in enumerate(str1)
    )
    return ''.join(letters)
# ============================================================================================================================ #
def get_num_combinations(list_a):
    count_x=1
    for i in range(len(list_a)):
        count_x=count_x*len(list_a[i])
    return count_x
# ============================================================================================================================ #
def get_ith_combination(list_a, num_i): #list_a is a list of lists to perform combination
    ith_combination=[]
    length_list=[]
    mod_list=[]

    ith_list=[]
    for one_list in list_a:
        length_list.append(len(one_list))
        mod_list.append(1)
        ith_list.append(0)

    for i in range(len(length_list)):
        for j in range(len(length_list)-i):
            mod_list[i]=length_list[-j-1]*mod_list[i]
    for i in range(len(length_list)-1):
        ith_list[i]= int((num_i % mod_list[i]) / mod_list[i+1]) + 1
        if (num_i % mod_list[i]) % mod_list[i+1]==0:
            ith_list[i]=ith_list[i]-1
        if num_i % mod_list[i]==0:
            ith_list[i]=length_list[i]
    ith_list[-1]=num_i % mod_list[-1]
    if ith_list[-1]==0:
        ith_list[-1]=length_list[-1]
    for j in range(len(ith_list)):
        ith_combination.append(list_a[j][ith_list[j]-1])
    return ith_combination
# ============================================================================================================================ #
def cart_prod( alistoflists): # return combinations of elements, each from a different set.
    # [[a,b,c],[d,e,f],[g,h], ...] converts to
    # [(a,d,g),(b,d,g),(c,d,g),(a,d,h),(b,d,h), ...]
    prod_list=alistoflists[0]
    if len(alistoflists)>1:
        for i in range(len(alistoflists)-1):
            prod_list=list(itertools.product(prod_list,alistoflists[i+1]))
            if i==0:
                pass
            else:
                set_list=prod_list
                prod_list=[]
                for j in range(len(set_list)):
                    inbracket=list(set_list[j][0])
                    nobracket=set_list[j][1]
                    inbracket.append(nobracket)
                    prod_list.append(tuple(inbracket))
    else:
        temp_list=[]
        for i in range(len(prod_list)):
            temp_list.append((prod_list[i],))
        prod_list=temp_list
    return prod_list
# ============================================================================================================================ #
def randomList(a): 
    b = [] 
    for i in range(len(a)): 
        element = random.choice(a) 
        a.remove(element) 
        b.append(element) 
    return b
# ============================================================================================================================ #
def iftuplestrinlist(a,b): # Check if all strs in one tuple are contained in a list
    allstrsin=False
    countall=0
    for i in range(len(a)):
        if a[i] in b:
            countall+=1
    if countall==len(a):
        allstrsin=True
    return allstrsin
# ============================================================================================================================ #
def find_nth(string_a, string_b, n): # Find the nth occurence of string_b in string_a (, find_nth('aaaaa','aa',3) = -1).
    s = string_a.find(string_b)
    while s >= 0 and n > 1:
        s = string_a.find(string_b, s+len(string_b))
        n -= 1
    return s
# ============================================================================================================================ #
def simplify_enzymeid(enzyme_a):
    # For better texture output
    if enzyme_a=='ecKAYLA_backward' or enzyme_a=='ecKAYLA_forward':
        return enzyme_a
    else:
        enzyme_b=enzyme_a
        thethirddotindex=find_nth(enzyme_b,'.',3)
        theshortlineindex=enzyme_b.index('_')
        start=thethirddotindex
        while start<theshortlineindex:
            enzyme_b=replace_n(enzyme_b,thethirddotindex,'')
            start+=1
        return enzyme_b
# ============================================================================================================================ #
def remove_hydrogen_nodes_in_rule(one_ruleH):
    while (one_ruleH.find('\\')!=-1):
        location=one_ruleH.index('\\')
        one_ruleH=replace_n(one_ruleH,location,'') # Remove '\'
    while (one_ruleH.find('/')!=-1):
        location=one_ruleH.index('/')
        one_ruleH=replace_n(one_ruleH,location,'') # Remove '/'
    while (one_ruleH.find('[H]')!=-1):
        location=one_ruleH.index('[H]')
        one_ruleH=replace_n(one_ruleH,location,'') # Remove '['
        one_ruleH=replace_n(one_ruleH,location,'') # Remove 'H'
        one_ruleH=replace_n(one_ruleH,location,'') # Remove ']'
    while (one_ruleH.find('[H:')!=-1):
        location=one_ruleH.index('[H:')
        if one_ruleH[location+4]=="]":
            one_ruleH=replace_n(one_ruleH,location,'') # Remove '['
            one_ruleH=replace_n(one_ruleH,location,'') # Remove 'H'
            one_ruleH=replace_n(one_ruleH,location,'') # Remove ':'
            one_ruleH=replace_n(one_ruleH,location,'') # Remove '*'
            one_ruleH=replace_n(one_ruleH,location,'') # Remove ']'
        elif one_ruleH[location+5]=="]":
            one_ruleH=replace_n(one_ruleH,location,'') # Remove '['
            one_ruleH=replace_n(one_ruleH,location,'') # Remove 'H'
            one_ruleH=replace_n(one_ruleH,location,'') # Remove ':'
            one_ruleH=replace_n(one_ruleH,location,'') # Remove '*'
            one_ruleH=replace_n(one_ruleH,location,'') # Remove '*'
            one_ruleH=replace_n(one_ruleH,location,'') # Remove ']'
        else:
            print ("You just find [H:***]? Is that even possible?")
    while (one_ruleH.find('()')!=-1):
        location=one_ruleH.index('()')
        one_ruleH=replace_n(one_ruleH,location,'') # Remove '('
        one_ruleH=replace_n(one_ruleH,location,'') # Remove ')'
    return one_ruleH
# ============================================================================================================================ #
def z_mkdir(new_dir):
    if os.path.isdir(new_dir):
        pass
    else:
        head, tail = os.path.split(new_dir)
        if head and not os.path.isdir(head):
            z_mkdir(head)
        if tail:
            os.mkdir(new_dir)
    return
# ============================================================================================================================ #
def multipletrialnames(result_name, num_digits=4):
    num_digits=int(num_digits)
    z_mkdir("../results")
    result_name = result_name + num_digits*"0"
    
    while os.path.isdir("../results/" + result_name)== True:
        result_name_list=[result_name[:-1*num_digits], result_name[-1*num_digits:]]
        result_name=result_name_list[0]+str(int('1'+result_name_list[1])+1)[-1*num_digits:]
    z_mkdir("../results/" + result_name)
    pathway_r_files = "pwy_r"
    z_mkdir("../results/" + result_name + "/" +  pathway_r_files)
    z_mkdir("../results/" + result_name + "/" +  pathway_r_files + "/pathways")
    return result_name

# ============================================================================================================================ #
from collections import OrderedDict,Counter
class OrderedCounter(Counter, OrderedDict):
     'Counter that remembers the order elements are first encountered'
     def __repr__(self):
         return '%s(%r)' % (self.__class__.__name__, OrderedDict(self))
     def __reduce__(self):
         return self.__class__, (OrderedDict(self),)

# ============================================================================================================================ #
class cartesian_product:
    def __init__(self, list_of_lists=[], unique_values=False):
        self._list_of_lists = list_of_lists
        self._unique_values = unique_values
        self._num_lists = len(self._list_of_lists)
        self._lengths = [len(x) for x in self._list_of_lists]
        self._total_size = 1
        for l in self._lengths:
            self._total_size *= l
        if (self._num_lists == 0 or self._total_size == 0):
            self._indices = None # this is like a flag indicating that the iteration is over
        else:
            self._indices = [0] * len(self._list_of_lists)
    def __iter__(self):
        return self
    def __len__(self):
        return self._total_size
    def get_current_item(self):
        if (self._indices == None):
            return []
        return [self._list_of_lists[i][self._indices[i]] for i in range(self._num_lists)]
    def increment(self):
        i = self._num_lists - 1
        self._indices[i] += 1
        while (i > 0 and self._indices[i] == self._lengths[i]):
            self._indices[i-1] += 1
            self._indices[i] = 0
            i -= 1
        if (self._indices[0] == self._lengths[0]):
            self._indices = None
    def __next__(self):
        if (self._unique_values):
            # if the current item contains repetitions, skip to the next one
            while (0 < len(set(self.get_current_item())) < self._num_lists):
                self.increment()
        if (self._indices == None):
            raise StopIteration()
        current_item = self.get_current_item()
        self.increment()
        return current_item

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
def test():
    print(cartesian_product([[0,1],[0,2],[0,3]], unique_values=False))
    for i in cartesian_product([[0,1],[0,2],[0,3]], unique_values=False):
        print(i)

# ============================================================================================================================ #
if (__name__ == '__main__'):
    test()
    multipletrialnames(result_name="1230", num_digits=5)