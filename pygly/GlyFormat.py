from new_monosaccharide import Monosaccharide, Substituent
from Composition import Composition
import copy
from new_link import Linkage
from new_glycan import Glycan
from rec2 import _all_members as number_all_members
from rec2 import _all_links as number_all_links

class MonoFormatter:

    def addSubstituent(self, monosaccharide, sstr):
        raise "Incomplete implementation"
        s = Substituent()
        #set substituent number and linkage here
        s.set_substituent(sstr)
        
        new_monosaccharide = copy.copy(monosaccharide)
        new_monosaccharide.add_substituent(s)
        return new_monosaccharide

class GlycoCTMonoFormatter(MonoFormatter):

    _glycoCTAnomerMapping = {     Monosaccharide.alphaAnomer:   'a',
                                   Monosaccharide.betaAnomer:   'b',
                                Monosaccharide.unknownAnomer:   'x',
                             Monosaccharide.uncyclizedAnomer:   'o',
                             }

    _glycoCTConfigMapping = {         Monosaccharide.dConfig:   'd',
                                      Monosaccharide.lConfig:   'l',
                                Monosaccharide.unknownConfig:    '',
                                         None:    '',
                             }

    _glycoCTStemMapping =   {Monosaccharide.groStem: 'gro',
                             Monosaccharide.eryStem: 'ery',
                             Monosaccharide.ribStem: 'rib',
                             Monosaccharide.araStem: 'ara',
                             Monosaccharide.allStem: 'all',
                             Monosaccharide.altStem: 'alt',
                             Monosaccharide.glcStem: 'glc',
                             Monosaccharide.manStem: 'man',
                             Monosaccharide.treStem: 'tre',
                             Monosaccharide.xylStem: 'xyl',
                             Monosaccharide.lyxStem: 'lyx',
                             Monosaccharide.gulStem: 'gul',
                             Monosaccharide.idoStem: 'ido',
                             Monosaccharide.galStem: 'gal',
                             Monosaccharide.talStem: 'tal',
                             }

    _glycoCTClassMapping =  {Monosaccharide.classTRI:   'TRI',
                             Monosaccharide.classTETRA: 'TET',
                             Monosaccharide.classPENTA: 'PEN',
                             Monosaccharide.classHEXA:  'HEX',
                             Monosaccharide.classHEPTA: 'HEP',
                             Monosaccharide.classOCTA:  'OCT',
                             Monosaccharide.classNONA:  'NON',
                             }

    _glycoCTModMapping =   {Monosaccharide.dMod:                'd',
                            Monosaccharide.ketoMod:          'keto',
                            Monosaccharide.enMod:              'en',
                            Monosaccharide.aMod:                'a',
                            Monosaccharide.aldiMod:          'aldi',
                            Monosaccharide.sp2Mod:            'sp2',
                            Monosaccharide.spMod:              'sp',
                            Monosaccharide.geminalMod:    'geminal',
                            }
    
    def toMonosaccharide(self, mstr, substituents = []):
        m = Monosaccharide()
        
        MODS = mstr.split('|')
        MONO = MODS[0].split('-')

        #set mono number
        #m.set_num(int(mono_num))

        #set ring
        ring =  MONO[-1]
        ringStart,ringEnd = ring.split(':')
        m.set_ring_start(int(ringStart))
        m.set_ring_end(int(ringEnd))
        MONO.remove(ring)

        #set class
        superclass = MONO[-1]
        MONO.remove(superclass)
        m.set_superclass(superclass)

        #set anomer
        anomer = MONO[0]
        MONO.remove(anomer)
        m.set_anomer(anomer)

        #set stem
        configs = []
        stems = []
        for st in MONO:
            configs.append(st[0])
            stems.append(st[1:])
        m.set_config(configs)
        m.set_stem(stems)

        #set mods
        mod_list = []
        mod_position_list = []
        for mod in MODS[1:]:
            num,modi = mod.split(':')
            mod_position_list.append(int(num))
            mod_list.append(modi)
        m.set_mod_positions(mod_position_list)
        m.set_mods(mod_list)

        #set substituents
        #s_num = int(mono_num) + 1
        for s in substituents:
            #print 'add_substituent in new_monosaccharide > from_str', substituent.name()
            m.add_substituent(s)#,s_num)
            #s_num += 1

        #set composition
        composition = {}
        C = Composition().unit(superclass)['C']
        H = Composition().unit(superclass)['H']
        N = Composition().unit(superclass)['N']
        O = Composition().unit(superclass)['O']

        for mod in mod_list:
            C += Composition().unit(mod)['C']
            H += Composition().unit(mod)['H']
            N += Composition().unit(mod)['N']
            O += Composition().unit(mod)['O']

        for substituent_string in substituents:
            C += Composition().unit(substituent_string)['C']
            H += Composition().unit(substituent_string)['H']
            N += Composition().unit(substituent_string)['N']
            O += Composition().unit(substituent_string)['O']

        m.set_composition(C, H, O, N)

        return m

    def toGlycoCTString(self, glycoCTmono):
        if glycoCTmono._glycoct_string != '':
            return glycoCTmono._glycoct_string
        else:
            mono_string = ''
            mono_string += GlycoCTMonoFormatter._glycoCTAnomerMapping[glycoCTmono.anomer()]

            if ((glycoCTmono.config() != []) and (glycoCTmono.stem() != [])):
                for c,s in zip(glycoCTmono.config(),glycoCTmono.stem()):
                    mono_string = mono_string  + '-' + GlycoCTMonoFormatter._glycoCTConfigMapping[c] + GlycoCTMonoFormatter._glycoCTStemMapping[s]
                
            mono_string = mono_string + '-' + GlycoCTMonoFormatter._glycoCTClassMapping[glycoCTmono.superclass()] + '-' + str(glycoCTmono.ring_start()) + ':' + str(glycoCTmono.ring_end())
            
            mod_index = 0
            for mod,position in zip(glycoCTmono.mods(),glycoCTmono.mod_positions()):
                mono_string = mono_string + '|' + str(position) + ':' + GlycoCTMonoFormatter._glycoCTModMapping[mod]
            return mono_string

class LinearCodeMonoFormatter(MonoFormatter):
    #links
    link_forNeuAc = Linkage(child=None,plink_pos=5,clink_pos=1,plink_type='d',clink_type='n',pnum=None,cnum=None,num=None)
    link_forHexNAc = Linkage(child=None,plink_pos=2,clink_pos=1,plink_type='d',clink_type='n',pnum=None,cnum=None,num=None)

    #substituents
    nacetyl_forNeuAc = Substituent(sstr='n-acetyl',num=None,link=copy.copy(link_forNeuAc))#link=copy.copy(link_forNeuAc))
    nacetyl_forHexNAc = Substituent(sstr='n-acetyl',num=None,link=copy.copy(link_forHexNAc))#link=copy.copy(link_forHexNAc))

    GN = Monosaccharide().setAll(-1,'glc',8,15,6,1,'d','x',1,5,'HEX',None,'GlcNAc','GN',3,[copy.copy(nacetyl_forHexNAc)])
    AN = Monosaccharide().setAll(-1,'gal',8,15,6,1,'d','x',1,5,'HEX',None,'GalNAc','AN',4,[copy.copy(nacetyl_forHexNAc)])
    NN = Monosaccharide().setAll(-1,'gro,gal',11,19,9,1,'d,d','a',2,6,'NON','1a,2keto,3d','NN',6,[copy.copy(nacetyl_forNeuAc)])
    M  = Monosaccharide().setAll(-1,'man',6,12,6,0,'d','x',1,5,'HEX',None,'Man','M',5,[])
    A  = Monosaccharide().setAll(-1,'gal',6,12,6,0,'d','x',1,5,'HEX',None,'Gal','A',2,[])
    F  = Monosaccharide().setAll(-1,'gal',6,12,5,0,'l','x',1,5,'HEX','6d','Gal','A',2,[])
    G  = Monosaccharide().setAll(-1,'glc',6,12,6,0,'d','x',1,5,'HEX',None,'Glc','G',3,[])

    _linearMapping = {'GN': GN,
                      'AN': AN,
                       'M': M,
                       'A': A,
                       'F': F,
                      'NN': NN,
                       'G': G,
                      }

    def toMonosaccharide(self, mstr, substituents = []):
        new_mono = copy.copy(LinearCodeMonoFormatter()._linearMapping[mstr]) #COPY
        print 'new_mono in toMonosaccharide', new_mono
        if substituents != []:
            for s in substituents:
                if (((new_mono.stem() == ['glc']) or (new_mono.stem() == ['gal'])) and (s == 'n-acetyl')):
                    print 'test 1 in substituent addition in LinearCodeMonoFormatter'
                    new_mono = LinearCodeMonoFormatter().addSubstituent(new_mono, s.copy())
                elif ((new_mono.stem() == ['gro,gal']) and (s == 'n-acetyl')):
                    print 'test 2 in substituent addition in LinearCodeMonoFormatter'
                    new_mono = LinearCodeMonoFormatter().addSubstituent(new_mono, s.copy())
                else:
                    print 'test 3 in LinearCodeMonoFormat, WARNING: substituent not added'
        return new_mono

    def toLinearCode(self, m):
        if m._lincode in self._linearMapping:
            return m._lincode
        else:
            return '?'

class GlycanFormatter:
    pass
    #raise "Not implemented"

class LinearGlyFormatter(GlycanFormatter):

    mono_symbols = ['GN', 'M', 'A', 'F', 'NN','AN']
    _linearMapper = LinearCodeMonoFormatter._linearMapping

    reverse_dict = {')':'(',
                    '(':')',
                    'NG':'GN',
                    'NA':'AN',
                    'M':'M',
                    'A':'A',
                    'F':'F',
                    'NN':'NN',
                    'JN':'NJ',#not in dict, but in lots of glycans in Cartoonist
                    '4b':'b4',
                    '2b':'b2',
                    '3a':'a3',
                    '6a':'a6',
                    '??':'??',
                    }

    def toLinearString(self,glycan):
        print 'add code'

    def toGlycan(self, linear_string):
        new_branch = False
        temp = ''
        temp_memory = None
        rev_str = ''
        
        g = Glycan()
        rev_str = self.reverseString(linear_string)
        branch_at = []
        parent = None

        #FOR TEST PURPOSES ONLY
        check_next_res = False
        unknown = False
        #
        skip = False
        temp2 = ''
        for n,c in enumerate(rev_str):
            if skip == False:
                temp+=c
                if temp in self.mono_symbols:
                    if temp == 'A':
                        try:
                            temp2 = temp + rev_str[n+1] ####This is a not-so-good fix that enables this loop to catch "AN" and "A"
                        except:
                            pass
                        if temp2 in self.mono_symbols:
                            #print "FOUND!!!!!!!!!!!!!!!!!!!!!",temp2
                            temp = temp2
                            skip = True
                    m = copy.deepcopy(Monosaccharide().copy(self._linearMapper[temp]))
                    #m = copy.copy(self._linearMapper[temp])#####COPY                    
                    #m = self._MonosaccharideDictionary.fromLinearCode(temp)
                    g.add_mono(m,parent)
                    if parent != None:
                        parent = None
                    temp_memory = m
                    temp = ''

                if temp == 'NJ':
                    unknown = True

                if temp.startswith('a') or temp.startswith('b') or temp.startswith('?'):
                    if len(temp) == 2:
                        #print temp, 'in links'
                        saved_link = temp
                        #print temp
                        if new_branch == True:
                            #need to incorporate link info
                            #print 'add/fix',temp
                            new_branch = False
                        temp = ''
                        
                if temp == '(':
                    #new branch
                    #print temp, 'branch point!'
                    branch_at.append(temp_memory)
                    temp = ''
                elif temp == ')':
                    #end of branch
                    #print temp, 'branch point!'
                    new_branch = True
                    parent = branch_at.pop()
                    #branch_at.append(g)
                    temp = ''

                last_character = temp
            else:
                skip = False

        if last_character != ')':
            branch_at.append(g)

        if unknown == True:
            g = None

        #number residues and links
        number_all_members(g.root(),1,True)
        number_all_links(g.root(),1,True)

        return g

    def reverseString(self,linear_string):
        linear_string = reversed(linear_string)
        new_str = ''

        #Reverse string
        reverse_temp = ''
        for c in linear_string:
            reverse_temp += c
            if reverse_temp in self.reverse_dict:
                new_str += self.reverse_dict[reverse_temp]
                reverse_temp = ''

        return new_str

class GlycoCTGlyFormatter(GlycanFormatter):
    
    def toGlycoCTString(self,glycan):
        RES_LIST =  []
        res_num = 1
        LIN_LIST = []
        l_count = 1

        count = 1
        for ms in glycan.all_members():
            count += 1
            RES_LIST.append(ms.glycoct_line())

        link_count = 1
        for l in sorted(glycan.all_links(),key=lambda link: link.cnum()):
            l.set_num(link_count)
            link_count += 1
            LIN_LIST.append(l.__str__())

        glycan_string = 'RES'
        for m in RES_LIST:
            glycan_string += ('\n'+ m)

        glycan_string += '\nLIN'
        for l in LIN_LIST:
            glycan_string += ('\n'+l)

        return glycan_string

    def toGlycan(self, glycan_string):
        G = Glycan()
        memory = ''
        R_ = {}
        L_list = []

        #FIRST PASS - SEPARATE RES AND LIN, AND CREATE A DICTIONARY OF RES INDEXED BY RES_NUM
        for l in glycan_string.split('\n'):
            if l.strip() == 'RES':
                memory = l
            elif l.strip() == 'LIN':
                memory = l
            elif l.strip() == '':
                pass
            else:
                if memory == 'RES':
                    num,ct_str = l.split(':',1)
                    n = num[:-1]
                    s_or_b = num[-1]
                    R_[int(n)] = (s_or_b,ct_str)
                elif memory == 'LIN':
                    L_list.append(l)

        #PARSE LINKS
        MONOS = {}
        p_num = p_type = p_link_pos = c_link_pos = c_num = c_type = ''

        MONOS[1] = Monosaccharide().from_str(1,R_[1][1])
        G = Glycan()

        for l in L_list:
            #PARSE LINK STRING
            link_num, link = l.split(':',1)
            parent,left_over = link.split('(')
            link_pos,child = left_over.split(')')

            p_type = parent[-1]
            p_num = int(parent[:-1])
            p_link_pos = int(link_pos.split('+')[0])
            c_link_pos = int(link_pos.split('+')[1])
            c_type = child[-1]
            c_num = int(child[:-1])

            l = (int(link_num),p_num,p_type,p_link_pos,c_link_pos,c_num,c_type)

            #MAKE THE CHILD OBJECT (MONO OR SUBSTITUENT)
            if c_type == 'n': #Substituent

                #Create linkage
                lin = Linkage()
                lin.set_link(child=None,plink_position=p_link_pos,clink_position=c_link_pos,plink_type=p_type,clink_type=c_type,pnum=p_num,cnum=int(c_num),num=link_num)

                #Create substituent and add the linkage to the substituent
                s = Substituent()
                s.set_substituent(R_[int(c_num)][1],n=c_num,link=lin)#########n=-1 before
                s.set_link(lin)
                
                #add substituent to parent monosaccharide
                MONOS[int(p_num)].add_substituent(s)                #add_substituent(self,s,num=None)
               
            elif c_type == 'd': #Monosaccharide

                #Create Monosaccharide
                m = Monosaccharide().from_str(mono_num=int(c_num),mono_string=R_[int(c_num)][1],substituent=None)

                #Add monosaccharide to the MONOS dictionary
                MONOS[int(c_num)] = m

        #for n in sorted(MONOS):
        #    print n, MONOS[n]

        G.add_mono(MONOS[1],None)
        for l in L_list:

            #PARSE LINK STRING
            link_num, link = l.split(':',1)
            parent,left_over = link.split('(')
            link_pos,child = left_over.split(')')

            p_type = parent[-1]
            p_num = int(parent[:-1])
            p_link_pos = int(link_pos.split('+')[0])
            c_link_pos = int(link_pos.split('+')[1])
            c_type = child[-1]
            c_num = int(child[:-1])

            if c_type == 'd': #Monosaccharide
                #(child=m,plink_position=p_link_pos,clink_position=c_link_pos,plink_type=p_type,clink_type=c_type,pnum=p_num,num=link_num)
                G.add_mono(MONOS[int(c_num)],MONOS[int(p_num)],p_link_pos,c_link_pos,p_type,c_type,p_num,c_num,link_num)#int(c_num))

        return G

if __name__ == '__main__':

    #TEST GlycoCTMonoFormatter()
    monofmt = GlycoCTMonoFormatter()
    m = monofmt.toMonosaccharide('a-dgro-dgal-NON-2:6|1:a|2:keto|3:d',['n-acetyl'])
    print 'test toGlycoCTString',monofmt.toGlycoCTString(m)

    #TEST LinearCodeMonoFormatter
    monofmt2 = LinearCodeMonoFormatter()
    m2 = monofmt2.toMonosaccharide('GN')
    print '\ntest toLinearCode', monofmt2.toLinearCode(m2)

    print 'finished mono test\n'

    ct_string = """RES
1b:b-dglc-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4s:n-acetyl
5b:b-dman-HEX-1:5
6b:a-dman-HEX-1:5
7b:x-dglc-HEX-1:5
8b:x-dgal-HEX-1:5
9s:n-acetyl
10b:a-dman-HEX-1:5
11b:x-dglc-HEX-1:5
12b:x-dgal-HEX-1:5
13b:x-dgro-dgal-NON-2:6|1:a|2:keto|3:d
14s:n-acetyl
15s:n-acetyl
LIN
1:1d(2+1)2n
2:1o(4+1)3d
3:3d(2+1)4n
4:3o(4+1)5d
5:5o(3+1)6d
6:6o(-1+1)7d
7:7o(-1+1)8d
8:7d(2+1)9n
9:5o(6+1)10d
10:10o(-1+1)11d
11:11o(-1+1)12d
12:12o(-1+2)13d
13:13d(5+1)14n
14:11d(2+1)15n
"""

    glyfmt = GlycoCTGlyFormatter()
    new_g = glyfmt.toGlycan(ct_string)
    print 'test',new_g
    newg_2 = glyfmt.toGlycoCTString(new_g)
    print newg_2


    glyfmt2 = LinearGlyFormatter()
    example1 = "Ab4GNb2Ma3(Ma6)Mb4GNb4GN"
    example2 = "NNa6Ab4GNb2(NNa3Ab4GNb4)Ma3(NNa6Ab4GNb2Ma6)Mb4GNb4(Fa6)GN"
    g = glyfmt2.toGlycan(example1)

    

    print example1
    print g.print_g()
    print g.mass()

    g2 = glyfmt.toGlycoCTString(g)
    print g2

    for l in open(file):
        g = GlycanReader(file,glyfmt)






    
