import re
import sys
import ConfigParser
from Monosaccharide import *

class WURCS2_MONO_PARSER():

    mono_pattern = re.compile(r"^([0-9a-zA-Z]{3,9})(-\d[abx])?(_[0-9?]-[0-9?])?(_([0-9?]|[0-9]\|[0-9])\*.*)?")

    def __init__(self):
        self.cache = {}
        self.load()

    def load(self):
        self.skelconfig = ConfigParser.SafeConfigParser()
        self.skelconfig.read('skeleton.ini')
        self.subsconfig = ConfigParser.SafeConfigParser()
        self.subsconfig.read('substituent.ini')

    skel_config_get_default_section = ""
    def skel_config_get(self,optionx,skelton_code = skel_config_get_default_section):
        try:
            return self.skelconfig.get(self.skel_config_get_default_section, optionx)
        except ConfigParser.NoOptionError:
            return None

    def parsing(self, mono_string):
        parsed = re.match(self.mono_pattern, mono_string)

        if not parsed:
            # Bad formatted WURCS mono
            raise ValueError

        skel = parsed.group(1)
        anomer = parsed.group(2)
        ring = parsed.group(3)
        substs = parsed.group(4)

        #print skel, anomer, ring, substs
        m = Monosaccharide()

        if skel in self.skelconfig.sections():
            self.skel_config_get_default_section = skel

            anomer_s = self.skel_config_get("anomer")
            if anomer_s:
                m.set_anomer(eval(anomer_s))

            config = self.skel_config_get("config")
            if config:
                if " " in config:
                    for c in config.split():
                        m.set_config(eval(c))
                else:
                    m.set_config(eval(config))

            stem = self.skel_config_get("stem")
            if stem:
                if " " in stem:
                    for s in stem.split():
                        m.set_stem(eval(s))
                else:
                    m.set_stem(eval(stem))

            superclass = self.skel_config_get("superclass")
            if superclass:
                m.set_superclass(eval(superclass))

            mods = self.skel_config_get("mods")
            if mods:
                mods_list = mods.split(" ")
                for i in range(0, len(mods_list), 2):
                    m.add_mod(mods_list[i], eval(mods_list[i + 1]))

        else:
            # skeleton code not supported
            raise ValueError

        if (not anomer_s) and anomer:
            if anomer[2] == "a":
                m.set_anomer(1)
            elif anomer[2] == "b":
                m.set_anomer(2)

        if ring and not ("?" in ring):
            r_s = int(ring[1])
            r_e = int(ring[-1])
            m.set_ring_start(r_s)
            m.set_ring_end(r_e)

        if substs:
            sub_list = substs[1:].split("_")
            for sub in sub_list:
                pp,sub_name = sub.split("*")
                try:
                    type = self.subsconfig.get(sub_name,"type")
                    pt = self.subsconfig.get(sub_name,"parent_type")
                    ct = self.subsconfig.get(sub_name,"child_type")

                    if len(pp)==1:
                        if pp=="?":
                            pp = -1
                        m.add_substituent(eval(type),parent_pos=pp,parent_type=eval(pt),child_pos=1,child_type=eval(ct))
                    else:
                        pp1=pp[0]
                        pp2=pp[-1]
                        m.add_substituent(eval(type),parent_pos=pp1,parent_pos2=pp2,parent_type=eval(pt),child_pos=1,child_type=eval(ct))

                except ConfigParser.NoSectionError:
                    # Substituent not supported now
                    print sub_name
                    # raise error

                except ConfigParser.NoOptionError:
                    # That won't happen
                    # raise error
                    pass

        return m

    def get(self, mono_string):
        try:
            result = self.cache[mono_string]

        except KeyError:
            result = self.parsing(mono_string)
            self.cache[mono_string] = result

        return result.clone()

		

if __name__ == "__main__":
	a = WURCS2_MONO_PARSER()
	try:
		print a.get("a2122h-1x_1-5_2*NCC/3=O")
	except:
		print "Not supported yet!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"

