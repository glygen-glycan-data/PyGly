
from __future__ import print_function

import re
import os
from pkg_resources import resource_stream
try:
    import configparser as ConfigParser
except ImportError:
    import ConfigParser

try:
    from StringIO import StringIO
except ImportError:
    pass

from . Monosaccharide import *

from . GlycanFormatterExceptions import WURCS20ParseError

class InvalidMonoError(WURCS20ParseError):
    def __init__(self, monostr):
        self.message = "WURCS2.0 parser: Invalid monosaccharide: %s" % (monostr,)

    def __str__(self):
        return self.message


class UnsupportedMonoError(WURCS20ParseError):
    def __init__(self, monostr):
        self.message = "WURCS2.0 parser: Unsupported monosaccharide: %s" % (monostr,)

    def __str__(self):
        return self.message


class UnsupportedSkeletonCodeError(UnsupportedMonoError):
    def __init__(self, skel):
        self.message = "WURCS2.0 parser: Unsupported skeleton code: %s" % (skel,)


class UnsupportedSubstituentError(UnsupportedMonoError):
    def __init__(self, sub):
        self.message = "WURCS2.0 parser: Unsupported substituent: %s" % (sub,)

def readconfig(inifilename):
    iniFile = [ s.decode('utf8') for s in resource_stream(__name__, inifilename).read().splitlines() ]
    cfg = ConfigParser.ConfigParser()
    if hasattr(cfg,'read_file'):
        cfg.read_file(iniFile,inifilename)
    else:                                                                                                                
        cfg.readfp(StringIO(u'\n'.join(iniFile)),inifilename)
    return cfg

class WURCS20MonoFormat:
    mono_pattern = re.compile(
        r"^([0-9a-zA-Z]{3,9})((-\d[abx])(_[0-9?]-[0-9?])?)?((_([0-9?]|[0-9]\|[0-9]|[0-9]-[0-9])(\*[^_]+)?)*)$")

    def __init__(self):
        self.cache = {}
        self.load()

    def load(self):
        self.skelconfig = readconfig('wurcs20_skeleton.ini')
        self.subsconfig = readconfig('wurcs20_substituent.ini')

    def getsubst(self,sub_name):
        try:
            sub_type = self.subsconfig.get(sub_name, "type")
            subst = Substituent(eval(sub_type))
        except (ConfigParser.NoSectionError, ValueError):
            raise UnsupportedSubstituentError(sub_name)
        return subst

    skel_config_get_default_section = ""

    def skel_config_get(self, optionx, skelton_code=skel_config_get_default_section):
        try:
            return self.skelconfig.get(self.skel_config_get_default_section, optionx)
        except ConfigParser.NoOptionError:
            return None

    def parsing(self, mono_string):
        parsed = re.search(self.mono_pattern, mono_string)

        # print(mono_string)
        if not parsed:
            # Bad formatted WURCS mono
            if mono_string != "<Q>":
                raise UnsupportedMonoError(mono_string)
            skel = mono_string
            anomer = ""
            ring = ""
            substs = ""
        else: 
            skel = parsed.group(1)
            anomer = parsed.group(3)
            ring = parsed.group(4)
            substs = parsed.group(5)

        # print(mono_string)
        # print(skel, anomer, ring, substs)
        m = Monosaccharide()
        m.set_external_descriptor(mono_string)

        if skel in self.skelconfig.sections():
            self.skel_config_get_default_section = skel

            anomer_s = self.skel_config_get("anomer")
            if anomer_s:
                m.set_anomer(eval(anomer_s))
                if m.anomer() == Anomer.uncyclized:
                    m.set_ring_start(0)
                    m.set_ring_end(0)

            stem = self.skel_config_get("stem")
            if stem:
                m.set_stem(*(map(lambda s: eval(s), stem.split())))

            config = self.skel_config_get("config")
            if config:
                m.set_config(*(map(lambda c: eval(c), config.split())))
            else:
                if stem:
                    m.set_config(*([None] * len(m.stem())))

            if m.stem() and m.config():
                assert len(m.stem()) == len(m.config()), "Inconsistent config/stem for %s" % (skel,)
            elif not m.stem():
                assert not m.config(), "Problem with config/stem for %s"%(skel, )

            superclass = self.skel_config_get("superclass")
            if superclass:
                m.set_superclass(eval(superclass))

            mods = self.skel_config_get("mods")
            if mods:
                mods_list = mods.split()
                for i in range(0, len(mods_list), 2):
                    m.add_mod(mods_list[i], eval(mods_list[i + 1]))

        else:
            # skeleton code not supported
            raise UnsupportedSkeletonCodeError(skel)

        if anomer:
            match = re.search(r"^-(\d)(.)$", anomer)
            if not match:
                raise InvalidMonoError(mono_string)
            if anomer_s and eval(anomer_s) != Anomer.missing:
                raise InvalidMonoError(mono_string)
            if match.group(2) == "a":
                m.set_anomer(Anomer.alpha)
            elif match.group(2) == "b":
                m.set_anomer(Anomer.beta)
            elif match.group(2) == "x":
                m.set_anomer(Anomer.missing)
            else:
                raise InvalidMonoError(mono_string)

        if ring:
            match = re.search(r"^_(\d|\?)-(\d|\?)$", ring)
            if not match:
                raise InvalidMonoError(mono_string)
            try:
                r_s = int(match.group(1))
                m.set_ring_start(r_s)
            except ValueError:
                pass
            try:
                r_e = int(match.group(2))
                m.set_ring_end(r_e)
            except ValueError:
                pass

        if substs:
            sub_list = substs[1:].split("_")
            for sub in sub_list:
                try:
                    pp, sub_name = sub.split("*", 1)
                except ValueError:
                    # missing *, no sub_name, coresponds to pp =
                    pp = sub
                    sub_name = "anhydro"
                    # Note, we appear to require \d-\d for sub, otherwise it is a no-op?
                    # cases with \d only do not appear to change the mass...
                    if not re.search(r'^\d-\d$', sub):
                        raise InvalidMonoError(mono_string)
                try:
                    sub_type = self.subsconfig.get(sub_name, "type")
                    pt = self.subsconfig.get(sub_name, "parent_type")
                    ct = self.subsconfig.get(sub_name, "child_type")
                    cp = self.subsconfig.get(sub_name, "child_pos")
                    sub_object = Substituent(eval(sub_type))
                    sub_object.set_external_descriptor(sub_name)

                    if cp != None:
                        cp = int(cp)

                    if re.search(r"^\d-\d$", pp) and sub_name != "anhydro":
                        for pp in pp.split("-"):
                            pp = int(pp)
                            m.add_substituent(sub_object, parent_pos=pp, parent_type=eval(pt), child_pos=cp,
                                              child_type=eval(ct))

                    elif '|' not in pp and '-' not in pp:
                        if pp == "?" or pp == "-1":
                            pp = -1
                        else:
                            pp = int(pp)
                        m.add_substituent(sub_object, parent_pos=pp, parent_type=eval(pt), child_pos=cp,
                                          child_type=eval(ct))
                    elif '|' in pp:
                        pp = list(map(int, pp.split('|')))
                        m.add_substituent(sub_object, parent_pos=pp, parent_type=eval(pt), child_pos=cp,
                                          child_type=eval(ct))
                    elif '-' in pp and sub_name == "anhydro":
                        pp = list(map(int, pp.split('-')))
                        if len(pp) != 2:
                            raise UnsupportedSubstituentError(sub)
                        m.add_substituent(sub_object, parent_pos=pp[0], parent_type=eval(pt),
                                                  child_pos=cp, child_type=eval(ct)).child()
                        pt = self.subsconfig.get(sub_name + "1", "parent_type")
                        ct = self.subsconfig.get(sub_name + "1", "child_type")
                        cp = self.subsconfig.get(sub_name + "1", "child_pos")
                        if cp != None:
                            cp = int(cp)
                        m.add_substituent(sub_object, parent_pos=pp[1], parent_type=eval(pt), child_pos=cp,
                                          child_type=eval(ct))
                    else:
                        raise UnsupportedSubstituentError(sub)
                except (ConfigParser.NoSectionError, ValueError):
                    # Substituent not supported or badly formatted
                    raise UnsupportedSubstituentError(sub)

        return m

    def get(self, mono_string):
        try:
            result = self.cache[mono_string]
        except KeyError:
            result = self.parsing(mono_string)
            self.cache[mono_string] = result

        return result.clone()


if __name__ == "__main__":

    import sys

    mf = WURCS20MonoFormat()
    args = False
    for w in sys.argv[1:]:
        args = True
        print(w, file=sys.stderr)
        print(mf.get(w))
    if not args:
        for l in sys.stdin:
            print(l.strip(), file=sys.stderr)
            print(mf.get(l.strip()))
