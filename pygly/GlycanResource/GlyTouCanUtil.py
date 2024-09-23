
try:
    from pygly.GlycanFormatter import WURCS20Format, GlycoCTFormat, \
                                      GlycanParseError, ZeroPlusLinkCountError, \
                                      UndeterminedLinkCountError, CircularError, \
                                      LinkCountError

    from pygly.WURCS20MonoFormatter import WURCS20MonoFormat, \
                                           UnsupportedSkeletonCodeError, \
                                           UnsupportedSubstituentError, \
                                           InvalidMonoError
except ImportError:
    from .. GlycanFormatter import WURCS20Format, GlycoCTFormat, \
                                GlycanParseError, ZeroPlusLinkCountError, \
                                UndeterminedLinkCountError, CircularError, \
                                LinkCountError

    from .. WURCS20MonoFormatter import WURCS20MonoFormat, \
                                     UnsupportedSkeletonCodeError, \
                                     UnsupportedSubstituentError, \
                                     InvalidMonoError

from collections import defaultdict

class GlyTouCanUtil(object):
    _wurcs_mono_format = WURCS20MonoFormat()
    _wurcs_format = WURCS20Format()
    _glycoct_format = GlycoCTFormat()
    _alphamap = None

    def getUnsupportedCodes(self, acc):
        codes = set()
        substs = set()
        invalid = set()
        other = set()
        sequence = self.getseq(acc, 'wurcs')
        if not sequence:
            return codes, substs, invalid, other
        monos = sequence.split('/[', 1)[1].split(']/')[0].split('][')
        for m in monos:
            try:
                g = self._wurcs_mono_format.parsing(m)
            except UnsupportedSkeletonCodeError as e:
                codes.add(e.message.rsplit(None, 1)[-1])
            except UnsupportedSubstituentError as e:
                substs.add(e.message.rsplit(None, 1)[-1])
            except InvalidMonoError as e:
                invalid.add(e.message.rsplit(None, 1)[-1])
            except GlycanParseError:
                pass
        try:
            g = self._wurcs_format.toGlycan(sequence)
        except ZeroPlusLinkCountError:
            other.add("0+ link count")
        except UndeterminedLinkCountError:
            other.add("undetermined link count")
        except CircularError:
            other.add("circular")
        except LinkCountError:
            other.add("bad link count")
        except GlycanParseError:
            pass
        return codes, substs, invalid, other

    def getGlycan(self, acc, format=None):
        if not format or (format == 'wurcs'):
            sequence = self.getseq(acc, 'wurcs')
            if sequence:
                try:
                    return self._wurcs_format.toGlycan(sequence)
                except GlycanParseError:
                    pass  # traceback.print_exc()
        if not format or (format == 'glycoct'):
            sequence = self.getseq(acc, 'glycoct')
            if sequence:
                try:
                    return self._glycoct_format.toGlycan(sequence)
                except GlycanParseError:
                    pass
        return None

    def glycoct(self, acc, fetch=None):
        g = self.getGlycan(acc,fetch)
        if not g:
            return None
        try:
            return g.glycoct()
        except (LookupError,ValueError):
            pass
        except RuntimeError:
            pass
        return None

    def umw(self, accession=None, format=None):
        g = self.getGlycan(accession,format)                                                                             
        if not g:
            return None
        try:
            return g.underivitized_molecular_weight()
        except (LookupError,ValueError):
            pass
        return None

    def wurcs2glycoct(self, acc):
        sequence = self.getseq(acc,'wurcs')
        if sequence:
            sequence1 = urllib.parse.quote_plus(sequence)
            url = 'https://api.glycosmos.org/glycanformatconverter/2.3.2-snapshot/wurcs2glycoct/'+sequence1
            try:
                data = json.loads(urllib.request.urlopen(url).read())
                if 'GlycoCT' in data:
                    return data['GlycoCT']
            except ValueError:
                pass
        return None

    def subsumptionbyapi(self, acc):
        sequence = self.getseq(acc,'wurcs')
        if sequence:
            sequence1 = urllib.parse.quote_plus(sequence)
            url = 'https://api.glycosmos.org/subsumption/0.2.0/'+sequence1
            data = urllib.request.urlopen(url).read()
            seen = set()
            lasts = None
            for triple in sorted([tuple([s.strip() for s in list(map(str,list(map(t.get,("S","P","O")))))]) for t in json.loads(data)]):
                if triple in seen:
                    continue
                seen.add(triple)
                if triple[0] != lasts:
                    if lasts != None:
                        print("")
                    print(triple[0])
                    lasts = triple[0]
                if triple[2] == sequence:
                    print(">>  "+"\t".join(triple[1:]))
                else:
                    print("    "+"\t".join(triple[1:]))

    def findskel(self, skel, maxcount=None):
        if maxcount != None:
            maxcount = int(maxcount)
        
        for acc, format, wurcs in self.allseq(format='wurcs'):
            glycoct = self.getseq(acc,format='glycoct')
            if not glycoct:
                continue
            monos = wurcs.split('/[', 1)[1].split(']/')[0].split('][')
            if maxcount != None and len(monos) > maxcount:
                continue
            for mono in monos:
                msk = re.search(r'^(.*?)([-_].*)?$',mono).group(1)
                assert msk
                m = re.search(r"^%s$"%(skel,),msk)
                if m:
                    yield acc, m.group(0)

    def multiseq(self):
        counts = defaultdict(set)
        for acc,fmt,seq in self.allseq():
            counts[(acc,fmt)].add(seq)
        for k,v in list(counts.items()):
            if len(v) > 1:
                yield k

    symbol2wurcs_definition = """
	NeuAc   AUd21122h_5*NCC/3=O   1
	NeuGc   AUd21122h_5*NCCO/3=O  2
	Fuc     u1221m                6
	Hex     uxxxxh                5
	HexNAc  uxxxxh_2*NCC/3=O      4
	dHex    uxxxxm                3
        HexA    uxxxxA                7
        HexN    uxxxxh_2*N            8
	Pent    uxxxh                 10
	P       *OPO/3O/3=O           -1
	Phospho *OPO/3O/3=O           -1
	S       *OSO/3=O/3=O          -1
	Sulpho  *OSO/3=O/3=O          -1
	"""
    symbol2wurcs=None
    wurcsorder=None

    def makecompseq(self,**incomp):
        if self.symbol2wurcs == None:
            self.symbol2wurcs = {}
            self.wurcsorder = {}
            for l in self.symbol2wurcs_definition.splitlines():
                if not l.strip():
                    continue
                sl = l.split()
                self.symbol2wurcs[sl[0]] = sl[1]
                self.wurcsorder[sl[1]] = int(sl[2])
        comp = {}
        subst = {}
        for k,v in incomp.items():
            if k not in self.symbol2wurcs:
                raise ValueError("Bad symbol in composition")
            skel = self.symbol2wurcs[k]
            v = int(v)
            if self.wurcsorder[skel] < 0:
                subst[skel] = v
            else:
                comp[skel] = v
        skels = sorted(comp,key=self.wurcsorder.get)
        total = sum(comp.values())
        if total == 0:
            raise ValueError("No monosaccharides in composition")
        uniq = len(skels)
        wurcsseq = "WURCS=2.0/%s,%s,%s/" % (uniq, total, "0+")
        wurcsseq += "".join(map(lambda sk: "[%s]" % (sk,), skels)) + "/"
        inds = []
        for i, sk in enumerate(skels):
            inds.extend([str(i + 1)] * comp[sk])
        wurcsseq += "-".join(inds)
        wurcsseq += "/"
        return self.fixcompwurcs(wurcsseq,subst)

    def fixcompwurcs(self, wurcsseq, subst=[]):
        if not self._alphamap:
            self._alphamap = dict()
            for i, c in enumerate(range(ord('a'), ord('z') + 1)):
                self._alphamap[i + 1] = chr(c)
                self._alphamap[chr(c)] = (i + 1)
            for i, c in enumerate(range(ord('A'), ord('Z') + 1)):
                self._alphamap[i + 1 + 26] = chr(c)
                self._alphamap[chr(c)] = (i + 1 + 26)
        prefix, counts, rest = wurcsseq.split('/', 2)
        unodes, nodes, edges = counts.split(',')
        nodes = int(nodes)
        assert '0+' in edges
        edges = (nodes - 1)
        ambignode = "|".join(["%s?" % (self._alphamap[i],) for i in range(1, nodes + 1)])
        ambigedge = "%s}-{%s" % (ambignode, ambignode)
        ambigedges = [ambigedge] * edges
        if hasattr(subst,'items'):
            subst = list(subst.items())
        for sub,cnt in subst:
            for i in range(cnt):
                ambigedges.insert(0,"%s}%s"%(ambignode,sub))
        return "%s/%s,%d,%d/%s%s" % (prefix, unodes, nodes, len(ambigedges), rest, "_".join(ambigedges))

    def anyglycan2wurcs(self, glycan):
        sequence = ""
        if isinstance(glycan, Glycan.Glycan):
            if not self.glycoct_format:
                self.glycoct_format = GlycoCTFormat()
            sequence = self.glycoct2wurcs(self.glycoct_format.toStr(glycan))
            if '0+' in sequence:
                sequence = self.fixcompwurcs(sequence)
        else:
            sequence = re.sub(r'\n\n+', r'\n', glycan)
            if sequence.strip().startswith('RES'):
                sequence = self.glycoct2wurcs(glycan)
        return sequence

    def glycoct2wurcs(self, seq):
        requestURL = "https://api.glycosmos.org/glycanformatconverter/2.3.2-snapshot/glycoct2wurcs/"
        encodedseq = urllib.parse.quote(seq, safe='')
        requestURL += encodedseq
        req = urllib.request.Request(requestURL)
        # self.wait()
        response = urllib.request.urlopen(req).read()

        result = json.loads(response)

        try:
            wurcs = result["WURCS"]
        except:
            raise ValueError("GlycoCT 2 WURCS conversion failed")

        return wurcs.strip()
