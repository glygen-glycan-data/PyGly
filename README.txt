
Installation
============

Installation for Linux only, at this time.

1. Unpack the module wherever you intend to work on it. Do not install
   in Python's site-packages! The bin and data directories contain
   useful stuff for managing the glycan structure databases, but are
   not appropriate for site-packages.

     tar zxf PyGly-1.0.0.tgz

2. (Optional). If you have the source only version of the module and
   support scripts, you'll need to execute scripts to download and build
   and index taxonomy and GlycomeDb databases.

     cd PyGly-1.0.0
     ./bin/init.sh

3. (Optional). Install the module so that other python programs can find it, without
   needing to use PYTHONPATH or manipulate sys.path.

     cd PyGly-1.0.0
     python setup.py --user develop

4. Dump a representative, with respect to topological (branching)
   equivalence, of all N-linked glycan structures with exactly 3 Hex
   monosaccharides and molecular weight at most 1500Da (underivitized)
   from human GlycomeDB:

     cd PyGly-1.0.0
     ./bin/GlyDbIndex.sh dump  data/GlycomeDb-Human.gdb toporepr 0 nlinked 1 Hex 3 maxmw 1500
     ./bin/GlyDbIndex.py count data/GlycomeDb-Human.gdb toporepr 0 nlinked 1 Hex 3 maxmw 1500

