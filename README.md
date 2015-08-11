# endysi

Endysi (ENsemble DYnamics SImulator) is a program for creating, 
simulating, and analysing models of competetive endogenous RNA (ceRNA) 
networks.  It was developed to facilitate my PhD research in this area. 
At the moment, it is restricted to this type of network only, but does 
allow one to create networks of any size with fixed connectivity.  (In
general, a ceRNA network is a bipartite network composed of M miRNA 
nodes, N ceRNA nodes, each connected by K edges.)  In the future, I 
might extend it for other types of biochemical networks (such as, 
signal transduction networks, gene regulatory networks, and metabolic
networks).  

This project is primarily for personal use; it is published here to 
allow for easy deployment of the latest changes to various machines 
used for research.  As such, there are a few things to mention to any 
one trying to run it. 

* First, Endysi WILL NOT WORK ON WINDOWS! I run Linux on all my machines,
and didn't bother with Windows support.  It SHOULD run on any *nix based
system (e.g., OSX or *BSD), but I've only ever run it on Linux-based 
systems.

* Endysi uses BioNetGen (www.bionetgen.org) on the backend.  It expects
the downloaded BioNetGen distribution to be in ~/apps/BioNetGen. No 
doubt I should change this, but hey, so far I'm the only one using this
thing.  BioNetGen itself depends on Perl, but should be installed by
default on any Linux distro.

* Python dependencies: 
    0. Written in Python 2.7, but I have run it in Python 3 (just for kicks).
    1. matplotlib
    2. scipy
    3. numpy
    4. pexpect
    5. networkx

* I personally keep my development/running environment separate from the system
environment using pyenv and plenv for Python and Perl, respectively. But 
this is not a requirement; it simply allows me to maintain control over
the running environment without messing with system libraries. 

### Installing and running

1. To install, simply clone the repo wherever you want and make sure you 
have the above dependencies installed.  

2. cd into the directory above and run:
`python2 endysi.py --help`
This will tell you about the command line arguments.  For example:
`python2 endysi.py -m 5 -n 10 -k 2 -s 100 --method ode`
will create 100 networks with 5 miRNAs and 10 ceRNAs with each miRNA 
targeting 2 ceRNAs and each ceRNA targeted by 2 miRNAs; these networks 
will then be simulated as systems of ODEs.  The other option for 
`--method` is 'ssa' (i.e., a Gillespie simulation).

3. When it runs, the above command will create a new directory for the 
raw data and results in ~/research/results/ceRNA/endysi.  


