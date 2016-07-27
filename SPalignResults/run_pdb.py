import subprocess as sub
import sys

# tm_pdb = ('Tm',"pdb1i4n_A.ent")
# tt_pdb = ('Tt',"pdb1vc4_A.ent")
# ss_pdb = ('Ss',"pdb2c3z_A.ent")

if   sys.argv[1] == 'Tm':
    template_name, current_template = 'Tm',"pdb1i4n_A.ent"
elif sys.argv[1] == 'Tt':
    template_name, current_template = 'Tt',"pdb1vc4_A.ent"
elif sys.argv[1] == 'Ss':
    template_name, current_template = 'Ss',"pdb2c3z_A.ent"
else:
    print "blah! wrong input: Tm Tt or Ss are acceptable only!"
    print "nohup python run_pdb.py Tm(Tt,Ss) &"
    sys.exit(1)

our_pds = [s.strip() for s in sub.check_output('ls ./pdbA/ | grep .ent',shell=True).strip().split('\n')]


def get_cmd(pdb1,pdb2):
    return "./SPalignNS -pair ./pdbA/%s ./pdbA/%s"%(pdb1,pdb2)

for pdb in our_pds:
    result = sub.check_output(get_cmd(current_template,pdb),shell=True)
    fp = open(pdb.replace('.','_')+'_%s.aln'%template_name,'w')
    fp.write(result)
    fp.close()
    # with open(pdb+'.aln','w') as fp:
    #   fp.write(result)

