s = str(sys.argv[1:][0])
name = 'AF-'+s+'-F1-model_v1.pdb'
from pymol import cmd
cmd.load(name)
show surface
cmd.save(s+'_mesh.obj')
quit
